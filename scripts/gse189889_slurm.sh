#!/usr/bin/env bash
# gse189889_slurm.sh
# GSE189889 — Acral Melanoma scRNA-seq (Li et al. Clin Cancer Res 2022)
#
#   AM1:        10x 3' v3  (8 SRR)  — soloUMIlen=12, 3M whitelist
#   AM2-AM7:    10x 5' v2  (8 SRR each)  — soloUMIlen=10, 737K whitelist
#   AM8_Node:   10x 5' v2  (4 SRR)
#   AM8_Toe:    10x 5' v2  (8 SRR)
#   Total: 10 Slurm jobs, 68 SRR runs
#
# Runs on Nixon Lab Slurm (defq: c5ad.12xlarge, 24 vCPU, 96 GB RAM)
#
# Usage:
#   bash gse189889_slurm.sh submit
#   bash gse189889_slurm.sh run <sample_id> <srr_csv> <chemistry>   # called by sbatch
#   bash gse189889_slurm.sh status
#   bash gse189889_slurm.sh sync-refs
#
# Chemistry strand note (verified 2026-06-07 against GEO soft file + L7):
#   Both 3' v3 and 5' v2 use --soloStrand Forward.
#   GEO characteristics: AM1 kit=3-prime; AM2-AM8 kit=5-prime (all confirmed).
#   The orthogonal_validation_plan.md had AM1/AM2-AM8 strands swapped — this script is correct.
#
# Disk note: 10 parallel jobs, 8 SRRs per sample ~48 GB FASTQ each → ~480 GB EFS peak during
#   download phase. EFS handles this; FASTQs are deleted per-sample immediately after STAR.

set -uo pipefail

# ── Paths ──────────────────────────────────────────────────────────────────────
EFS="/efs/users/jmarston"
WORKDIR="${EFS}/acral-melanoma"
STAR_INDEX="${EFS}/atherosclerosis/reference/star_index_hybrid"
WL_3M="${EFS}/atherosclerosis/reference/3M-february-2018.txt"      # 3' v3 (AM1)
WL_737K="${EFS}/atherosclerosis/reference/737K-august-2016.txt"    # 5' v2 (AM2-AM8)
GTF="${EFS}/atherosclerosis/reference/loci.hg38.v0_2.gtf"
FRAG_MAP="${EFS}/atherosclerosis/reference/fragment_to_tu.hg38.v0_2.tsv"
VIRAL_WL="${EFS}/atherosclerosis/reference/human_tropic_viruses.tsv"

CONDA="${EFS}/bin/miniforge3"
STAR_BIN="${CONDA}/envs/STAR/bin/STAR"
SAMTOOLS="${CONDA}/envs/riboseq/bin/samtools"    # L43b: samtools in riboseq env, not STAR env
FASTERQ="${EFS}/bin/local/fasterq-dump"          # L31: CentOS-native binary
CTQ_PY="${EFS}/ctq/.venv/bin/python"             # L37: dedicated ctq venv (numpy<2, torch<2.3)
CTQ_SCRIPT="${EFS}/ctq/scripts/run_sctequant.py"
CTQ_BIN="${EFS}/ctq/.venv/bin/ctq"

FASTQ_DIR="${WORKDIR}/fastq"
TMP_DIR="${WORKDIR}/fastq/.tmp"
STAR_DIR="${WORKDIR}/starsolo"
CTQ_DIR="${WORKDIR}/sctequant"
VIRAL_DIR="${WORKDIR}/viral"
LOG_DIR="${WORKDIR}/logs"
LOCK_DIR="${WORKDIR}/locks"

S3_OUT="s3://jez-research-data/acral-melanoma"

STAR_THREADS=14    # leaves 2 vCPU headroom on 16-vCPU defq node
N_CTQ_MAX=3

# ── Logging ────────────────────────────────────────────────────────────────────
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [${SAMPLE_ID:-MAIN}] $*" | tee -a "${LOG_DIR}/${SAMPLE_ID:-main}.log"; }

# ── Semaphore (prevent >N concurrent ctq jobs on shared EFS) ──────────────────
sem_acquire() {
    local name=$1 max=$2
    while true; do
        count=$(ls "${LOCK_DIR}/${name}."*.lck 2>/dev/null | wc -l)
        [ "${count}" -lt "${max}" ] && touch "${LOCK_DIR}/${name}.$$.lck" && return 0
        sleep 10
    done
}
sem_release() { local name=$1; rm -f "${LOCK_DIR}/${name}.$$.lck" 2>/dev/null || true; }

# ── Download one SRR and detect R1/R2 by read length (L5d) ───────────────────
# Outputs: sets globals DETECTED_R1 and DETECTED_R2 (gz paths in FASTQ_DIR)
download_srr() {
    local srr=$1 sample=$2
    DETECTED_R1="" DETECTED_R2=""

    local out_r1="${FASTQ_DIR}/${sample}_${srr}_R1.fastq.gz"
    local out_r2="${FASTQ_DIR}/${sample}_${srr}_R2.fastq.gz"

    if [[ -f "${out_r1}" && -f "${out_r2}" ]]; then
        log "  ${srr}: already present — skip download"
        DETECTED_R1="${out_r1}"
        DETECTED_R2="${out_r2}"
        return 0
    fi

    log "  ${srr}: fasterq-dump start"
    mkdir -p "${TMP_DIR}"
    # --include-technical required: SRA marks CB+UMI (R1) as "technical" — L5b
    "${FASTERQ}" "${srr}" \
        --split-files --include-technical \
        --outdir "${TMP_DIR}" --temp "${TMP_DIR}" \
        --threads 8 --progress \
        2>>"${LOG_DIR}/${sample}.log"

    local tmp_r1="" tmp_r2="" _n _f _len

    if [[ -f "${TMP_DIR}/${srr}_4.fastq" ]]; then
        # 4-file layout: detect R1/R2 by read length, not by position (L5d)
        for _n in 1 2 3 4; do
            _f="${TMP_DIR}/${srr}_${_n}.fastq"
            [[ -f "$_f" ]] || continue
            _len=$(awk 'NR==2{print length($0); exit}' "$_f")
            if [[ $_len -ge 16 && $_len -le 40 && -z "${tmp_r1}" ]]; then
                tmp_r1="$_f"
                log "    4-file: _${_n}.fastq = R1 CB+UMI (${_len}bp)"
            elif [[ $_len -gt 50 && -z "${tmp_r2}" ]]; then
                tmp_r2="$_f"
                log "    4-file: _${_n}.fastq = R2 cDNA (${_len}bp)"
            else
                log "    4-file: _${_n}.fastq = index/other (${_len}bp) — deleting"
                rm -f "$_f"
            fi
        done
        if [[ -z "${tmp_r1}" || -z "${tmp_r2}" ]]; then
            log "ERROR: 4-file layout: could not identify R1/R2 by read length for ${srr}"
            return 1
        fi
    elif [[ -f "${TMP_DIR}/${srr}_3.fastq" ]]; then
        log "    3-file layout: _2=R1 (CB+UMI), _3=R2 (cDNA)"
        tmp_r1="${TMP_DIR}/${srr}_2.fastq"
        tmp_r2="${TMP_DIR}/${srr}_3.fastq"
        rm -f "${TMP_DIR}/${srr}_1.fastq"
    else
        log "    2-file layout: _1=R1 (CB+UMI), _2=R2 (cDNA)"
        tmp_r1="${TMP_DIR}/${srr}_1.fastq"
        tmp_r2="${TMP_DIR}/${srr}_2.fastq"
    fi

    log "    compressing ${srr} R1 (pigz)"
    pigz -p 4 -c "${tmp_r1}" > "${out_r1}"
    log "    compressing ${srr} R2 (pigz)"
    pigz -p 4 -c "${tmp_r2}" > "${out_r2}"
    rm -f "${tmp_r1}" "${tmp_r2}" 2>/dev/null || true

    DETECTED_R1="${out_r1}"
    DETECTED_R2="${out_r2}"
    log "  ${srr}: done — R1 $(du -sh ${out_r1}|cut -f1)  R2 $(du -sh ${out_r2}|cut -f1)"
}

# ── Per-sample pipeline ────────────────────────────────────────────────────────
run_sample() {
    local SAMPLE_ID=$1
    local SRR_CSV=$2     # comma-separated SRR list, e.g. SRR17075055,SRR17075056,...
    local CHEMISTRY=$3   # "3prime" or "5prime"

    local STAR_OUT="${STAR_DIR}/${SAMPLE_ID}"
    local BAM="${STAR_OUT}/${SAMPLE_ID}_Aligned.sortedByCoord.out.bam"
    local BARCODES_FILTERED="${STAR_OUT}/${SAMPLE_ID}_Solo.out/GeneFull_Ex50pAS/filtered/barcodes.tsv"
    local BARCODES_RAW="${STAR_OUT}/${SAMPLE_ID}_Solo.out/GeneFull_Ex50pAS/raw/barcodes.tsv"
    local BARCODES=""
    local CTQ_OUT="${CTQ_DIR}/${SAMPLE_ID}"
    local VIRAL_OUT="${VIRAL_DIR}/${SAMPLE_ID}"

    # Chemistry-specific STAR parameters (both use Forward — L7 + VaD template confirmed)
    local WHITELIST UMILEN
    if [[ "${CHEMISTRY}" == "3prime" ]]; then
        WHITELIST="${WL_3M}"
        UMILEN=12
        log "Chemistry: 3' v3 (soloUMIlen=12, 3M whitelist, Forward)"
    else
        WHITELIST="${WL_737K}"
        UMILEN=10
        log "Chemistry: 5' v2 (soloUMIlen=10, 737K whitelist, Forward)"
    fi

    mkdir -p "${LOG_DIR}" "${STAR_OUT}" "${CTQ_OUT}" "${VIRAL_OUT}" "${LOCK_DIR}" "${FASTQ_DIR}"

    if [[ -f "${CTQ_OUT}/stage2/matrix.mtx" && -f "${VIRAL_OUT}/viral_counts_t1.tsv" ]]; then
        log "SKIP — already complete"
        return 0
    fi

    # ── STAGE 1: DOWNLOAD — all SRRs for this sample ─────────────────────────
    local all_r1="" all_r2="" first_r1=""

    for srr in $(echo "${SRR_CSV}" | tr ',' ' '); do
        download_srr "${srr}" "${SAMPLE_ID}"
        if [[ -z "${DETECTED_R1}" || -z "${DETECTED_R2}" ]]; then
            log "ERROR: download_srr failed for ${srr}"
            return 1
        fi
        all_r1="${all_r1:+${all_r1},}${DETECTED_R1}"
        all_r2="${all_r2:+${all_r2},}${DETECTED_R2}"
        [[ -z "${first_r1}" ]] && first_r1="${DETECTED_R1}"
    done

    local n_srr
    n_srr=$(echo "${SRR_CSV}" | tr ',' '\n' | wc -l)
    log "DOWNLOAD done — ${n_srr} SRRs assembled for ${SAMPLE_ID}"

    # ── FASTQ VERIFY (check first SRR's R1 only — sufficient to detect swap) ─
    local r1_len
    r1_len=$(zcat "${first_r1}" | awk 'NR%4==2{print length($0); exit}')
    log "FASTQ verify — R1 length: ${r1_len} bp (expect 26-28 for 5' v2, 28 for 3' v3)"

    if [[ "${r1_len}" -lt 20 ]]; then
        log "ERROR: R1 too short (${r1_len} bp) — wrong layout. Aborting."
        return 1
    fi

    # R1 pos-1 T fraction: genuine CB reads show <70% T; cDNA reads show ≥70% (L5 / L5b)
    local r1_t_count
    r1_t_count=$(zcat "${first_r1}" | awk 'NR%4==2 && NR<=400{print substr($0,1,1)}' | grep -c "^T$" || true)
    if [[ "${r1_t_count}" -ge 70 ]]; then
        log "ERROR: R1 pos-1 T=${r1_t_count}/100 (≥70%) — likely cDNA reads, R1/R2 may be swapped. Aborting."
        return 1
    fi
    log "FASTQ verify passed — R1 pos-1 T: ${r1_t_count}/100"

    # ── STAGE 2: STAR ALIGNMENT ──────────────────────────────────────────────
    # Remove stale zero-byte BAM from a previously failed run
    if [[ -f "${BAM}" && ! -s "${BAM}" ]]; then
        rm -f "${BAM}"
        rm -rf "${STAR_OUT}/${SAMPLE_ID}_Solo.out"
        log "Removed stale zero-byte BAM"
    fi

    if [[ ! -f "${BAM}" ]]; then
        # outTmpDir on EFS — compute node root EBS is tiny (L5c); use EFS for sort temp
        local star_tmp_id="${SLURM_JOB_ID:-$$}"
        local star_tmp_dir="${WORKDIR}/tmp/STAR_${SAMPLE_ID}_${star_tmp_id}"
        rm -rf "${star_tmp_dir}" 2>/dev/null || true
        mkdir -p "${WORKDIR}/tmp"

        log "ALIGN start — soloUMIlen=${UMILEN} soloFeatures=GeneFull_Ex50pAS soloStrand=Forward"
        log "  R2 files: ${all_r2}"
        log "  R1 files: ${all_r1}"

        "${STAR_BIN}" \
            --runThreadN "${STAR_THREADS}" \
            --genomeDir "${STAR_INDEX}" \
            --readFilesIn "${all_r2}" "${all_r1}" \
            --readFilesCommand zcat \
            --soloType CB_UMI_Simple \
            --soloCBwhitelist "${WHITELIST}" \
            --soloCBstart 1 --soloCBlen 16 \
            --soloUMIstart 17 --soloUMIlen "${UMILEN}" \
            --soloBarcodeReadLength 0 \
            --soloStrand Forward \
            --soloFeatures GeneFull_Ex50pAS \
            --soloCellFilter EmptyDrops_CR \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI AS NM CB UB \
            --outSAMmultNmax -1 \
            --outFileNamePrefix "${STAR_OUT}/${SAMPLE_ID}_" \
            --outTmpDir "${star_tmp_dir}" \
            >> "${LOG_DIR}/${SAMPLE_ID}.log" 2>&1

        local star_exit=$?
        if [[ ${star_exit} -ne 0 ]]; then
            log "ERROR: STAR failed (exit ${star_exit})"
            return 1
        fi

        local uniq_pct
        uniq_pct=$(grep "Uniquely mapped reads %" "${STAR_OUT}/${SAMPLE_ID}_Log.final.out" 2>/dev/null | awk '{print $NF}')
        log "ALIGN done — uniquely mapped: ${uniq_pct}"

        # Delete all per-SRR FASTQs immediately after successful alignment
        for srr in $(echo "${SRR_CSV}" | tr ',' ' '); do
            rm -f "${FASTQ_DIR}/${SAMPLE_ID}_${srr}_R1.fastq.gz" \
                  "${FASTQ_DIR}/${SAMPLE_ID}_${srr}_R2.fastq.gz"
        done
        rm -rf "${star_tmp_dir}" 2>/dev/null || true
        log "FASTQs deleted (${SAMPLE_ID})"
    else
        log "ALIGN skip — BAM present"
        # Clean up any stale FASTQs even on skip
        for srr in $(echo "${SRR_CSV}" | tr ',' ' '); do
            rm -f "${FASTQ_DIR}/${SAMPLE_ID}_${srr}_R1.fastq.gz" \
                  "${FASTQ_DIR}/${SAMPLE_ID}_${srr}_R2.fastq.gz" 2>/dev/null || true
        done
    fi

    # Resolve barcode path (GeneFull_Ex50pAS — L8)
    if [[ -f "${BARCODES_FILTERED}" ]]; then
        BARCODES="${BARCODES_FILTERED}"
        log "Barcodes (filtered): $(wc -l < ${BARCODES}) cells"
    elif [[ -f "${BARCODES_RAW}" ]]; then
        BARCODES="${BARCODES_RAW}"
        log "WARNING: using raw barcodes ($(wc -l < ${BARCODES}))"
    else
        log "ERROR: no barcodes.tsv under GeneFull_Ex50pAS — check STAR output"
        return 1
    fi

    # ── INDEX BAM (required before ctq and viral — L43c) ─────────────────────
    if [[ ! -f "${BAM}.bai" ]]; then
        log "INDEX BAM"
        "${SAMTOOLS}" index -@ 4 "${BAM}" 2>>"${LOG_DIR}/${SAMPLE_ID}.log"
    fi
    if [[ ! -f "${BAM}.bai" ]]; then
        log "ERROR: BAM indexing failed"
        return 1
    fi

    # ── STAGE 3: scTEQuant stages 1+2 ────────────────────────────────────────
    if [[ ! -f "${CTQ_OUT}/stage2/matrix.mtx" ]]; then
        rm -rf "${CTQ_OUT}/stage1/_tmp_chrom_gtfs_s1" 2>/dev/null || true  # stale temp (L15)
        local n_cells
        n_cells=$(wc -l < "${BARCODES}")
        log "scTEQuant start — ${n_cells} cells"
        sem_acquire "ctq" "${N_CTQ_MAX}"
        if "${CTQ_PY}" "${CTQ_SCRIPT}" \
            --bam "${BAM}" --gtf "${GTF}" \
            --fragment-to-tu "${FRAG_MAP}" \
            --barcodes "${BARCODES}" \
            --output-dir "${CTQ_OUT}" --n-cells "${n_cells}" \
            >> "${LOG_DIR}/${SAMPLE_ID}.log" 2>&1; then
            log "scTEQuant done"
        else
            log "scTEQuant FAILED — continuing to viral"  # non-fatal: allows BAM to be inspected
        fi
        sem_release "ctq"
    else
        log "scTEQuant skip — stage2/matrix.mtx present"
    fi

    # ── STAGE 4: VIRAL — human-tropic panel ──────────────────────────────────
    if [[ ! -f "${VIRAL_OUT}/viral_counts_t1.tsv" ]]; then
        log "VIRAL start — human-tropic panel (EBV, KSHV, HPV, CMV, HHV priority)"
        "${CTQ_BIN}" viral count-reads \
            --bam "${BAM}" --barcodes "${BARCODES}" \
            --min-reads 1 --min-cells 1 \
            --output "${VIRAL_OUT}/viral_counts_raw_t1.tsv" \
            --report "${VIRAL_OUT}/viral_detection_t1.tsv" \
            >> "${LOG_DIR}/${SAMPLE_ID}.log" 2>&1

        # Post-filter to human_tropic_viruses.tsv whitelist
        "${CTQ_PY}" -c "
import pandas as pd
df = pd.read_csv('${VIRAL_OUT}/viral_counts_raw_t1.tsv', sep='\t', index_col=0)
wl_df = pd.read_csv('${VIRAL_WL}', sep='\t')
wl = wl_df['accession'].tolist()
keep = [c for c in df.columns if c in wl]
df[keep].to_csv('${VIRAL_OUT}/viral_counts_t1.tsv', sep='\t')
print(f'Viral filter: {len(df.columns)} total -> {len(keep)} human-tropic retained')
if 'name' in wl_df.columns:
    name_map = dict(zip(wl_df['accession'], wl_df['name']))
    print('Retained:', [name_map.get(a, a) for a in keep])
" >> "${LOG_DIR}/${SAMPLE_ID}.log" 2>&1
        log "VIRAL done"
    else
        log "VIRAL skip — viral_counts_t1.tsv present"
    fi

    # ── STAGE 5: S3 UPLOAD + BAM CLEANUP ─────────────────────────────────────
    # Gate BAM deletion on confirmed downstream outputs (L36 pattern)
    if [[ -f "${CTQ_OUT}/stage2/matrix.mtx" && -f "${VIRAL_OUT}/viral_counts_t1.tsv" ]]; then
        log "S3 upload start"
        aws s3 sync "${STAR_OUT}/${SAMPLE_ID}_Solo.out/" \
            "${S3_OUT}/starsolo/${SAMPLE_ID}/" --no-progress \
            >> "${LOG_DIR}/${SAMPLE_ID}.log" 2>&1
        aws s3 sync "${CTQ_OUT}/" \
            "${S3_OUT}/sctequant/${SAMPLE_ID}/" --no-progress \
            >> "${LOG_DIR}/${SAMPLE_ID}.log" 2>&1
        aws s3 sync "${VIRAL_OUT}/" \
            "${S3_OUT}/viral/${SAMPLE_ID}/" --no-progress \
            >> "${LOG_DIR}/${SAMPLE_ID}.log" 2>&1
        log "S3 upload done"
        rm -f "${BAM}" "${BAM}.bai"
        log "CLEANUP done — BAM deleted"
    else
        log "WARNING: BAM retained — ctq or viral output missing; re-run to retry without re-downloading"
    fi

    log "SAMPLE COMPLETE: ${SAMPLE_ID}"
}

# ── Mode dispatch ──────────────────────────────────────────────────────────────
MODE="${1:-submit}"

if [[ "${MODE}" == "run" ]]; then
    # Called by sbatch
    SAMPLE_ID="$2"
    SRR_CSV="$3"
    CHEMISTRY="$4"
    run_sample "${SAMPLE_ID}" "${SRR_CSV}" "${CHEMISTRY}"

elif [[ "${MODE}" == "submit" ]]; then
    mkdir -p "${LOG_DIR}" "${LOCK_DIR}" "${TMP_DIR}" "${FASTQ_DIR}"

    # Pre-flight: all reference files must exist before submitting
    for f in "${WL_3M}" "${WL_737K}" "${GTF}" "${FRAG_MAP}" "${VIRAL_WL}"; do
        if [[ ! -f "$f" ]]; then
            echo "ERROR: missing reference file: $f" >&2
            echo "Run: bash ${EFS}/gse189889_slurm.sh sync-refs" >&2
            exit 1
        fi
    done

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Submitting GSE189889 acral melanoma (10 jobs)"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Output: ${WORKDIR}"

    # Format: "SAMPLE_ID SRR_CSV CHEMISTRY"
    # No declare -A — bash 3.2 / CentOS 7 incompatible (L35)
    for args in \
        "AM1 SRR17075055,SRR17075056,SRR17075057,SRR17075058,SRR17075059,SRR17075060,SRR17075061,SRR17075062 3prime" \
        "AM2 SRR17075063,SRR17075064,SRR17075065,SRR17075066,SRR17075067,SRR17075068,SRR17075069,SRR17075070 5prime" \
        "AM3 SRR17075071,SRR17075072,SRR17075073,SRR17075074,SRR17075075,SRR17075076,SRR17075077,SRR17075078 5prime" \
        "AM4 SRR17075079,SRR17075080,SRR17075081,SRR17075082,SRR17075083,SRR17075084,SRR17075085,SRR17075086 5prime" \
        "AM5 SRR17075087,SRR17075088,SRR17075089,SRR17075090,SRR17075091,SRR17075092,SRR17075093,SRR17075094 5prime" \
        "AM6 SRR17075095,SRR17075096,SRR17075097,SRR17075098,SRR17075099,SRR17075100,SRR17075101,SRR17075102 5prime" \
        "AM7 SRR17075103,SRR17075104,SRR17075105,SRR17075106,SRR17075107,SRR17075108,SRR17075109,SRR17075110 5prime" \
        "AM8_Node SRR17075111,SRR17075112,SRR17075113,SRR17075114 5prime" \
        "AM8_Toe SRR17075115,SRR17075116,SRR17075117,SRR17075118,SRR17075119,SRR17075120,SRR17075121,SRR17075122 5prime"
    do
        sample=$(echo "$args" | awk '{print $1}')
        srr_csv=$(echo "$args" | awk '{print $2}')
        chemistry=$(echo "$args" | awk '{print $3}')

        if [[ -f "${CTQ_DIR}/${sample}/stage2/matrix.mtx" && \
              -f "${VIRAL_DIR}/${sample}/viral_counts_t1.tsv" ]]; then
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] SKIP ${sample} — already complete"
            continue
        fi

        JOB_ID=$(sbatch \
            --partition=defq \
            --cpus-per-task=16 \
            --exclusive \
            --job-name="gse189889_${sample}" \
            --output="${LOG_DIR}/slurm_${sample}_%j.log" \
            --error="${LOG_DIR}/slurm_${sample}_%j.err" \
            --wrap="bash ${EFS}/gse189889_slurm.sh run ${sample} ${srr_csv} ${chemistry}" \
            | awk '{print $NF}')

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Submitted ${sample} (${srr_csv}) → job ${JOB_ID}"
    done

    echo ""
    echo "Monitor: squeue -u jmarston"
    echo "Status:  bash ${EFS}/gse189889_slurm.sh status"
    echo "Logs:    ${LOG_DIR}/"

elif [[ "${MODE}" == "sync-refs" ]]; then
    echo "Checking reference files for GSE189889..."
    echo ""

    # All references should already exist from atherosclerosis/gse282111 projects.
    # The only potential new download is 737K whitelist (5' v2) — usually already present.
    if [[ ! -f "${WL_737K}" ]]; then
        echo "Downloading 737K whitelist (5' v2 / 3' v2) from S3..."
        aws s3 cp s3://jez-research-data/reference/whitelists/737K-august-2016.txt.gz \
            /tmp/737K.txt.gz --no-progress 2>/dev/null \
        || curl -L "https://cf.10xgenomics.com/supp/cell-exp/737K-august-2016.txt.gz" \
            -o /tmp/737K.txt.gz
        gunzip -c /tmp/737K.txt.gz > "${WL_737K}"
        rm -f /tmp/737K.txt.gz
        echo "737K whitelist ready: $(wc -l < ${WL_737K}) barcodes (expect 737280)"
    fi

    all_ok=true
    for label_colon_path in \
        "3M whitelist (3' v3):${WL_3M}" \
        "737K whitelist (5' v2):${WL_737K}" \
        "ctq GTF:${GTF}" \
        "frag-to-TU map:${FRAG_MAP}" \
        "viral whitelist:${VIRAL_WL}" \
        "STAR index:${STAR_INDEX}"
    do
        label=$(echo "$label_colon_path" | cut -d: -f1)
        fpath=$(echo "$label_colon_path" | cut -d: -f2-)
        if [[ -e "$fpath" ]]; then
            echo "  OK      ${label}"
        else
            echo "  MISSING ${label}: ${fpath}" >&2
            all_ok=false
        fi
    done

    if [[ "${all_ok}" == "true" ]]; then
        echo ""
        echo "All reference files ready. Run: bash ${EFS}/gse189889_slurm.sh submit"
    else
        echo "" >&2
        echo "ERROR: missing references above must be resolved before submitting." >&2
        exit 1
    fi

elif [[ "${MODE}" == "status" ]]; then
    echo "=== Slurm queue ==="
    squeue -u jmarston -o "%.10i %.20j %.8T %.10M %.6C" 2>/dev/null || true

    echo ""
    echo "=== GSE189889 sample status ==="
    for args in \
        "AM1 3prime" \
        "AM2 5prime" \
        "AM3 5prime" \
        "AM4 5prime" \
        "AM5 5prime" \
        "AM6 5prime" \
        "AM7 5prime" \
        "AM8_Node 5prime" \
        "AM8_Toe 5prime"
    do
        sample=$(echo "$args" | awk '{print $1}')
        chem=$(echo "$args" | awk '{print $2}')
        bam="${STAR_DIR}/${sample}/${sample}_Aligned.sortedByCoord.out.bam"
        if [[ -f "${CTQ_DIR}/${sample}/stage2/matrix.mtx" && \
              -f "${VIRAL_DIR}/${sample}/viral_counts_t1.tsv" ]]; then
            echo "  ${sample} (${chem}): COMPLETE"
        elif [[ -f "${bam}" ]]; then
            echo "  ${sample} (${chem}): aligned — ctq/viral pending"
        elif ls "${FASTQ_DIR}/${sample}_"*.fastq.gz > /dev/null 2>&1; then
            echo "  ${sample} (${chem}): downloading"
        else
            echo "  ${sample} (${chem}): not started"
        fi
    done

    echo ""
    echo "=== Recent log activity ==="
    for f in "${LOG_DIR}"/*.log; do
        [[ -f "$f" ]] || continue
        sname=$(basename "$f" .log)
        last=$(grep -v "join\|spot\|[0-9]\.[0-9]*%" "$f" 2>/dev/null | tail -1)
        echo "  ${sname}: ${last:-no log yet}"
    done

else
    echo "Usage: $0 {submit|run <sample> <srr_csv> <chemistry>|status|sync-refs}" >&2
    exit 1
fi
