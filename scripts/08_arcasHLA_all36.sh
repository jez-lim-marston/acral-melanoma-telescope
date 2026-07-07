#!/usr/bin/env bash
# 08_arcasHLA_all36.sh
# Type HLA class I + II for ALL 36 RNA-seq samples with arcasHLA for methodological
# consistency (replacing the partial OptiType call set that is missing 8 samples).
#
# Inputs (auto-detected, in priority order):
#   1. FASTQ_DIR env var pointing to a folder containing ${SRR}_1.fastq.gz and
#      ${SRR}_2.fastq.gz (or unsuffixed .fastq.gz → will be re-paired with fasterq-dump).
#   2. BAM_DIR env var pointing to a folder containing ${SRR}.bam (arcasHLA extract
#      works directly from a GRCh38 BAM).
#   3. Neither set → stream-download paired FASTQs from ENA per sample.
#
# Usage:
#   FASTQ_DIR=/Volumes/MyDrive/acral_fastq  bash scripts/08_arcasHLA_all36.sh
#   BAM_DIR=/Volumes/MyDrive/acral_bams      bash scripts/08_arcasHLA_all36.sh
#   bash scripts/08_arcasHLA_all36.sh                           # ENA fallback
#
# Output:
#   results/hla_arcas/${SRR}/${SRR}.genotype.json    per-sample (class I + II)
#   results/hla_arcas/arcas_calls.tsv                merged 36-row table
#   results/hla_arcas/arcas_vs_optitype_concordance.tsv  sanity check on 28 overlap

set -euo pipefail

PROJECT=/Users/jezmarston/Desktop/Bulk_Acral_Claude
BT2_INDEX=/Users/jezmarston/Desktop/HERVK_113_115_Population_Analysis/data/grch38_bt2/GRCh38
OUTDIR=${PROJECT}/results/hla_arcas
TMPDIR=${PROJECT}/data/hla_tmp
LOGDIR=${PROJECT}/logs
THREADS=4

mkdir -p "${OUTDIR}" "${TMPDIR}" "${LOGDIR}"

# All 36 Liang 2017 RNA-seq accessions, mapped to PATIENT IDs.
ALL_SRRS=(
  SRR3882763:1   SRR3882768:2   SRR3882773:3   SRR3882778:4   SRR3882783:5
  SRR3882788:6   SRR3882793:7   SRR3882800:9   SRR3882805:10  SRR3882810:11
  SRR3882815:12  SRR3882820:13  SRR3882825:14  SRR3882830:15  SRR3882835:16
  SRR3882840:17  SRR3882845:18  SRR3882849:19  SRR3882853:20  SRR3882857:21
  SRR3882861:22  SRR3882864:23  SRR3882868:24  SRR3882873:25a SRR3882876:25b
  SRR3882880:26  SRR3882883:27  SRR3882888:28  SRR3882892:29a SRR3882897:29c
  SRR3882901:30  SRR3882905:31  SRR3882910:32  SRR3882913:33  SRR3882918:34a
  SRR3882920:34b
)

# ---------- env setup (one-time) ----------------------------------------------
ENV_NAME=arcas
if ! mamba env list | grep -q "^${ENV_NAME} "; then
  echo "[setup] creating ${ENV_NAME} env"
  mamba create -y -n "${ENV_NAME}" -c bioconda -c conda-forge \
    arcas-hla bowtie2 samtools=1.19 sra-tools pigz
fi
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${ENV_NAME}"

if [[ ! -f "${OUTDIR}/.reference_ok" ]]; then
  arcasHLA reference --update
  touch "${OUTDIR}/.reference_ok"
fi

# ---------- helpers -----------------------------------------------------------
ena_base() {
  local SRR=$1
  echo "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${SRR:0:6}/00${SRR: -1}/${SRR}"
}

prepare_extracted_reads() {
  # Populates ${SDIR}/${SRR}.extracted.{1,2}.fq.gz from whichever input source is available.
  local SRR=$1 SDIR=$2

  # Case A: existing BAM → arcasHLA extract directly
  if [[ -n "${BAM_DIR:-}" && -f "${BAM_DIR}/${SRR}.bam" ]]; then
    echo "[${SRR}] using BAM from ${BAM_DIR}"
    [[ -f "${BAM_DIR}/${SRR}.bam.bai" ]] || samtools index "${BAM_DIR}/${SRR}.bam"
    caffeinate -i -s arcasHLA extract --paired -t ${THREADS} \
      -o "${SDIR}" "${BAM_DIR}/${SRR}.bam"
    return
  fi

  # Case B: existing paired FASTQ → skip download, align, extract
  local R1 R2
  if [[ -n "${FASTQ_DIR:-}" ]]; then
    if   [[ -f "${FASTQ_DIR}/${SRR}_1.fastq.gz" && -f "${FASTQ_DIR}/${SRR}_2.fastq.gz" ]]; then
      R1="${FASTQ_DIR}/${SRR}_1.fastq.gz"; R2="${FASTQ_DIR}/${SRR}_2.fastq.gz"
    elif [[ -f "${FASTQ_DIR}/${SRR}_1.fastq"    && -f "${FASTQ_DIR}/${SRR}_2.fastq"    ]]; then
      R1="${FASTQ_DIR}/${SRR}_1.fastq";    R2="${FASTQ_DIR}/${SRR}_2.fastq"
    fi
  fi

  # Case C: download from ENA
  if [[ -z "${R1:-}" ]]; then
    echo "[${SRR}] downloading from ENA"
    R1=${TMPDIR}/${SRR}_1.fastq.gz
    R2=${TMPDIR}/${SRR}_2.fastq.gz
    local URL; URL=$(ena_base "${SRR}")
    caffeinate -i -s curl -L -C - --max-time 14400 -o "${R1}" "${URL}/${SRR}_1.fastq.gz"
    caffeinate -i -s curl -L -C - --max-time 14400 -o "${R2}" "${URL}/${SRR}_2.fastq.gz"
  fi

  # Align → sort → index → arcasHLA extract
  local BAM=${TMPDIR}/${SRR}.sorted.bam
  caffeinate -i -s bash -c "
    bowtie2 --very-sensitive -p ${THREADS} -x ${BT2_INDEX} \
      -1 '${R1}' -2 '${R2}' 2> '${SDIR}/${SRR}.bt2.log' \
      | samtools sort -@ ${THREADS} -o '${BAM}' -
    samtools index '${BAM}'
  "
  caffeinate -i -s arcasHLA extract --paired -t ${THREADS} \
    -o "${SDIR}" "${BAM}"

  # Cleanup ENA-downloaded intermediates only
  if [[ "${R1}" == "${TMPDIR}"* ]]; then rm -f "${R1}" "${R2}"; fi
  rm -f "${BAM}" "${BAM}.bai"
}

# ---------- per-sample loop ---------------------------------------------------
for PAIR in "${ALL_SRRS[@]}"; do
  SRR=${PAIR%%:*}
  PID=${PAIR##*:}
  SDIR=${OUTDIR}/${SRR}
  mkdir -p "${SDIR}"

  if [[ -f "${SDIR}/${SRR}.genotype.json" ]]; then
    echo "[skip] ${SRR} (Pt ${PID}) already genotyped"
    continue
  fi

  echo "=============================================================="
  echo "[${SRR}] Pt ${PID} — start $(date)"
  echo "=============================================================="

  prepare_extracted_reads "${SRR}" "${SDIR}"

  caffeinate -i -s arcasHLA genotype \
    "${SDIR}/${SRR}.extracted.1.fq.gz" \
    "${SDIR}/${SRR}.extracted.2.fq.gz" \
    -g A,B,C,DPB1,DQA1,DQB1,DRB1 \
    -t ${THREADS} \
    -o "${SDIR}" 2>&1 | tee "${LOGDIR}/arcas_${SRR}.log"

  rm -f "${SDIR}/${SRR}.extracted.1.fq.gz" "${SDIR}/${SRR}.extracted.2.fq.gz"
  echo "[${SRR}] Pt ${PID} — done $(date)"
done

# ---------- merge JSONs → arcas_calls.tsv ------------------------------------
python3 - <<'PYEOF'
import json, glob, os, csv
OUT = "/Users/jezmarston/Desktop/Bulk_Acral_Claude/results/hla_arcas"
GENES = ["A","B","C","DPB1","DQA1","DQB1","DRB1"]
rows = []
for jf in sorted(glob.glob(f"{OUT}/*/*.genotype.json")):
    srr = os.path.basename(jf).split(".")[0]
    d = json.load(open(jf))
    row = {"sample_id": srr}
    for g in GENES:
        alleles = d.get(g, [])
        row[f"{g}1"] = alleles[0] if len(alleles) > 0 else ""
        row[f"{g}2"] = alleles[1] if len(alleles) > 1 else ""
    rows.append(row)
fields = ["sample_id"] + [f"{g}{i}" for g in GENES for i in (1,2)]
with open(f"{OUT}/arcas_calls.tsv","w",newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
    w.writeheader(); w.writerows(rows)
print(f"wrote {OUT}/arcas_calls.tsv ({len(rows)} rows)")
PYEOF

# ---------- concordance vs original OptiType on 28 overlap -------------------
python3 - <<'PYEOF'
import pandas as pd, os
OUT = "/Users/jezmarston/Desktop/Bulk_Acral_Claude/results/hla_arcas"
opti_path = "/Users/jezmarston/Desktop/Weill Cornell/Dermatology Projects/Acral Melanoma Projects/Acral Melanom Bulk RNA-seq Paper/HLA_Metadata.xlsx"
arcas = pd.read_csv(f"{OUT}/arcas_calls.tsv", sep="\t")
opti  = pd.read_excel(opti_path).rename(columns={"Unnamed: 0":"sample_id"})
# arcasHLA emits e.g. "A*02:01:01:01"; OptiType emits "A*02:01". Truncate arcas to 2-field.
def two_field(x):
    if not isinstance(x,str) or x == "": return ""
    parts = x.split(":")
    return ":".join(parts[:2])
for c in arcas.columns:
    if c != "sample_id": arcas[c] = arcas[c].map(two_field)
m = arcas.merge(opti, on="sample_id", suffixes=("_arcas","_opti"))
rows = []
for _, r in m.iterrows():
    for locus in ["A","B","C"]:
        a = sorted([r[f"{locus}1_arcas"], r[f"{locus}2_arcas"]])
        o = sorted([r[f"{locus}1_opti"],  r[f"{locus}2_opti"]])
        rows.append({"sample_id": r["sample_id"], "locus": locus,
                     "arcas": "/".join(a), "opti": "/".join(o),
                     "match": a == o})
conc = pd.DataFrame(rows)
conc.to_csv(f"{OUT}/arcas_vs_optitype_concordance.tsv", sep="\t", index=False)
print(conc.groupby("locus")["match"].agg(["sum","count"]))
print(f"wrote {OUT}/arcas_vs_optitype_concordance.tsv")
PYEOF

echo "[done] See ${OUTDIR}/arcas_calls.tsv and arcas_vs_optitype_concordance.tsv"
