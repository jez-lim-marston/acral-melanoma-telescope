#!/usr/bin/env bash
# 02_align_telescope.sh
#
# For each of the 36 Liang 2017 (PRJNA304068) acral melanoma RNA-seq samples:
#   1. Download paired-end FASTQs from ENA (EBI FTP, resumable)
#   2. Align with bowtie2 --very-sensitive -k 100 (multi-mapping; required for Telescope)
#   3. Sort and index the BAM
#   4. Run telescope assign (best_random reassignment)
#   5. Delete FASTQ and BAM to conserve disk (stream-and-delete strategy)
#
# Output per sample:
#   results/telescope/${SRR}/${SRR}_telescope_report.tsv   ← main count table
#   results/telescope/${SRR}/${SRR}_telescope_stats.tsv
#
# Disk peak: ~1 FASTQ pair (~12 GB) + BAM (~4 GB) at any one time.
# Wall time:  ~3–4 h per sample; samples are processed sequentially by default.
#             Set PARALLEL_JOBS > 1 and ensure sufficient RAM (8 GB per job).
#
# Usage:
#   bash scripts/02_align_telescope.sh
#
#   To process a single sample:
#   SRR=SRR3882763  bash scripts/02_align_telescope.sh
#
# Requirements (PATH):
#   bowtie2    (≥ 2.4)   samtools   (≥ 1.16)
#   telescope  (≥ 1.0.3, github.com/mlbendall/telescope)
#   caffeinate (macOS, prevents sleep)   -- remove on Linux

set -euo pipefail

PROJECT="$(cd "$(dirname "$0")/.." && pwd)"
REF="${PROJECT}/data/reference"
OUTDIR="${PROJECT}/results/telescope"
TMPDIR="${PROJECT}/data/tmp_align"
LOGDIR="${PROJECT}/logs"
THREADS=4
PARALLEL_JOBS=1    # increase if RAM allows (each job uses ~8 GB)

BT2_INDEX="${REF}/bowtie2_index/GRCh38"
TELESCOPE_GTF="${REF}/HERVK_telescope.gtf"

mkdir -p "${OUTDIR}" "${TMPDIR}" "${LOGDIR}"

# 36 RNA-seq SRR accessions from PRJNA304068 (Liang et al. 2017)
# Format: SRR:PATIENT_ID
ALL_SRRS=(
  SRR3882763:1    SRR3882768:2    SRR3882773:3    SRR3882778:4    SRR3882783:5
  SRR3882788:6    SRR3882793:7    SRR3882800:9    SRR3882805:10   SRR3882810:11
  SRR3882815:12   SRR3882820:13   SRR3882825:14   SRR3882830:15   SRR3882835:16
  SRR3882840:17   SRR3882845:18   SRR3882849:19   SRR3882853:20   SRR3882857:21
  SRR3882861:22   SRR3882864:23   SRR3882868:24   SRR3882873:25a  SRR3882876:25b
  SRR3882880:26   SRR3882883:27   SRR3882888:28   SRR3882892:29a  SRR3882897:29c
  SRR3882901:30   SRR3882905:31   SRR3882910:32   SRR3882913:33   SRR3882918:34a
  SRR3882920:34b
)

# ── Helpers ───────────────────────────────────────────────────────────────────
log() { echo "[$(date '+%H:%M:%S')] $*"; }

ena_url() {
  local SRR=$1 SUFFIX=$2
  local DIR6="${SRR:0:6}"
  local PAD="00${SRR: -1}"   # e.g. SRR3882763 → 003
  echo "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${DIR6}/${PAD}/${SRR}/${SRR}${SUFFIX}"
}

download_fastq() {
  local SRR=$1 OUTDIR=$2
  local R1="${OUTDIR}/${SRR}_1.fastq.gz"
  local R2="${OUTDIR}/${SRR}_2.fastq.gz"

  for SUFFIX in "_1.fastq.gz" "_2.fastq.gz"; do
    local DEST="${OUTDIR}/${SRR}${SUFFIX}"
    local URL; URL=$(ena_url "${SRR}" "${SUFFIX}")
    if [[ ! -f "${DEST}" ]]; then
      log "[${SRR}] Downloading ${URL}"
      caffeinate -i -s curl -L -C - --max-time 14400 --retry 3 \
        -o "${DEST}" "${URL}"
    fi
  done
  echo "${R1}" "${R2}"
}

process_sample() {
  local SRR=$1 PID=$2
  local SDIR="${OUTDIR}/${SRR}"
  local REPORT="${SDIR}/${SRR}_telescope_report.tsv"

  if [[ -f "${REPORT}" ]]; then
    log "[skip] ${SRR} (Pt ${PID}) — Telescope report already exists"
    return
  fi

  mkdir -p "${SDIR}"
  log "════════════════════════════════════════════"
  log "[${SRR}] Pt ${PID} — start $(date)"

  # 1. Download FASTQs
  local FQDIR="${TMPDIR}/${SRR}"
  mkdir -p "${FQDIR}"
  download_fastq "${SRR}" "${FQDIR}"
  local R1="${FQDIR}/${SRR}_1.fastq.gz"
  local R2="${FQDIR}/${SRR}_2.fastq.gz"

  # 2. Bowtie2 alignment (multi-mapping, -k 100)
  local BAM="${TMPDIR}/${SRR}.sorted.bam"
  log "[${SRR}] bowtie2 align → sort"
  caffeinate -i -s bash -c "
    bowtie2 \
      --very-sensitive \
      -k 100 \
      -p ${THREADS} \
      --no-unal \
      -x '${BT2_INDEX}' \
      -1 '${R1}' -2 '${R2}' \
      2> '${SDIR}/${SRR}.bt2.log' \
    | samtools sort -@ ${THREADS} -o '${BAM}' -
  "
  samtools index "${BAM}"
  log "[${SRR}] BAM ready: ${BAM}"

  # Delete FASTQs to free disk
  rm -f "${R1}" "${R2}"
  rmdir --ignore-fail-on-non-empty "${FQDIR}" 2>/dev/null || true

  # 3. Telescope assign
  log "[${SRR}] telescope assign"
  caffeinate -i -s telescope assign \
    --reassign_mode best_random \
    --theta_prior 200000 \
    --max_iter 200 \
    --updated_sam \
    "${BAM}" \
    "${TELESCOPE_GTF}" \
    --outdir "${SDIR}" \
    --exp_tag "${SRR}" \
    2>&1 | tee "${LOGDIR}/telescope_${SRR}.log"

  # Delete BAM to free disk
  rm -f "${BAM}" "${BAM}.bai"

  log "[${SRR}] Pt ${PID} — done $(date)"
}

# ── Main loop ─────────────────────────────────────────────────────────────────
# If SRR env var is set, process only that sample (for testing)
if [[ -n "${SRR:-}" ]]; then
  PID=$(printf '%s\n' "${ALL_SRRS[@]}" | grep "^${SRR}:" | cut -d: -f2)
  process_sample "${SRR}" "${PID:-unknown}"
  exit 0
fi

# Sequential processing (default; safest for disk and RAM)
for PAIR in "${ALL_SRRS[@]}"; do
  SRR="${PAIR%%:*}"
  PID="${PAIR##*:}"
  process_sample "${SRR}" "${PID}"
done

log "All 36 samples processed. Telescope reports in ${OUTDIR}/"
log "Next step: bash scripts/03_merge_telescope_counts.py"
