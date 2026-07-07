#!/usr/bin/env bash
# 04_kallisto_gene_counts.sh
#
# Quantifies transcript-level expression for all 36 PRJNA304068 RNA-seq samples
# using kallisto, then aggregates to a per-sample TPM matrix (Gene_AM.csv).
#
# Pipeline per sample:
#   1. Download paired-end FASTQs from ENA (if not already present)
#   2. kallisto quant → transcript-level TPM estimates
#   3. Delete FASTQs (stream-and-delete; Gene_AM.csv is ~20 MB)
#
# Final merge:
#   4. Concatenate per-sample abundance.tsv → Gene_AM.csv (transcript × sample)
#      (tx → gene aggregation is performed inside AM_Analysis_Revised.Rmd via biomaRt)
#
# Input:
#   data/reference/kallisto_index.idx  (built by 01_prepare_references.sh)
#
# Output:
#   results/kallisto/{SRR}/abundance.tsv  — per-sample transcript TPMs
#   Acral R/Gene_AM.csv                   — merged transcript × sample matrix
#
# Usage:
#   bash scripts/04_kallisto_gene_counts.sh
#
#   Single sample test:
#   SRR=SRR3882763  bash scripts/04_kallisto_gene_counts.sh
#
# Requirements (PATH):
#   kallisto  (≥ 0.46; quantification is single-threaded for reproducibility)
#   curl  caffeinate  python3

set -euo pipefail

PROJECT="$(cd "$(dirname "$0")/.." && pwd)"
REF="${PROJECT}/data/reference"
OUTDIR="${PROJECT}/results/kallisto"
TMPDIR="${PROJECT}/data/tmp_fastq_kallisto"
LOGDIR="${PROJECT}/logs"
THREADS=4

KALLISTO_INDEX="${REF}/kallisto_index.idx"

mkdir -p "${OUTDIR}" "${TMPDIR}" "${LOGDIR}"

# 36 RNA-seq SRR accessions from PRJNA304068 (Liang et al. 2017)
# Format: SRR:PATIENT_ID  (same order as metadata)
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
  local PAD="00${SRR: -1}"
  echo "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${DIR6}/${PAD}/${SRR}/${SRR}${SUFFIX}"
}

process_sample() {
  local SRR=$1 PID=$2
  local SDIR="${OUTDIR}/${SRR}"
  local DONE="${SDIR}/abundance.tsv"

  if [[ -f "${DONE}" ]]; then
    log "[skip] ${SRR} (Pt ${PID}) — kallisto output already exists"
    return
  fi

  mkdir -p "${SDIR}"
  log "────────────────────────────────────────────"
  log "[${SRR}] Pt ${PID} — start"

  # Download FASTQs
  local FQDIR="${TMPDIR}/${SRR}"
  mkdir -p "${FQDIR}"
  local R1="${FQDIR}/${SRR}_1.fastq.gz"
  local R2="${FQDIR}/${SRR}_2.fastq.gz"

  for SUFFIX in "_1.fastq.gz" "_2.fastq.gz"; do
    local DEST="${FQDIR}/${SRR}${SUFFIX}"
    if [[ ! -f "${DEST}" ]]; then
      local URL; URL=$(ena_url "${SRR}" "${SUFFIX}")
      log "[${SRR}] Downloading ${URL}"
      caffeinate -i -s curl -L -C - --max-time 14400 --retry 3 \
        -o "${DEST}" "${URL}"
    fi
  done

  # kallisto quant (fragment length defaults inferred from paired-end data)
  log "[${SRR}] kallisto quant"
  caffeinate -i -s kallisto quant \
    --index "${KALLISTO_INDEX}" \
    --output-dir "${SDIR}" \
    --threads "${THREADS}" \
    --bootstrap-samples 0 \
    "${R1}" "${R2}" \
    2>&1 | tee "${LOGDIR}/kallisto_${SRR}.log"

  # Delete FASTQs
  rm -f "${R1}" "${R2}"
  rmdir --ignore-fail-on-non-empty "${FQDIR}" 2>/dev/null || true

  log "[${SRR}] done"
}

# ── Main ──────────────────────────────────────────────────────────────────────
if [[ ! -f "${KALLISTO_INDEX}" ]]; then
  echo "ERROR: Kallisto index not found at ${KALLISTO_INDEX}"
  echo "       Run 01_prepare_references.sh first."
  exit 1
fi

if [[ -n "${SRR:-}" ]]; then
  PID=$(printf '%s\n' "${ALL_SRRS[@]}" | grep "^${SRR}:" | cut -d: -f2)
  process_sample "${SRR}" "${PID:-unknown}"
else
  for PAIR in "${ALL_SRRS[@]}"; do
    SRR="${PAIR%%:*}"
    PID="${PAIR##*:}"
    process_sample "${SRR}" "${PID}"
  done
fi

# ── Merge per-sample TPMs into Gene_AM.csv ────────────────────────────────────
log "Merging kallisto outputs → Gene_AM.csv"
python3 - <<'PYEOF'
import os, sys
import pandas as pd

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

SAMPLE_ORDER = [
    "SRR3882763","SRR3882768","SRR3882773","SRR3882778","SRR3882783",
    "SRR3882788","SRR3882793","SRR3882800","SRR3882805","SRR3882810",
    "SRR3882815","SRR3882820","SRR3882825","SRR3882830","SRR3882835",
    "SRR3882840","SRR3882845","SRR3882849","SRR3882853","SRR3882857",
    "SRR3882861","SRR3882864","SRR3882868","SRR3882873","SRR3882876",
    "SRR3882880","SRR3882883","SRR3882888","SRR3882892","SRR3882897",
    "SRR3882901","SRR3882905","SRR3882910","SRR3882913","SRR3882918",
    "SRR3882920",
]

frames = {}
missing = []
for srr in SAMPLE_ORDER:
    path = os.path.join(PROJECT, "results", "kallisto", srr, "abundance.tsv")
    if not os.path.isfile(path):
        missing.append(srr)
        continue
    df = pd.read_csv(path, sep="\t", usecols=["target_id", "tpm"])
    frames[srr] = df.set_index("target_id")["tpm"]

if missing:
    print(f"WARNING: Missing kallisto output for {len(missing)} samples: {missing}",
          file=sys.stderr)

if not frames:
    sys.exit("ERROR: No kallisto outputs found.")

found = [s for s in SAMPLE_ORDER if s in frames]
mat = pd.DataFrame({s: frames[s] for s in found})
mat.index.name = "Gene"    # column name used by AM_Analysis_Revised.Rmd
mat = mat.reset_index()

OUT_CSV = os.path.join(PROJECT, "Acral R", "Gene_AM.csv")
mat.to_csv(OUT_CSV, index=False)
print(f"Wrote {OUT_CSV}  ({mat.shape[0]} transcripts × {len(found)} samples)")
PYEOF

log "Done. Gene_AM.csv written to: ${PROJECT}/Acral R/Gene_AM.csv"
log "Next step: open AM_Analysis_Revised.Rmd in RStudio and knit."
