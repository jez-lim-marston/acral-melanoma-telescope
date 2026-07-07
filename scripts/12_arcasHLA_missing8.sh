#!/usr/bin/env bash
# arcasHLA typing pipeline for the 8 missing-HLA patients (2, 3, 6, 9, 10,
# 11, 16, 28) via dbGaP-authorized download of PRJNA304068 RNA-seq runs.
#
# Per-sample workflow (stream-and-delete; CLAUDE.md §3):
#   prefetch --ngc -> fasterq-dump -> bowtie2 -> samtools sort/index
#   -> arcasHLA extract -> arcasHLA genotype -> rm intermediates
#
# Peak transient disk per sample: ~16 GB (FASTQ + BAM). ~3-4 h / sample.
# Total ~24-30 h for all 8. Wrap with `caffeinate -i -s` to prevent sleep.
#
# Usage:
#   caffeinate -i -s bash scripts/12_arcasHLA_missing8.sh
#
# Resume: safe to re-run; each SRR is skipped if its arcasHLA .json exists.

set -euo pipefail

# --- Config ------------------------------------------------------------------
ROOT=/Users/jezmarston/Desktop/Bulk_Acral_Claude
NGC=/Users/jezmarston/Documents/prj_24438_D27244.ngc
SRR_LIST="${ROOT}/data/missing8_SRRs.txt"
BT2_INDEX="${ROOT}/data/grch38_bt2/GRCh38"
SCRATCH="${HOME}/hla_scratch"
OUT_DIR="${ROOT}/results/arcasHLA"
SRA_DIR="${HOME}/ncbi_dbGaP/sra"
THREADS=4

mkdir -p "${SCRATCH}/fq" "${SCRATCH}/bam" "${OUT_DIR}" "${SRA_DIR}"

# --- Prerequisites -----------------------------------------------------------
# 1. Bowtie2 index symlink
if [[ ! -d "${ROOT}/data/grch38_bt2" ]]; then
  SISTER="/Users/jezmarston/Desktop/HERVK_113_115_Population_Analysis/data/grch38_bt2"
  [[ -d "${SISTER}" ]] || { echo "ERROR: bowtie2 index not found at ${SISTER}" >&2; exit 1; }
  ln -s "${SISTER}" "${ROOT}/data/grch38_bt2"
  echo "Symlinked bowtie2 index from sister project"
fi

# 2. arcasHLA env; use `mamba run -n arcashla` wrapper for every call so PATH
# doesn't need to persist across subshells.
ARCAS="mamba run -n arcashla"
if ! mamba env list | awk '{print $1}' | grep -qx "arcashla"; then
  echo "arcashla mamba env not found; creating..."
  mamba create -y -n arcashla -c bioconda arcas-hla pigz samtools bowtie2 biopython
fi

# 3. IMGTHLA reference (one-time). bioconda arcas-hla 0.6.0 ships without the
# IMGTHLA database; we clone and unzip hla.dat manually, then let arcasHLA
# reference --update build the kallisto index.
ARCAS_SHARE="$(mamba run -n arcashla bash -c 'echo $CONDA_PREFIX')/share/arcas-hla-0.6.0-2"
if [[ ! -f "${ARCAS_SHARE}/dat/IMGTHLA/hla.dat" ]]; then
  echo "Cloning IMGTHLA database (one-time, ~5 min)..."
  cd "${ARCAS_SHARE}/dat"
  rm -rf IMGTHLA
  git clone --depth 1 https://github.com/ANHIG/IMGTHLA.git
  if [[ -f IMGTHLA/hla.dat.zip && ! -f IMGTHLA/hla.dat ]]; then
    (cd IMGTHLA && unzip -o hla.dat.zip)
  fi
  cd - >/dev/null
fi
if [[ ! -f "${ARCAS_SHARE}/dat/ref/hla.fasta" ]]; then
  echo "Building arcasHLA kallisto index..."
  ${ARCAS} arcasHLA reference --update
fi

# --- Per-sample loop ---------------------------------------------------------
while IFS= read -r SRR; do
  [[ -z "${SRR}" ]] && continue
  json="${OUT_DIR}/${SRR}/${SRR}.genotype.json"
  if [[ -f "${json}" ]]; then
    echo "SKIP ${SRR}: already genotyped (${json})"
    continue
  fi

  echo
  echo "========================================================"
  echo "[$(date '+%H:%M:%S')] Starting ${SRR}"
  echo "========================================================"

  # 1. Download SRA (resumable, dbGaP-authorized)
  echo "[1/5] prefetch ${SRR}"
  prefetch --ngc "${NGC}" -O "${SRA_DIR}" "${SRR}"

  # 2. Decode to paired FASTQ
  echo "[2/5] fasterq-dump ${SRR}"
  fasterq-dump --ngc "${NGC}" --split-files --threads "${THREADS}" \
    -O "${SCRATCH}/fq" -t "${SCRATCH}" \
    "${SRA_DIR}/${SRR}/${SRR}.sra" || \
    fasterq-dump --ngc "${NGC}" --split-files --threads "${THREADS}" \
      -O "${SCRATCH}/fq" -t "${SCRATCH}" "${SRR}"

  # 3. Align to GRCh38 (sort + index in one pipeline) using env's bowtie2
  # and samtools so versions match arcasHLA expectations
  echo "[3/5] bowtie2 -> samtools sort ${SRR}"
  ${ARCAS} bash -c "bowtie2 -x '${BT2_INDEX}' \
    -1 '${SCRATCH}/fq/${SRR}_1.fastq' \
    -2 '${SCRATCH}/fq/${SRR}_2.fastq' \
    -p ${THREADS} --very-sensitive 2>> '${SCRATCH}/${SRR}.bowtie2.log' \
    | samtools sort -@ 2 -o '${SCRATCH}/bam/${SRR}.bam' -"
  ${ARCAS} samtools index "${SCRATCH}/bam/${SRR}.bam"

  # Delete FASTQ now that BAM exists
  rm -f "${SCRATCH}/fq/${SRR}_1.fastq" "${SCRATCH}/fq/${SRR}_2.fastq"

  # 4. arcasHLA extract chr6 reads
  echo "[4/5] arcasHLA extract ${SRR}"
  mkdir -p "${OUT_DIR}/${SRR}"
  ${ARCAS} arcasHLA extract "${SCRATCH}/bam/${SRR}.bam" \
    -o "${OUT_DIR}/${SRR}" -t "${THREADS}" --paired

  # 5. arcasHLA genotype
  echo "[5/5] arcasHLA genotype ${SRR}"
  ${ARCAS} arcasHLA genotype \
    "${OUT_DIR}/${SRR}/${SRR}.extracted.1.fq.gz" \
    "${OUT_DIR}/${SRR}/${SRR}.extracted.2.fq.gz" \
    -g A,B,C,DPB1,DQB1,DRB1 \
    -o "${OUT_DIR}/${SRR}" -t "${THREADS}"

  # Cleanup: drop BAM + extracted fastqs (keep JSON)
  rm -f "${SCRATCH}/bam/${SRR}.bam" "${SCRATCH}/bam/${SRR}.bam.bai"
  rm -f "${OUT_DIR}/${SRR}"/*.extracted.*.fq.gz

  # Optionally drop the SRA to save external disk (comment out to cache)
  rm -rf "${SRA_DIR}/${SRR}"

  echo "[$(date '+%H:%M:%S')] Done ${SRR}. JSON: ${json}"
done < "${SRR_LIST}"

# --- Summary -----------------------------------------------------------------
echo
echo "========================================================"
echo "All samples processed. Typing results:"
ls -1 "${OUT_DIR}"/*/*.genotype.json 2>/dev/null | wc -l
echo "========================================================"
