#!/usr/bin/env bash
# 01_prepare_references.sh
#
# Downloads and prepares all reference files needed for the HERV/TE and gene
# quantification pipelines:
#   • GRCh38 primary assembly (GENCODE v45)
#   • Bowtie2 index (for Telescope multi-mapping alignment)
#   • Telescope annotation GTF (RepeatMasker UCSC → ERVK-filtered)
#   • Kallisto transcript index (for gene-level expression)
#
# Run once before starting per-sample processing.
# Estimated wall time: 3–5 h (index builds dominate; parallelised where possible)
#
# Usage:
#   bash scripts/01_prepare_references.sh
#
# Requirements (PATH):
#   bowtie2-build  (≥ 2.4)
#   kallisto       (≥ 0.46)
#   samtools       (≥ 1.16)
#   wget / curl

set -euo pipefail

PROJECT="$(cd "$(dirname "$0")/.." && pwd)"
REF="${PROJECT}/data/reference"
LOGDIR="${PROJECT}/logs"
mkdir -p "${REF}" "${LOGDIR}"

log() { echo "[$(date '+%H:%M:%S')] $*" | tee -a "${LOGDIR}/01_prepare_references.log"; }

# ── A. GRCh38 primary assembly ────────────────────────────────────────────────
GENOME_GZ="${REF}/GRCh38.primary_assembly.genome.fa.gz"
GENOME_FA="${REF}/GRCh38.primary_assembly.genome.fa"

if [[ ! -f "${GENOME_FA}" ]]; then
  log "Downloading GRCh38 primary assembly (GENCODE v45) ..."
  caffeinate -i -s curl -L -C - --max-time 14400 \
    -o "${GENOME_GZ}" \
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/GRCh38.primary_assembly.genome.fa.gz"
  log "Decompressing ..."
  gunzip -k "${GENOME_GZ}"
  log "Genome ready: ${GENOME_FA}"
else
  log "Genome already present; skipping download."
fi

# ── B. Bowtie2 index (for Telescope / HERV multi-mapping) ────────────────────
BT2_PREFIX="${REF}/bowtie2_index/GRCh38"
BT2_DONE="${REF}/bowtie2_index/.done"

if [[ ! -f "${BT2_DONE}" ]]; then
  mkdir -p "$(dirname "${BT2_PREFIX}")"
  log "Building bowtie2 index (4 threads; ~2–3 h) ..."
  caffeinate -i -s bowtie2-build \
    --threads 4 \
    "${GENOME_FA}" \
    "${BT2_PREFIX}" \
    2>&1 | tee "${LOGDIR}/bt2_build.log"
  touch "${BT2_DONE}"
  log "Bowtie2 index ready: ${BT2_PREFIX}"
else
  log "Bowtie2 index already built; skipping."
fi

# ── C. Telescope annotation GTF (ERVK-only, from UCSC RepeatMasker) ──────────
# The Bendall lab retro.hg38.v1.gtf (curated, named loci) is preferred where
# available.  Fall back to constructing one from UCSC rmsk.txt.gz.
TELESCOPE_GTF="${REF}/HERVK_telescope.gtf"

if [[ ! -f "${TELESCOPE_GTF}" ]]; then
  log "Downloading UCSC RepeatMasker annotation ..."
  RMSK_GZ="${REF}/rmsk.txt.gz"
  caffeinate -i -s curl -L -C - --max-time 7200 \
    -o "${RMSK_GZ}" \
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz"

  log "Filtering ERVK loci and converting to Telescope GTF ..."
  python3 - <<'PYEOF'
import gzip, re, sys

RMSK_GZ   = "data/reference/rmsk.txt.gz"
OUT_GTF   = "data/reference/HERVK_telescope.gtf"

# UCSC rmsk.txt column layout (0-based):
# 0  bin, 1 swScore, 2 milliDiv, 3 milliDel, 4 milliIns,
# 5  genoName, 6 genoStart (0-based), 7 genoEnd, 8 strand,
# 9  repName, 10 repClass, 11 repFamily, …

kept = 0
with gzip.open(RMSK_GZ, "rt") as fi, open(OUT_GTF, "w") as fo:
    for line in fi:
        cols = line.rstrip().split("\t")
        if len(cols) < 13:
            continue
        rep_class  = cols[11]
        rep_family = cols[12] if len(cols) > 12 else ""
        if rep_class != "LTR" or rep_family != "ERVK":
            continue
        chrom   = cols[5]
        start   = int(cols[6])          # 0-based
        end     = int(cols[7])          # half-open
        strand  = cols[8] if cols[8] in ("+", "-") else "+"
        rep_name = cols[9]
        locus_id = f"{rep_name}_{chrom}_{start}_{end}"
        attrs = (
            f'gene_id "{locus_id}"; '
            f'transcript_id "{locus_id}"; '
            f'family_id "ERVK"; '
            f'class_id "LTR"; '
            f'repName "{rep_name}";'
        )
        gtf_start = start + 1   # GTF is 1-based
        fo.write(
            f"{chrom}\tRepeatMasker\texon\t{gtf_start}\t{end}\t.\t{strand}\t.\t{attrs}\n"
        )
        kept += 1

print(f"Wrote {kept} ERVK loci to {OUT_GTF}")
PYEOF

  log "Telescope GTF ready: ${TELESCOPE_GTF}"
else
  log "Telescope GTF already present; skipping."
fi

# ── D. Kallisto transcript index (for gene-level expression) ──────────────────
# GENCODE v45 transcriptome (protein-coding + lncRNA)
TRANSCRIPTOME_GZ="${REF}/gencode.v45.transcripts.fa.gz"
KALLISTO_INDEX="${REF}/kallisto_index.idx"

if [[ ! -f "${KALLISTO_INDEX}" ]]; then
  if [[ ! -f "${TRANSCRIPTOME_GZ}" ]]; then
    log "Downloading GENCODE v45 transcript FASTA ..."
    caffeinate -i -s curl -L -C - --max-time 7200 \
      -o "${TRANSCRIPTOME_GZ}" \
      "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.transcripts.fa.gz"
  fi
  log "Building kallisto index ..."
  caffeinate -i -s kallisto index \
    -i "${KALLISTO_INDEX}" \
    "${TRANSCRIPTOME_GZ}" \
    2>&1 | tee "${LOGDIR}/kallisto_index.log"
  log "Kallisto index ready: ${KALLISTO_INDEX}"
else
  log "Kallisto index already present; skipping."
fi

log "All reference files are ready."
log "  Genome:            ${GENOME_FA}"
log "  Bowtie2 index:     ${BT2_PREFIX}.*.bt2"
log "  Telescope GTF:     ${TELESCOPE_GTF}"
log "  Kallisto index:    ${KALLISTO_INDEX}"
