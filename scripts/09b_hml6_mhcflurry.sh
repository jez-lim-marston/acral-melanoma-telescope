#!/usr/bin/env bash
# MHCflurry-2.0 binding predictions for HML6_20p11.21 ORFs.
# Runs entirely inside the dedicated "mhcflurry" mamba env so it doesn't
# matter what Python is active in the calling shell.
#
# One-time setup:
#   mamba create -y -n mhcflurry python=3.11 pip biopython
#   mamba run -n mhcflurry pip install mhcflurry
#   mamba run -n mhcflurry mhcflurry-downloads fetch
#
# Run:
#   cd /Users/jezmarston/Desktop/Bulk_Acral_Claude
#   bash scripts/09b_hml6_mhcflurry.sh

set -euo pipefail

PROJECT=/Users/jezmarston/Desktop/Bulk_Acral_Claude
OUTDIR=${PROJECT}/results
INPUT=${OUTDIR}/hml6_netmhcpan_input.fasta
ALLELES_FILE=${OUTDIR}/hla_unique_alleles_netmhcpan.txt
ENV_NAME=mhcflurry

# Pre-flight: is the env actually present and does it have mhcflurry?
if ! mamba env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  echo "ERROR: mamba env '${ENV_NAME}' not found. Create it with:" >&2
  echo "  mamba create -y -n ${ENV_NAME} python=3.11 pip biopython" >&2
  echo "  mamba run -n ${ENV_NAME} pip install mhcflurry" >&2
  echo "  mamba run -n ${ENV_NAME} mhcflurry-downloads fetch" >&2
  exit 1
fi
if ! mamba run -n "${ENV_NAME}" python -c "import mhcflurry" 2>/dev/null; then
  echo "ERROR: mhcflurry not importable in env '${ENV_NAME}'. Install with:" >&2
  echo "  mamba run -n ${ENV_NAME} pip install mhcflurry" >&2
  echo "  mamba run -n ${ENV_NAME} mhcflurry-downloads fetch" >&2
  exit 1
fi

if [[ ! -f "${INPUT}" ]]; then
  echo "ERROR: ${INPUT} not found — run §12c in R first." >&2; exit 1
fi
if [[ ! -f "${ALLELES_FILE}" ]]; then
  echo "ERROR: ${ALLELES_FILE} not found — run §12c in R first." >&2; exit 1
fi

# Reformat NetMHCpan allele syntax → MHCflurry syntax:
#   HLA-A02:01  →  HLA-A*02:01
ALLELES=$(paste -sd, "${ALLELES_FILE}" | \
  sed -E 's/(HLA-[ABC])([0-9]{2}:[0-9]{2})/\1\*\2/g')

echo "Alleles (MHCflurry format): ${ALLELES}"
echo "Input FASTA: ${INPUT}"

# Tile the FASTA into 8–11 mer peptides and write a peptide × allele
# Cartesian-product CSV. MHCflurry 2.2 requires all allele info to be IN
# the input CSV when --peptides/--alleles flags are omitted.
PEPTIDES=${OUTDIR}/hml6_peptides_alleles.csv
mamba run -n "${ENV_NAME}" python - "${INPUT}" "${PEPTIDES}" "${ALLELES}" <<'PYEOF'
import sys
from Bio import SeqIO
input_fa, out_file, allele_csv = sys.argv[1], sys.argv[2], sys.argv[3]
alleles = allele_csv.split(",")

peptides = set()
for rec in SeqIO.parse(input_fa, "fasta"):
    seq = str(rec.seq)
    for k in (8, 9, 10, 11):
        for i in range(len(seq) - k + 1):
            p = seq[i:i+k]
            if "*" in p or "X" in p: continue
            peptides.add(p)

peptides = sorted(peptides)
with open(out_file, "w") as fh:
    fh.write("peptide,allele\n")
    for a in alleles:
        for p in peptides:
            fh.write(f"{p},{a}\n")
print(f"Wrote {len(peptides)} peptides × {len(alleles)} alleles = "
      f"{len(peptides) * len(alleles)} rows to {out_file}")
PYEOF

# Run MHCflurry inside the env (input CSV carries peptide+allele columns).
OUT=${OUTDIR}/hml6_mhcflurry_predictions.csv
mamba run -n "${ENV_NAME}" mhcflurry-predict \
  "${PEPTIDES}" \
  --out "${OUT}" \
  --no-flanking

echo ""
echo "Predictions: ${OUT}"
echo ""

# Binder summary (env-isolated)
mamba run -n "${ENV_NAME}" python - "${OUT}" <<'PYEOF'
import sys, pandas as pd
df = pd.read_csv(sys.argv[1])
print("Total predictions:", len(df))
print("Strong binders (percentile_rank < 0.5):",
      (df['mhcflurry_affinity_percentile'] < 0.5).sum())
print("Weak binders    (percentile_rank < 2.0):",
      (df['mhcflurry_affinity_percentile'] < 2.0).sum())
print()
print("Binder distribution by allele (rank < 2%):")
print(df[df['mhcflurry_affinity_percentile'] < 2.0]
        .groupby('allele').size()
        .sort_values(ascending=False)
        .to_string())
PYEOF
