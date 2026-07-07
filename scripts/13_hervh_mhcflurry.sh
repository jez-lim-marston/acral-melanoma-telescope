#!/usr/bin/env bash
# MHCflurry 2.2 binding predictions for HERVH_6q21a ORFs.
# Mirrors scripts/09b_hml6_mhcflurry.sh; see that script for env setup notes.
#
# Run:
#   cd /Users/jezmarston/Desktop/Bulk_Acral_Claude
#   bash scripts/13_hervh_mhcflurry.sh

set -euo pipefail

PROJECT=/Users/jezmarston/Desktop/Bulk_Acral_Claude
OUTDIR=${PROJECT}/results
INPUT=${OUTDIR}/hervh_netmhcpan_input.fasta
ALLELES_FILE=${OUTDIR}/hla_unique_alleles_netmhcpan.txt
ENV_NAME=mhcflurry

if ! mamba env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  echo "ERROR: mamba env '${ENV_NAME}' not found." >&2
  echo "  See scripts/09b_hml6_mhcflurry.sh for env setup." >&2
  exit 1
fi
if ! mamba run -n "${ENV_NAME}" python -c "import mhcflurry" 2>/dev/null; then
  echo "ERROR: mhcflurry not importable in env '${ENV_NAME}'." >&2; exit 1
fi
if [[ ! -f "${INPUT}" ]]; then
  echo "ERROR: ${INPUT} not found - run §12c2 in R first." >&2; exit 1
fi
if [[ ! -f "${ALLELES_FILE}" ]]; then
  echo "ERROR: ${ALLELES_FILE} not found - run §12c in R first." >&2; exit 1
fi

# Reformat NetMHCpan allele syntax -> MHCflurry syntax:
#   HLA-A02:01  ->  HLA-A*02:01
ALLELES=$(paste -sd, "${ALLELES_FILE}" | \
  sed -E 's/(HLA-[ABC])([0-9]{2}:[0-9]{2})/\1\*\2/g')

echo "Alleles (MHCflurry format): ${ALLELES}"
echo "Input FASTA: ${INPUT}"

PEPTIDES=${OUTDIR}/hervh_peptides_alleles.csv
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
print(f"Wrote {len(peptides)} peptides x {len(alleles)} alleles = "
      f"{len(peptides) * len(alleles)} rows to {out_file}")
PYEOF

OUT=${OUTDIR}/hervh_mhcflurry_predictions.csv
mamba run -n "${ENV_NAME}" mhcflurry-predict \
  "${PEPTIDES}" \
  --out "${OUT}" \
  --no-flanking

echo
echo "Predictions: ${OUT}"
echo

mamba run -n "${ENV_NAME}" python - "${OUT}" <<'PYEOF'
import sys, pandas as pd
df = pd.read_csv(sys.argv[1])
print("Total predictions:", len(df))
print("Strong binders (percentile_rank < 0.5):",
      (df['mhcflurry_affinity_percentile'] < 0.5).sum())
print("Weak binders    (percentile_rank < 2.0):",
      (df['mhcflurry_affinity_percentile'] < 2.0).sum())
print()
print("HLA-A*02:01 strong binders:")
a0201 = df[(df['allele'] == 'HLA-A*02:01') & (df['mhcflurry_affinity_percentile'] < 0.5)]
print(f"  {len(a0201)} peptides")
if len(a0201):
    print(a0201[['peptide','mhcflurry_affinity_percentile']].head(10).to_string())
PYEOF
