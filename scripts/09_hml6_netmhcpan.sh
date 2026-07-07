#!/usr/bin/env bash
# NetMHCpan-4.1 binding predictions for HML6_20p11.21 ORFs.
# Assumes netMHCpan is in PATH and the 'tmp' directory is writable.
set -euo pipefail

INPUT=/Users/jezmarston/Desktop/Bulk_Acral_Claude/results/hml6_netmhcpan_input.fasta
OUT=/Users/jezmarston/Desktop/Bulk_Acral_Claude/results/hml6_netmhcpan_predictions.tsv
ALLELES=/Users/jezmarston/Desktop/Bulk_Acral_Claude/results/hla_unique_alleles_netmhcpan.txt

# Collapse the allele file into a comma-separated list
ALLELE_ARG=$(paste -sd, "${ALLELES}")

# 9-mer predictions with rank-based binder classification
netMHCpan -f "${INPUT}" \
          -BA \
          -l 9 \
          -a "${ALLELE_ARG}" \
          -xls \
          -xlsfile "${OUT}"

echo "NetMHCpan predictions written to ${OUT}"
echo "Strong binders: Rank < 0.5; Weak binders: Rank < 2.0."
