# Locus-Specific Retrotransposable Element Expression in Acral Melanoma

Code accompanying:

> Marston JL, Fei T, Reyes-Gopar H, Bendall ML, Nixon DF (2026). **Locus-specific HERV
> expression identifies an aggressive, NK-depleted, checkpoint-refractory acral melanoma
> phenotype.** *Frontiers in Medicine* 13:1834744. https://doi.org/10.3389/fmed.2026.1834744

## Overview

This repository contains the analysis pipeline used to quantify locus-specific human
endogenous retrovirus (HERV) and LINE-1 expression in acral melanoma (AM) bulk RNA-seq,
using [Telescope](https://github.com/mlbendall/telescope) with a RepeatMasker-derived
annotation, and to relate locus-level expression to tumor purity, immune infiltrate,
HLA typing, and clinical outcomes (ICI response, overall survival).

```
FASTQ (paired-end)
  -> bowtie2 --very-sensitive -k 100 -p 4
BAM (multi-mapping, sorted, indexed)
  -> telescope assign --reassign_mode best_random
per-locus read counts
  -> edgeR / DESeq2
Differential loci (acral vs. non-acral vs. normal)
  -> Cox / logistic regression on clinical outcomes
Clinical-outcome-linked HERV-K loci
```

## Scripts

| Script | Purpose |
|---|---|
| `01_prepare_references.sh` | Build GRCh38 + RepeatMasker HML-2 reference/annotation |
| `02_align_telescope.sh` | bowtie2 alignment + Telescope locus-level quantification |
| `03_merge_telescope_counts.py` | Merge per-sample Telescope reports into a locus-by-sample count matrix |
| `04_kallisto_gene_counts.sh` | Gene-level quantification (kallisto) for cross-validation |
| `08_arcasHLA_all36.sh`, `12_arcasHLA_missing8.sh` | HLA typing (arcasHLA) from RNA-seq |
| `09_hml6_netmhcpan.sh`, `09b_hml6_mhcflurry.sh`, `13_hervh_mhcflurry.sh` | Peptide-HLA binding prediction for candidate HERV-derived antigens |
| `10_build_figure_panels.sh`, `11_build_latex_figures.sh` | Manuscript figure generation |
| `14_build_supplementary_tables.py` | Supplementary table assembly |
| `gse189889_slurm.sh` | Slurm batch submission for the GSE189889 cohort |

## Dependencies

- [Telescope](https://github.com/mlbendall/telescope) (build from source; the PyPI
  `telescope-ngs` package is broken)
- bowtie2, samtools
- arcasHLA
- kallisto
- R (edgeR/DESeq2) for differential expression and survival analysis

## Data availability

Sample-level RNA-seq data used in this study are from the publicly available
NCBI SRA BioProject **PRJNA304068** (33 patients, 36 samples), with clinical
metadata from Liang et al. (Supplementary Tables in *Genome Research* 2017;
27:24-32). See the manuscript's Data Availability Statement and Supplementary
Material for full details.

## Citation

If you use this pipeline, please cite the paper above.
