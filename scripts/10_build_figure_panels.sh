#!/usr/bin/env bash
# Build figure_panels_main.pdf, figure_panels_supplementary.pdf, and
# figure_panels_all.pdf by concatenating per-panel PDFs produced by the
# notebook sections. Re-run this after the notebook completes to refresh
# the submission-ready figure bundles.
#
# Usage:
#   bash scripts/10_build_figure_panels.sh

set -euo pipefail

FIGDIR=/Users/jezmarston/Desktop/Bulk_Acral_Claude/results/figures
OUTDIR=/Users/jezmarston/Desktop/Bulk_Acral_Claude/docs

require() {
  for f in "$@"; do
    if [[ ! -f "$f" ]]; then
      echo "WARN: missing $f; run the notebook to generate it" >&2
    fi
  done
}

# --- Main (one panel per page, manuscript order) -----------------------------
MAIN_PANELS=(
  # Figure 1 (prognostic signature)
  "${FIGDIR}/fig1A_univariate_cox_volcano.pdf"     # 1A
  "${FIGDIR}/fig1B_loocv_permutation.pdf"          # 1B
  "${FIGDIR}/fig1C_timedep_auc.pdf"                # 1C
  "${FIGDIR}/km_candidate_signature.pdf"           # 1D (KM tertile)
  "${FIGDIR}/mv_cox_forest.pdf"                    # 1E (forest)
  "${FIGDIR}/mv_cox_M3_ph_diagnostic.pdf"          # 1 supp (PH)

  # Figure 2 (immune + oncofetal)
  "${FIGDIR}/fig2A_immune_heatmap.pdf"             # 2A
  "${FIGDIR}/fig2B_nk_scatter.pdf"                 # 2B
  "${FIGDIR}/fig2C_oncofetal_bars.pdf"             # 2C
  "${FIGDIR}/oncofetal_axis_scatter.pdf"           # 2D

  # Figure 3 (ICI response)
  "${FIGDIR}/ici_response_core_boxplot.pdf"        # 3A
  "${FIGDIR}/ici_neoantigen_paradox.pdf"           # 3C

  # Figure 4 (HML6 antigen)
  "${FIGDIR}/fig4A_hml_tracks.pdf"                 # 4A
  "${FIGDIR}/fig4B_orf_motif_table.pdf"            # 4B
  "${FIGDIR}/hml_core_coexpression.pdf"            # 4C
  "${FIGDIR}/hml6_binder_coverage_per_patient.pdf" # 4D
)
require "${MAIN_PANELS[@]}"

# Build: only include files that actually exist (so it works incrementally)
existing_main=()
for f in "${MAIN_PANELS[@]}"; do
  [[ -f "$f" ]] && existing_main+=("$f")
done
pdfunite "${existing_main[@]}" "${OUTDIR}/figure_panels_main.pdf"
pages_main=$(pdfinfo "${OUTDIR}/figure_panels_main.pdf" | awk '/^Pages/ {print $2}')
echo "Built figure_panels_main.pdf: ${pages_main} pages, $(du -h "${OUTDIR}/figure_panels_main.pdf" | cut -f1)"

# --- Supplementary -----------------------------------------------------------
SUPP_PANELS=(
  "${FIGDIR}/pca_coloured_by_PATIENT_ID.pdf"
  "${FIGDIR}/supp_pc1_pc2_drivers.pdf"
  "${FIGDIR}/pca_HERV.pdf"
  "${FIGDIR}/pca_HERV-K.pdf"
  "${FIGDIR}/pca_L1.pdf"
  "${FIGDIR}/pca_gene.pdf"
  "${FIGDIR}/pc_covariate_diagnostic.pdf"
  "${FIGDIR}/stromal_by_cluster_k2.pdf"
  "${FIGDIR}/stromal_by_cluster_k4.pdf"
  "${FIGDIR}/km_cluster_data.driven.k2..Supercluster..pdf"
  "${FIGDIR}/km_cluster_data.driven.k4.pdf"
  "${FIGDIR}/km_cluster_original.manual.CLUSTER.pdf"
  "${FIGDIR}/purity_pathology_vs_estimate.pdf"
)

existing_supp=()
for f in "${SUPP_PANELS[@]}"; do
  [[ -f "$f" ]] && existing_supp+=("$f")
done
pdfunite "${existing_supp[@]}" "${OUTDIR}/figure_panels_supplementary.pdf"
pages_supp=$(pdfinfo "${OUTDIR}/figure_panels_supplementary.pdf" | awk '/^Pages/ {print $2}')
echo "Built figure_panels_supplementary.pdf: ${pages_supp} pages, $(du -h "${OUTDIR}/figure_panels_supplementary.pdf" | cut -f1)"

# --- Combined ----------------------------------------------------------------
pdfunite "${OUTDIR}/figure_panels_main.pdf" \
         "${OUTDIR}/figure_panels_supplementary.pdf" \
         "${OUTDIR}/figure_panels_all.pdf"
pages_all=$(pdfinfo "${OUTDIR}/figure_panels_all.pdf" | awk '/^Pages/ {print $2}')
echo "Built figure_panels_all.pdf: ${pages_all} pages, $(du -h "${OUTDIR}/figure_panels_all.pdf" | cut -f1)"
