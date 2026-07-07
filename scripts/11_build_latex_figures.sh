#!/usr/bin/env bash
# Build publication-quality composite figures with LaTeX (subcaption +
# standalone). Reads per-panel PDFs from results/figures/ and writes
# Figure_*.pdf to docs/figures_latex/.
#
# Requirements:
#   - pdflatex on PATH (MacTeX, BasicTeX, or tinytex)
#   - TeX packages: standalone, subcaption, caption, graphicx, helvet
#   - With tinytex these auto-install on first compile
#
# Usage:
#   bash scripts/11_build_latex_figures.sh               # all figures
#   bash scripts/11_build_latex_figures.sh Figure_1      # single figure

set -euo pipefail

TEXDIR=/Users/jezmarston/Desktop/Bulk_Acral_Claude/docs/figures_latex
cd "${TEXDIR}"

if ! command -v pdflatex >/dev/null 2>&1; then
  # Search known install locations in priority order
  for candidate in \
    "$HOME/Library/TinyTeX/bin/universal-darwin/pdflatex" \
    "$HOME/Library/TinyTeX/bin/x86_64-darwin/pdflatex" \
    "/usr/local/texlive/2026basic/bin/universal-darwin/pdflatex" \
    "/usr/local/texlive/2025/bin/universal-darwin/pdflatex" \
    "/Library/TeX/texbin/pdflatex"; do
    if [[ -x "$candidate" ]]; then
      export PATH="$(dirname "$candidate"):$PATH"
      break
    fi
  done
  if ! command -v pdflatex >/dev/null 2>&1; then
    echo "ERROR: pdflatex not on PATH and no known installation found." >&2
    echo "Install with R:  install.packages('tinytex'); tinytex::install_tinytex()" >&2
    exit 1
  fi
fi

# Default to all Figure_*.tex; allow single-figure arg (e.g. "Figure_1")
if [[ $# -gt 0 ]]; then
  TARGETS=("$@")
else
  TARGETS=(Figure_1 Figure_2 Figure_3 Figure_4 Figure_5 Figure_S1 Figure_S2 Figure_S3 Figure_S4)
fi

echo "Building $(( ${#TARGETS[@]} )) figure(s) in ${TEXDIR}"

for name in "${TARGETS[@]}"; do
  stem="${name%.tex}"
  tex="${stem}.tex"
  if [[ ! -f "${tex}" ]]; then
    echo "  SKIP ${stem}: ${tex} not found" >&2
    continue
  fi
  echo "  --> pdflatex ${tex}"
  # -halt-on-error: fail fast; -interaction=nonstopmode: no prompts
  pdflatex -halt-on-error -interaction=nonstopmode "${tex}" \
    > "${stem}.build.log" 2>&1 || {
      echo "  FAILED: ${stem}. See ${TEXDIR}/${stem}.build.log" >&2
      tail -30 "${stem}.build.log" >&2
      exit 1
    }
  rm -f "${stem}.aux" "${stem}.log"
  size=$(du -h "${stem}.pdf" | cut -f1)
  echo "      wrote ${stem}.pdf (${size})"
done

echo
echo "All figures built in: ${TEXDIR}"
ls -lh "${TEXDIR}"/Figure_*.pdf
