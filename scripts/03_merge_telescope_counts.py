#!/usr/bin/env python3
"""
03_merge_telescope_counts.py

Merges per-sample Telescope reports into a single locus × sample count matrix.

Input:
    results/telescope/{SRR}/{SRR}_telescope_report.tsv  (one per sample)

Output:
    Acral R/HERV_hiconf_AM.csv   — all Telescope loci, final_count column
    results/telescope_merged_full.tsv  — same, tab-separated, all Telescope columns

Filtering (high-confidence loci written to HERV_hiconf_AM.csv):
    Loci with final_count > 0 in at least MIN_SAMPLES samples are retained.
    This matches the filtering applied in the paper (no floor; all loci included
    for DESeq2 independent filtering downstream).

Usage:
    python3 scripts/03_merge_telescope_counts.py
    python3 scripts/03_merge_telescope_counts.py --min_samples 1  # include all loci
"""

import argparse, glob, os, re, sys
import pandas as pd

# ── Configuration ──────────────────────────────────────────────────────────────
PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

TELESCOPE_DIR = os.path.join(PROJECT, "results", "telescope")
OUT_CSV       = os.path.join(PROJECT, "Acral R", "HERV_hiconf_AM.csv")
OUT_TSV       = os.path.join(PROJECT, "results", "telescope_merged_full.tsv")

# Sample order used in the paper (matches Acral_Melanoma_Metadata_Cluster_HLA.txt)
SAMPLE_ORDER = [
    "SRR3882763", "SRR3882768", "SRR3882773", "SRR3882778", "SRR3882783",
    "SRR3882788", "SRR3882793", "SRR3882800", "SRR3882805", "SRR3882810",
    "SRR3882815", "SRR3882820", "SRR3882825", "SRR3882830", "SRR3882835",
    "SRR3882840", "SRR3882845", "SRR3882849", "SRR3882853", "SRR3882857",
    "SRR3882861", "SRR3882864", "SRR3882868", "SRR3882873", "SRR3882876",
    "SRR3882880", "SRR3882883", "SRR3882888", "SRR3882892", "SRR3882897",
    "SRR3882901", "SRR3882905", "SRR3882910", "SRR3882913", "SRR3882918",
    "SRR3882920",
]


def load_report(path: str, srr: str) -> pd.Series:
    """Load final_count column from a Telescope report TSV."""
    df = pd.read_csv(path, sep="\t", comment="#", index_col=0)
    if "final_count" not in df.columns:
        raise ValueError(f"No 'final_count' column in {path}")
    # Drop the '__no_feature' and '__ambiguous' pseudo-rows if present
    df = df[~df.index.str.startswith("__")]
    return df["final_count"].rename(srr)


def main(min_samples: int = 1) -> None:
    # ── Collect per-sample reports ────────────────────────────────────────────
    reports = {}
    missing = []
    for srr in SAMPLE_ORDER:
        pattern = os.path.join(TELESCOPE_DIR, srr, f"{srr}_telescope_report.tsv")
        matches = glob.glob(pattern)
        if not matches:
            missing.append(srr)
        else:
            reports[srr] = load_report(matches[0], srr)

    if missing:
        print(f"WARNING: Missing Telescope reports for {len(missing)} samples: {missing}",
              file=sys.stderr)

    if not reports:
        sys.exit("ERROR: No Telescope reports found. "
                 "Run 02_align_telescope.sh first.")

    # ── Build count matrix ─────────────────────────────────────────────────────
    found_order = [s for s in SAMPLE_ORDER if s in reports]
    mat = pd.DataFrame({s: reports[s] for s in found_order}).fillna(0).astype(int)
    mat.index.name = ""
    print(f"Loaded {mat.shape[0]} loci × {mat.shape[1]} samples")

    # ── Filter: loci with any reads ────────────────────────────────────────────
    detected = (mat > 0).sum(axis=1)
    mat = mat[detected >= min_samples]
    print(f"After filtering (≥{min_samples} samples with count > 0): {mat.shape[0]} loci")

    # ── Write outputs ──────────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    mat.to_csv(OUT_CSV)
    print(f"Wrote {OUT_CSV}")

    mat.to_csv(OUT_TSV, sep="\t")
    print(f"Wrote {OUT_TSV}")

    # Summary stats
    total_counts = mat.sum(axis=0)
    print("\nPer-sample total HERV counts (should be ≥ 1 000 for usable samples):")
    for s, c in total_counts.items():
        flag = " ← LOW" if c < 1000 else ""
        print(f"  {s}: {c:,}{flag}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--min_samples", type=int, default=1,
                        help="Minimum number of samples with count > 0 to retain a locus")
    args = parser.parse_args()
    main(min_samples=args.min_samples)
