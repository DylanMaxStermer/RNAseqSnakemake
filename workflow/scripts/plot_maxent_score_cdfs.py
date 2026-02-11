#!/usr/bin/env python3
"""
Plot empirical CDFs of MaxEntScan splice site scores by exon type.

Generates two plots:
  - donor_maxent_cdf.png  : donor (5' splice site) scores
  - acceptor_maxent_cdf.png : acceptor (3' splice site) scores

Each plot has one CDF line per exon-type category.
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# ---------------------------------------------------------------------------
# Input file definitions
# ---------------------------------------------------------------------------

SPLICE_SITES_BASE = (
    "/project/yangili1/dylan_stermer/GitHubTesting/RNAseq_Snakemake/"
    "results/ExonCharacteristics/GRCh38_GencodeRelease44Comprehensive/SpliceSites"
)
LNCRNA_BASE = (
    "/project/yangili1/dylan_stermer/GitHubTesting//"
    "maxEnt/output/scores"
)

CATEGORIES = [
    {
        "label": "Constitutive",
        "donor":    f"{SPLICE_SITES_BASE}/Constitutive/donor_sites_scored.bed",
        "acceptor": f"{SPLICE_SITES_BASE}/Constitutive/acceptor_sites_scored.bed",
    },
    {
        "label": "Exon flanking PE",
        "donor":    f"{SPLICE_SITES_BASE}/exon_flanking_PE_all/donor_sites_scored.bed",
        "acceptor": f"{SPLICE_SITES_BASE}/exon_flanking_PE_all/acceptor_sites_scored.bed",
    },
    {
        "label": "Leafcutter CE (GTF matched)",
        "donor":    f"{SPLICE_SITES_BASE}/leafcutter_CE_gtf_matched/donor_sites_scored.bed",
        "acceptor": f"{SPLICE_SITES_BASE}/leafcutter_CE_gtf_matched/acceptor_sites_scored.bed",
    },
    {
        "label": "Poison exon",
        "donor":    f"{SPLICE_SITES_BASE}/poisonExon_all/donor_sites_scored.bed",
        "acceptor": f"{SPLICE_SITES_BASE}/poisonExon_all/acceptor_sites_scored.bed",
    },
    {
        "label": "lncRNA",
        "donor":    f"{LNCRNA_BASE}/lncRNA_donor_sites_scored.bed",
        "acceptor": f"{LNCRNA_BASE}/lncRNA_acceptor_sites_scored.bed",
    },
]

OUTPUT_DIR = (
    "/project/yangili1/dylan_stermer/GitHubTesting/RNAseq_Snakemake/analysis"
)

# Column index of maxent_score (0-based)
SCORE_COL = 7

# Distinct colors for up to 6 categories
COLORS = ["#2196F3", "#E91E63", "#4CAF50", "#FF9800", "#9C27B0", "#00BCD4"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_scores(filepath: str) -> np.ndarray:
    """Load maxent_score column from a scored BED file, skipping comment lines."""
    if not os.path.exists(filepath):
        print(f"  WARNING: file not found — {filepath}", file=sys.stderr)
        return np.array([])
    df = pd.read_csv(filepath, sep="\t", comment="#", header=None)
    scores = pd.to_numeric(df.iloc[:, SCORE_COL], errors="coerce").dropna().values
    return scores


def empirical_cdf(scores: np.ndarray):
    """Return (x, y) arrays for an empirical CDF."""
    if len(scores) == 0:
        return np.array([]), np.array([])
    x = np.sort(scores)
    y = np.arange(1, len(x) + 1) / len(x)
    return x, y


def make_cdf_plot(site_type: str, output_path: str):
    """Create and save a CDF plot for donor or acceptor sites."""
    fig, ax = plt.subplots(figsize=(9, 5))

    for i, cat in enumerate(CATEGORIES):
        filepath = cat[site_type]
        scores = load_scores(filepath)
        n = len(scores)
        if n == 0:
            print(f"  Skipping {cat['label']} — no data")
            continue
        x, y = empirical_cdf(scores)
        color = COLORS[i % len(COLORS)]
        ax.plot(x, y, color=color, linewidth=1.8,
                label=f"{cat['label']} (n={n:,})")
        print(f"  {cat['label']}: n={n:,}, median={np.median(scores):.2f}")

    site_label = "Donor (5' splice site)" if site_type == "donor" else "Acceptor (3' splice site)"
    ax.set_xlabel("MaxEntScan score", fontsize=12)
    ax.set_ylabel("CDF", fontsize=12)
    ax.set_title(f"MaxEntScan score distribution — {site_label}", fontsize=13)
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 1)

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {output_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("=== Donor CDF ===")
    make_cdf_plot("donor", os.path.join(OUTPUT_DIR, "donor_maxent_cdf.png"))

    print("\n=== Acceptor CDF ===")
    make_cdf_plot("acceptor", os.path.join(OUTPUT_DIR, "acceptor_maxent_cdf.png"))

    print("\nDone.")


if __name__ == "__main__":
    main()
