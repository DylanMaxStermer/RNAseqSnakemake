#!/usr/bin/env python3
"""
Create binned CDF plots for splice site scores.

Bins events by cryptic_length magnitude and visualizes score distributions
as Cumulative Distribution Functions comparing canonical vs cryptic sites.
"""

import argparse
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def extract_9mer(sequence):
    """
    Extract 9-mer from exon|intron sequence.

    Rule: 3 exon bases + 6 intron bases

    Example: "GTTCCCACAG|GTATTCACAT" → "CAGGTATTC"
    - Exon part: "GTTCCCACAG" → last 3 bases: "CAG"
    - Intron part: "GTATTCACAT" → first 6 bases: "GTATTC"
    - Result: "CAG" + "GTATTC" = "CAGGTATTC"

    Args:
        sequence: String with format "EXONBASES|INTRONBASES"

    Returns:
        9-mer string or None
    """
    if pd.isna(sequence) or '|' not in sequence:
        return None

    exon_part, intron_part = sequence.split('|')

    # Extract last 3 bases from exon and first 6 bases from intron
    if len(exon_part) >= 3 and len(intron_part) >= 6:
        return exon_part[-3:] + intron_part[:6]
    else:
        return None


def load_and_process_data(analysis_file, reference_df):
    """
    Load analysis file and merge with reference scores.

    Args:
        analysis_file: Path to splice_site_analysis TSV
        reference_df: DataFrame with reference scores

    Returns:
        DataFrame with matched scores for canonical and cryptic sites
    """
    print(f"  Loading analysis file: {analysis_file}")

    # Read analysis file
    df = pd.read_csv(analysis_file, sep='\t')
    print(f"    Events: {len(df)}")

    # Extract 9-mers
    print("  Extracting 9-mers...")
    df['canonical_9mer'] = df['canonical_sequence'].apply(extract_9mer)
    df['cryptic_9mer'] = df['cryptic_sequence'].apply(extract_9mer)

    # Count successful extractions
    n_canonical = df['canonical_9mer'].notna().sum()
    n_cryptic = df['cryptic_9mer'].notna().sum()
    print(f"    Canonical 9-mers extracted: {n_canonical}/{len(df)}")
    print(f"    Cryptic 9-mers extracted: {n_cryptic}/{len(df)}")

    # Merge with reference scores for canonical
    print("  Matching canonical 9-mers to reference...")
    df = df.merge(
        reference_df[['SE.5SS.bed.seq', 'Balance', 'SE.5SS.bed.score', 'LUC7_score']],
        left_on='canonical_9mer',
        right_on='SE.5SS.bed.seq',
        how='left',
        suffixes=('', '_ref')
    )

    # Rename columns for clarity
    df = df.rename(columns={
        'Balance': 'canonical_Balance',
        'SE.5SS.bed.score': 'canonical_score',
        'LUC7_score': 'canonical_LUC7'
    })

    # Drop the extra SE.5SS.bed.seq column from merge
    df = df.drop(columns=['SE.5SS.bed.seq'], errors='ignore')

    # Merge with reference scores for cryptic
    print("  Matching cryptic 9-mers to reference...")
    df = df.merge(
        reference_df[['SE.5SS.bed.seq', 'Balance', 'SE.5SS.bed.score', 'LUC7_score']],
        left_on='cryptic_9mer',
        right_on='SE.5SS.bed.seq',
        how='left',
        suffixes=('', '_ref2')
    )

    # Rename columns for clarity
    df = df.rename(columns={
        'Balance': 'cryptic_Balance',
        'SE.5SS.bed.score': 'cryptic_score',
        'LUC7_score': 'cryptic_LUC7'
    })

    # Drop the extra SE.5SS.bed.seq column from merge
    df = df.drop(columns=['SE.5SS.bed.seq'], errors='ignore')

    # Count matches
    n_canonical_matched = df['canonical_score'].notna().sum()
    n_cryptic_matched = df['cryptic_score'].notna().sum()
    print(f"    Canonical matched: {n_canonical_matched}/{len(df)}")
    print(f"    Cryptic matched: {n_cryptic_matched}/{len(df)}")

    return df


def bin_data_by_cryptic_length(df, n_bins=4, use_quantile=True):
    """
    Bin data by cryptic_length or rank.

    Args:
        df: DataFrame with cryptic_length and rank columns
        n_bins: Number of bins (default 4)
        use_quantile: If True, use quantile-based binning (equal counts per bin).
                      If False, use rank-based binning (equal rank ranges).

    Returns:
        df_with_bins: DataFrame with added 'bin' column
        bin_info: List of dicts with bin metadata
    """
    # Sort by cryptic_length
    df = df.sort_values('cryptic_length').reset_index(drop=True)

    if use_quantile:
        # Quantile-based binning: equal counts per bin
        print(f"  Using quantile-based binning (equal counts per bin)")
        df['bin'] = pd.qcut(df['cryptic_length'], q=n_bins, labels=False, duplicates='raise')
    else:
        # Rank-based binning: divide rank into equal ranges
        print(f"  Using rank-based binning (equal rank ranges)")
        df = df.sort_values('rank').reset_index(drop=True)
        df['bin'] = pd.cut(df['rank'], bins=n_bins, labels=False, include_lowest=True)

    # Calculate bin metadata
    bin_info = []
    for bin_num in range(n_bins):
        bin_df = df[df['bin'] == bin_num]
        if len(bin_df) > 0:
            bin_info.append({
                'bin': bin_num,
                'min': int(bin_df['cryptic_length'].min()),
                'max': int(bin_df['cryptic_length'].max()),
                'count': len(bin_df),
                'label': f"Bin {bin_num + 1} ({int(bin_df['cryptic_length'].min())} to {int(bin_df['cryptic_length'].max())} nt)"
            })

    return df, bin_info


def compute_cdf(data):
    """
    Compute empirical CDF for a 1D array.

    Args:
        data: 1D numpy array or pandas Series

    Returns:
        x_values: Sorted unique values
        y_values: CDF values (0 to 1)
    """
    data_clean = data.dropna().values
    if len(data_clean) == 0:
        return np.array([]), np.array([])

    data_sorted = np.sort(data_clean)
    y = np.arange(1, len(data_sorted) + 1) / len(data_sorted)
    return data_sorted, y


def create_cdf_plots(df_longer, df_shorter, bin_info_longer, bin_info_shorter,
                     output_path, n_bins=4):
    """
    Create 2x3 grid of CDF plots.

    Layout:
    - Columns: Shorter (cryptic_length < 0) | Longer (cryptic_length > 0)
    - Rows: Balance | LUC7 | Splice Score

    For each subplot:
    - Plot CDFs for each bin
    - Canonical = solid line, Cryptic = dashed line
    - Each bin gets unique color (viridis colormap)
    - Same bin shares color for canonical and cryptic

    Args:
        df_longer: DataFrame for lengthening events (binned)
        df_shorter: DataFrame for shortening events (binned)
        bin_info_longer: Bin metadata for longer
        bin_info_shorter: Bin metadata for shorter
        output_path: Save path for figure
        n_bins: Number of bins
    """
    print(f"  Creating CDF plots: {output_path}")

    fig, axes = plt.subplots(3, 2, figsize=(16, 12))

    # Get color palette - blue, red, green, purple
    color_palette = ['#3498db', '#e74c3c', '#2ecc71', '#9b59b6']  # blue, red, green, purple
    colors = color_palette[:n_bins]

    # Row 0: Balance
    # Row 1: LUC7
    # Row 2: Splice Score

    metrics = [
        ('Balance', 'canonical_Balance', 'cryptic_Balance'),
        ('LUC7 Score', 'canonical_LUC7', 'cryptic_LUC7'),
        ('Splice Score', 'canonical_score', 'cryptic_score')
    ]

    for row_idx, (metric_name, canonical_col, cryptic_col) in enumerate(metrics):
        # Left column: Shorter events (cryptic_length < 0)
        ax_shorter = axes[row_idx, 0]
        for bin_num in range(n_bins):
            bin_data = df_shorter[df_shorter['bin'] == bin_num]
            if len(bin_data) > 0:
                # Get bin info for label
                bin_info = [b for b in bin_info_shorter if b['bin'] == bin_num][0]
                bin_label = f"Bin {bin_num+1} ({bin_info['min']} to {bin_info['max']} nt, n={bin_info['count']})"

                # Canonical - solid line
                x_can, y_can = compute_cdf(bin_data[canonical_col])
                if len(x_can) > 0:
                    ax_shorter.plot(x_can, y_can, color=colors[bin_num],
                                   linestyle='-', linewidth=2,
                                   label=f"{bin_label} - Canonical")
                # Cryptic - dashed line
                x_cry, y_cry = compute_cdf(bin_data[cryptic_col])
                if len(x_cry) > 0:
                    ax_shorter.plot(x_cry, y_cry, color=colors[bin_num],
                                   linestyle='--', linewidth=2,
                                   label=f"{bin_label} - Cryptic")

        ax_shorter.set_xlabel(metric_name, fontsize=11)
        ax_shorter.set_ylabel('CDF', fontsize=11)
        ax_shorter.set_title(f'{metric_name} - Shortening Events', fontsize=12, fontweight='bold')
        ax_shorter.grid(True, alpha=0.3)
        if row_idx == 0:  # Only show legend on first row
            ax_shorter.legend(fontsize=8, loc='best', ncol=2)

        # Right column: Longer events (cryptic_length > 0)
        ax_longer = axes[row_idx, 1]
        for bin_num in range(n_bins):
            bin_data = df_longer[df_longer['bin'] == bin_num]
            if len(bin_data) > 0:
                # Get bin info for label
                bin_info = [b for b in bin_info_longer if b['bin'] == bin_num][0]
                bin_label = f"Bin {bin_num+1} ({bin_info['min']} to {bin_info['max']} nt, n={bin_info['count']})"

                # Canonical - solid line
                x_can, y_can = compute_cdf(bin_data[canonical_col])
                if len(x_can) > 0:
                    ax_longer.plot(x_can, y_can, color=colors[bin_num],
                                  linestyle='-', linewidth=2,
                                  label=f"{bin_label} - Canonical")
                # Cryptic - dashed line
                x_cry, y_cry = compute_cdf(bin_data[cryptic_col])
                if len(x_cry) > 0:
                    ax_longer.plot(x_cry, y_cry, color=colors[bin_num],
                                  linestyle='--', linewidth=2,
                                  label=f"{bin_label} - Cryptic")

        ax_longer.set_xlabel(metric_name, fontsize=11)
        ax_longer.set_ylabel('CDF', fontsize=11)
        ax_longer.set_title(f'{metric_name} - Lengthening Events', fontsize=12, fontweight='bold')
        ax_longer.grid(True, alpha=0.3)
        if row_idx == 0:  # Only show legend on first row
            ax_longer.legend(fontsize=8, loc='best', ncol=2)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"    Saved CDF plot")


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description='Create binned CDF plots for splice site scores',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example (quantile-based binning):
  python create_score_cdf.py \\
    --reference conor_handedness.tsv \\
    --longer splice_site_analysis_longer.tsv \\
    --shorter splice_site_analysis_shorter.tsv \\
    --output cdf_analysis/ \\
    --bin 4 \\
    --quantile

Example (rank-based binning):
  python create_score_cdf.py \\
    --reference conor_handedness.tsv \\
    --longer splice_site_analysis_longer.tsv \\
    --shorter splice_site_analysis_shorter.tsv \\
    --output cdf_analysis/ \\
    --bin 4 \\
    --no-quantile
        """
    )

    # Same flags as create_score_heatmap.py except --cluster
    parser.add_argument('--reference', required=True,
                        help='Path to conor_handedness.tsv with reference scores')
    parser.add_argument('--longer', required=True,
                        help='Path to splice_site_analysis_longer.tsv')
    parser.add_argument('--shorter', required=True,
                        help='Path to splice_site_analysis_shorter.tsv')
    parser.add_argument('--output', required=True,
                        help='Output directory')
    parser.add_argument('--scores', nargs='+',
                        default=['Balance', 'SE.5SS.bed.score', 'LUC7_score'],
                        help='Score columns to visualize (default: Balance, SE.5SS.bed.score, LUC7_score)')
    parser.add_argument('--bin', type=int, default=4,
                        help='Number of bins for cryptic_length (default: 4)')
    parser.add_argument('--quantile', dest='quantile', action='store_true',
                        help='Use quantile-based binning (equal counts per bin) - default')
    parser.add_argument('--no-quantile', dest='quantile', action='store_false',
                        help='Use rank-based binning (equal rank ranges)')
    parser.set_defaults(quantile=True)

    args = parser.parse_args()

    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("BINNED CDF PLOT GENERATION")
    print("=" * 60)

    # Load reference scores
    print(f"\nLoading reference scores: {args.reference}")
    reference_df = pd.read_csv(args.reference, sep='\t')
    print(f"  Reference sequences: {len(reference_df)}")
    print(f"  Columns: {list(reference_df.columns)}")

    # Process longer events
    print("\n" + "-" * 60)
    print("PROCESSING LONGER EVENTS (cryptic_length > 0)")
    print("-" * 60)
    longer_df = load_and_process_data(args.longer, reference_df)
    longer_df, bin_info_longer = bin_data_by_cryptic_length(longer_df, n_bins=args.bin, use_quantile=args.quantile)
    print(f"  Binned into {len(bin_info_longer)} bins:")
    for info in bin_info_longer:
        print(f"    {info['label']} - {info['count']} events")

    # Process shorter events
    print("\n" + "-" * 60)
    print("PROCESSING SHORTER EVENTS (cryptic_length < 0)")
    print("-" * 60)
    shorter_df = load_and_process_data(args.shorter, reference_df)
    shorter_df, bin_info_shorter = bin_data_by_cryptic_length(shorter_df, n_bins=args.bin, use_quantile=args.quantile)
    print(f"  Binned into {len(bin_info_shorter)} bins:")
    for info in bin_info_shorter:
        print(f"    {info['label']} - {info['count']} events")

    # Create CDF plots
    print("\n" + "-" * 60)
    print("CREATING CDF PLOTS")
    print("-" * 60)
    cdf_plot_path = output_dir / 'cdf_plots.png'
    create_cdf_plots(longer_df, shorter_df, bin_info_longer, bin_info_shorter,
                     cdf_plot_path, n_bins=args.bin)

    # Save binned data
    binned_longer_path = output_dir / 'binned_longer.tsv'
    binned_shorter_path = output_dir / 'binned_shorter.tsv'
    longer_df.to_csv(binned_longer_path, sep='\t', index=False)
    shorter_df.to_csv(binned_shorter_path, sep='\t', index=False)
    print(f"  Saved binned data: {binned_longer_path}")
    print(f"  Saved binned data: {binned_shorter_path}")

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Output directory: {output_dir}")
    print(f"Number of bins: {args.bin}")
    print(f"Binning method: {'Quantile-based (equal counts)' if args.quantile else 'Rank-based (equal rank ranges)'}")
    print(f"Longer events: {len(longer_df)}")
    print(f"Shorter events: {len(shorter_df)}")
    print(f"\nGenerated files:")
    print(f"  - {cdf_plot_path}")
    print(f"  - {binned_longer_path}")
    print(f"  - {binned_shorter_path}")
    print("=" * 60)


if __name__ == '__main__':
    main()
