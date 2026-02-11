#!/usr/bin/env python3
"""
Create heatmaps of splice site scores for canonical and cryptic splice sites.

Integrates splice site analysis with reference scores to visualize how
Balance, Splice Score, and LUC7 scores correlate with cryptic site usage.
"""

import argparse
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
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
        9-mer string
    """
    if pd.isna(sequence) or '|' not in sequence:
        return None

    exon_part, intron_part = sequence.split('|')

    # Extract last 3 bases from exon and first 6 bases from intron
    if len(exon_part) >= 3 and len(intron_part) >= 6:
        return exon_part[-3:] + intron_part[:6]
    else:
        return None


def load_and_process_data(analysis_file, reference_df, cluster_by='none'):
    """
    Load analysis file and merge with reference scores.

    Args:
        analysis_file: Path to splice_site_analysis TSV
        reference_df: DataFrame with reference scores
        cluster_by: Clustering method ('none', 'balance', 'score', 'luc7')

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
    df = df.drop(columns=['SE.5SS.bed.seq_ref'], errors='ignore')

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
    df = df.drop(columns=['SE.5SS.bed.seq_ref2'], errors='ignore')

    # Count matches
    n_canonical_matched = df['canonical_score'].notna().sum()
    n_cryptic_matched = df['cryptic_score'].notna().sum()
    print(f"    Canonical matched: {n_canonical_matched}/{len(df)}")
    print(f"    Cryptic matched: {n_cryptic_matched}/{len(df)}")

    # Apply clustering/sorting
    if cluster_by == 'none':
        # Sort by rank (already ranked by |cryptic_length|)
        df = df.sort_values('rank')
        print("  Sorting by original rank (|ΔSJ|)")
    elif cluster_by == 'balance':
        df = df.sort_values('canonical_Balance', ascending=False)
        print("  Clustering by canonical Balance score")
    elif cluster_by == 'score':
        df = df.sort_values('canonical_score', ascending=False)
        print("  Clustering by canonical Splice Score")
    elif cluster_by == 'luc7':
        df = df.sort_values('canonical_LUC7', ascending=False)
        print("  Clustering by canonical LUC7 score")

    return df


def create_heatmap(df, output_path, event_type, cluster_by='none'):
    """
    Create heatmap visualization with 8 rows × 1 column (all vertical).

    Layout: For each metric, canonical appears first, then cryptic below it.
    - Row 0: Balance - Canonical
    - Row 1: Balance - Cryptic
    - Row 2: Splice Score - Canonical
    - Row 3: Splice Score - Cryptic
    - Row 4: LUC7 - Canonical
    - Row 5: LUC7 - Cryptic
    - Row 6: Rank - Canonical
    - Row 7: Rank - Cryptic

    Color scales:
    - Balance: RdBu_r, range -3 to 3
    - Splice Score: viridis, data-driven range
    - LUC7: plasma, data-driven range
    - Rank: Blues (light to dark blue), continuous

    Args:
        df: DataFrame with matched scores
        output_path: Path to save heatmap image
        event_type: String describing event type (e.g., "Exon Longer")
        cluster_by: Clustering method ('none', 'balance', 'score', 'luc7')
    """
    print(f"  Creating heatmap: {output_path}")

    n_events = len(df)

    # Create figure with 8 rows × 1 column (all vertical, thinner rows)
    fig, axes = plt.subplots(8, 1, figsize=(20, 10))

    # Prepare individual metric arrays
    canonical_balance = df['canonical_Balance'].fillna(0).values.reshape(1, -1)
    canonical_score = df['canonical_score'].fillna(df['canonical_score'].mean()).values.reshape(1, -1)
    canonical_luc7 = df['canonical_LUC7'].fillna(0).values.reshape(1, -1)

    cryptic_balance = df['cryptic_Balance'].fillna(0).values.reshape(1, -1)
    cryptic_score = df['cryptic_score'].fillna(df['cryptic_score'].mean()).values.reshape(1, -1)
    cryptic_luc7 = df['cryptic_LUC7'].fillna(0).values.reshape(1, -1)

    # Rank array (same for both canonical and cryptic)
    rank_array = df['rank'].values.reshape(1, -1)
    max_rank = df['rank'].max()

    # Calculate ranges
    score_min = min(canonical_score.min(), cryptic_score.min())
    score_max = max(canonical_score.max(), cryptic_score.max())
    luc7_min = min(canonical_luc7.min(), cryptic_luc7.min())
    luc7_max = max(canonical_luc7.max(), cryptic_luc7.max())

    # Row 0: Balance - Canonical
    sns.heatmap(canonical_balance, ax=axes[0],
                cmap='RdBu_r', center=0, vmin=-3, vmax=3,
                yticklabels=['Canonical'], xticklabels=False,
                cbar_kws={'label': 'Balance', 'ticks': [-3, -2, -1, 0, 1, 2, 3]})
    axes[0].set_ylabel('Balance\nCanonical', fontsize=10)

    # Row 1: Balance - Cryptic
    sns.heatmap(cryptic_balance, ax=axes[1],
                cmap='RdBu_r', center=0, vmin=-3, vmax=3,
                yticklabels=['Cryptic'], xticklabels=False,
                cbar_kws={'label': 'Balance', 'ticks': [-3, -2, -1, 0, 1, 2, 3]})
    axes[1].set_ylabel('Balance\nCryptic', fontsize=10)

    # Row 2: Splice Score - Canonical
    sns.heatmap(canonical_score, ax=axes[2],
                cmap='viridis', vmin=score_min, vmax=score_max,
                yticklabels=['Canonical'], xticklabels=False,
                cbar_kws={'label': 'Splice Score'})
    axes[2].set_ylabel('Splice Score\nCanonical', fontsize=10)

    # Row 3: Splice Score - Cryptic
    sns.heatmap(cryptic_score, ax=axes[3],
                cmap='viridis', vmin=score_min, vmax=score_max,
                yticklabels=['Cryptic'], xticklabels=False,
                cbar_kws={'label': 'Splice Score'})
    axes[3].set_ylabel('Splice Score\nCryptic', fontsize=10)

    # Row 4: LUC7 - Canonical
    sns.heatmap(canonical_luc7, ax=axes[4],
                cmap='plasma', vmin=luc7_min, vmax=luc7_max,
                yticklabels=['Canonical'], xticklabels=False,
                cbar_kws={'label': 'LUC7 Score'})
    axes[4].set_ylabel('LUC7\nCanonical', fontsize=10)

    # Row 5: LUC7 - Cryptic
    sns.heatmap(cryptic_luc7, ax=axes[5],
                cmap='plasma', vmin=luc7_min, vmax=luc7_max,
                yticklabels=['Cryptic'], xticklabels=False,
                cbar_kws={'label': 'LUC7 Score'})
    axes[5].set_ylabel('LUC7\nCryptic', fontsize=10)

    # Row 6: Rank - Canonical (light blue to dark blue)
    sns.heatmap(rank_array, ax=axes[6],
                cmap='Blues', vmin=1, vmax=max_rank,
                yticklabels=['Canonical'], xticklabels=False,
                cbar_kws={'label': 'Rank'})
    axes[6].set_ylabel('Rank\nCanonical', fontsize=10)

    # Row 7: Rank - Cryptic (light blue to dark blue)
    sns.heatmap(rank_array, ax=axes[7],
                cmap='Blues', vmin=1, vmax=max_rank,
                yticklabels=['Cryptic'], xticklabels=False,
                cbar_kws={'label': 'Rank'})
    axes[7].set_ylabel('Rank\nCryptic', fontsize=10)
    axes[7].set_xlabel('Events (ranked by |ΔSJ|, largest to smallest)', fontsize=10)

    # Overall title with clustering info
    if cluster_by == 'none':
        title = f'Splice Site Scores: {event_type}'
    elif cluster_by == 'balance':
        title = f'Splice Site Scores: {event_type}\n(Clustered by Canonical Balance)'
    elif cluster_by == 'score':
        title = f'Splice Site Scores: {event_type}\n(Clustered by Canonical Splice Score)'
    elif cluster_by == 'luc7':
        title = f'Splice Site Scores: {event_type}\n(Clustered by Canonical LUC7)'

    fig.suptitle(title, fontsize=14, fontweight='bold')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"    Saved heatmap with {n_events} events")


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description='Create splice site score heatmaps',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python create_score_heatmap.py \\
    --reference conor_handedness.tsv \\
    --longer splice_site_analysis_longer.tsv \\
    --shorter splice_site_analysis_shorter.tsv \\
    --output heatmap_analysis/

Example with clustering by canonical LUC7 score:
  python create_score_heatmap.py \\
    --reference conor_handedness.tsv \\
    --longer splice_site_analysis_longer.tsv \\
    --shorter splice_site_analysis_shorter.tsv \\
    --output heatmap_analysis/ \\
    --cluster luc7
        """
    )

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
    parser.add_argument('--cluster', choices=['none', 'balance', 'score', 'luc7'],
                        default='none',
                        help='Cluster events by canonical site metric (default: none, rank order)')

    args = parser.parse_args()

    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("SPLICE SITE SCORE HEATMAP GENERATION")
    print("=" * 60)

    # Load reference scores
    print(f"\nLoading reference scores: {args.reference}")
    reference_df = pd.read_csv(args.reference, sep='\t')
    print(f"  Reference sequences: {len(reference_df)}")
    print(f"  Columns: {list(reference_df.columns)}")

    # Process longer events
    print("\n" + "-" * 60)
    print("PROCESSING LONGER EVENTS (ΔSJ > 0)")
    print("-" * 60)
    longer_df = load_and_process_data(args.longer, reference_df, cluster_by=args.cluster)

    heatmap_longer_path = output_dir / 'heatmap_longer.png'
    create_heatmap(longer_df, heatmap_longer_path, 'Exon Longer (ΔSJ > 0)', cluster_by=args.cluster)

    scored_longer_path = output_dir / 'scored_longer.tsv'
    longer_df.to_csv(scored_longer_path, sep='\t', index=False)
    print(f"  Saved scored data: {scored_longer_path}")

    # Process shorter events
    print("\n" + "-" * 60)
    print("PROCESSING SHORTER EVENTS (ΔSJ < 0)")
    print("-" * 60)
    shorter_df = load_and_process_data(args.shorter, reference_df, cluster_by=args.cluster)

    heatmap_shorter_path = output_dir / 'heatmap_shorter.png'
    create_heatmap(shorter_df, heatmap_shorter_path, 'Exon Shorter (ΔSJ < 0)', cluster_by=args.cluster)

    scored_shorter_path = output_dir / 'scored_shorter.tsv'
    shorter_df.to_csv(scored_shorter_path, sep='\t', index=False)
    print(f"  Saved scored data: {scored_shorter_path}")

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Output directory: {output_dir}")
    print(f"Longer events: {len(longer_df)}")
    print(f"Shorter events: {len(shorter_df)}")
    print(f"\nGenerated files:")
    print(f"  - {heatmap_longer_path}")
    print(f"  - {heatmap_shorter_path}")
    print(f"  - {scored_longer_path}")
    print(f"  - {scored_shorter_path}")
    print("=" * 60)


if __name__ == '__main__':
    main()
