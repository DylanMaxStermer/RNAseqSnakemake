#!/usr/bin/env python3
"""
Splice Site Grid Visualization

Generates a grid visualization showing ranked splice site events with their
displacement from canonical splice sites, and extracts flanking sequences.

Layout per row:
[GREY padding (seq_range bp)] [BLACK region (cryptic_length)] [GREY padding (seq_range bp)]

Usage:
    python splice_site_grid_analysis.py \
        --productive path/to/productive_da_flanking \
        --reference path/to/Reference.fa \
        --output output_directory \
        --sequence-range 15
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import argparse
import os
import sys

try:
    import pyfaidx
except ImportError:
    print("ERROR: pyfaidx is required. Please activate the maxEnt conda environment.")
    sys.exit(1)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Simplified splice site grid visualization.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
    python splice_site_grid_analysis.py \\
        --productive productive_da_flanking \\
        --reference Reference.fa \\
        --output results/ \\
        --sequence-range 15
        """
    )

    parser.add_argument('--productive', required=True,
                        help='Path to productive_da_flanking file')
    parser.add_argument('--reference', required=True,
                        help='Path to reference genome FASTA file')
    parser.add_argument('--output', '-o', required=True,
                        help='Output directory')
    parser.add_argument('--chromosome', default='All',
                        help='Filter to specific chromosome (default: All)')
    parser.add_argument('--sequence-range', type=int, default=15,
                        help='bp upstream/downstream of splice site for padding (default: 15)')
    parser.add_argument('--IntronExon-range', type=int, default=None,
                        help='Max absolute cryptic_length to include (default: no filter)')
    parser.add_argument('--junction-type', default='altDonor',
                        choices=['altDonor', 'altAcceptor'],
                        help='Junction type to analyze (default: altDonor)')

    return parser.parse_args()


def parse_coordinate(coord_str):
    """
    Parse coordinate string 'chr:start-end' into (chr, start, end).

    Returns None if coordinate is invalid or NA.
    """
    if not coord_str or coord_str == 'NA' or pd.isna(coord_str):
        return None
    try:
        # Handle comma-separated (take first)
        if ',' in str(coord_str):
            coord_str = str(coord_str).split(',')[0]

        chrom, pos = str(coord_str).split(':')
        start, end = pos.split('-')
        return chrom, int(start), int(end)
    except (ValueError, AttributeError):
        return None


def calculate_intron_length(coord_str):
    """Calculate intron length from coordinate string."""
    parsed = parse_coordinate(coord_str)
    if parsed is None:
        return None
    return parsed[2] - parsed[1]  # end - start


def calculate_cryptic_length(row):
    """
    Calculate cryptic_length_absolute and signed cryptic_length.

    Sign determination:
    - if len(intron_coord) < len(spanning_productive) -> exon Lengthened (+)
    - if len(intron_coord) > len(spanning_productive) -> exon Shortened (-)

    Returns: (cryptic_length_absolute, cryptic_length)
    """
    # Get splice region size (cryptic_length_absolute)
    splice_region = row['productive_event_splice_region']
    parsed = parse_coordinate(splice_region)
    if parsed is None:
        return None, None

    cryptic_length_absolute = parsed[2] - parsed[1]  # end - start

    if cryptic_length_absolute == 0:
        return 0, 0

    # Compare intron lengths to determine sign
    intron_len = calculate_intron_length(row['intron_coord'])
    spanning_len = calculate_intron_length(row['spanning_productive'])

    if intron_len is None or spanning_len is None:
        return cryptic_length_absolute, None

    # Determine sign
    if intron_len < spanning_len:
        # Cryptic intron shorter = exon LENGTHENED = POSITIVE
        sign = 1
    elif intron_len > spanning_len:
        # Cryptic intron longer = exon SHORTENED = NEGATIVE
        sign = -1
    else:
        sign = 0

    cryptic_length = sign * cryptic_length_absolute
    return cryptic_length_absolute, cryptic_length


def get_donor_position(row, coord_col='intron_coord'):
    """
    Get the donor splice site position based on strand.

    For donor sites (5' splice site of intron):
    - Plus strand (+): Donor is at START of intron
    - Minus strand (-): Donor is at END of intron

    Returns: (chrom, position) or (None, None)
    """
    parsed = parse_coordinate(row[coord_col])
    if parsed is None:
        return None, None

    chrom, start, end = parsed
    strand = row['strand']

    if strand == '+':
        return chrom, start
    else:
        return chrom, end


def extract_splice_site_sequence(genome, chrom, position, seq_range, strand):
    """
    Extract sequence around a splice site position.

    Returns: (sequence_with_marker, coord_range) or (None, None)

    Format: upstream|downstream (e.g., "ACGTACGTACGTACG|TACGTACGTACGTAC")
    """
    if position is None or chrom is None:
        return None, None

    if chrom not in genome:
        return None, None

    chrom_len = len(genome[chrom])

    if strand == '+':
        up_start = position - seq_range
        up_end = position
        down_start = position
        down_end = position + seq_range

        if up_start < 0 or down_end > chrom_len:
            return None, None

        upstream = str(genome[chrom][up_start:up_end])
        downstream = str(genome[chrom][down_start:down_end])
        coord_range = f"{chrom}:{up_start}-{down_end}"

    else:  # Minus strand - reverse complement
        up_start = position
        up_end = position + seq_range
        down_start = position - seq_range
        down_end = position

        if down_start < 0 or up_end > chrom_len:
            return None, None

        upstream = str(genome[chrom][up_start:up_end].reverse.complement)
        downstream = str(genome[chrom][down_start:down_end].reverse.complement)
        coord_range = f"{chrom}:{down_start}-{up_end}"

    return upstream + '|' + downstream, coord_range


def extract_all_sequences(df, genome, seq_range):
    """Extract sequences around cryptic and canonical junctions."""
    print(f"  Extracting sequences ({seq_range} bp on each side)...")

    cryptic_sequences = []
    cryptic_ranges = []
    canonical_sequences = []
    canonical_ranges = []

    for idx, row in df.iterrows():
        if idx % 500 == 0 and idx > 0:
            print(f"    Processed {idx}/{len(df)} events...")

        strand = row['strand']

        # Get both donor positions
        chrom_cryptic, cryptic_pos = get_donor_position(row, 'intron_coord')
        chrom_canonical, canonical_pos = get_donor_position(row, 'spanning_productive')

        # Extract cryptic junction sequence
        if chrom_cryptic is not None:
            cryptic_seq, crypt_range = extract_splice_site_sequence(
                genome, chrom_cryptic, cryptic_pos, seq_range, strand
            )
        else:
            cryptic_seq, crypt_range = None, None

        # Extract canonical junction sequence
        if chrom_canonical is not None:
            canonical_seq, canon_range = extract_splice_site_sequence(
                genome, chrom_canonical, canonical_pos, seq_range, strand
            )
        else:
            canonical_seq, canon_range = None, None

        cryptic_sequences.append(cryptic_seq)
        cryptic_ranges.append(crypt_range)
        canonical_sequences.append(canonical_seq)
        canonical_ranges.append(canon_range)

    df['cryptic_sequence'] = cryptic_sequences
    df['cryptic_sequence_range'] = cryptic_ranges
    df['canonical_sequence'] = canonical_sequences
    df['canonical_sequence_range'] = canonical_ranges

    return df


def load_and_process_data(file_path, junction_type, chromosome, intron_exon_range):
    """Load data and calculate cryptic lengths."""
    print(f"  Loading data from: {file_path}")
    df = pd.read_csv(file_path, sep='\t')
    print(f"  Total rows: {len(df)}")

    # Filter by junction type
    df = df[df['junction_type'] == junction_type].copy()
    print(f"  After junction type filter ({junction_type}): {len(df)}")

    # Filter by chromosome if specified
    if chromosome != 'All':
        df = df[df['chrom'] == chromosome].copy()
        print(f"  After chromosome filter ({chromosome}): {len(df)}")

    # Calculate cryptic lengths
    print("  Calculating cryptic lengths...")
    results = df.apply(calculate_cryptic_length, axis=1)
    df['cryptic_length_absolute'] = [r[0] for r in results]
    df['cryptic_length'] = [r[1] for r in results]

    # Remove rows with None/NA cryptic_length
    df = df.dropna(subset=['cryptic_length'])
    df = df[df['cryptic_length'] != 0]  # Exclude zero-length events
    print(f"  After removing invalid/zero lengths: {len(df)}")

    # Filter for positive cryptic_length only (exon elongation events)
    before = len(df)
    df = df[df['cryptic_length'] > 0].copy()
    print(f"  After filtering for positive ΔSJ only: {len(df)} (removed {before - len(df)} negative events)")

    # Apply IntronExon-range filter
    if intron_exon_range is not None:
        before = len(df)
        df = df[df['cryptic_length_absolute'] <= intron_exon_range]
        print(f"  After IntronExon-range filter (max {intron_exon_range}): {len(df)} (removed {before - len(df)})")

    return df


def rank_events(df):
    """
    Rank events for visualization (positive ΔSJ only):
    - Order: largest at top, smallest at bottom
    - Sorted by cryptic_length descending
    """
    # Sort by cryptic_length descending (largest first)
    df_ranked = df.sort_values('cryptic_length', ascending=False).reset_index(drop=True)

    # Assign rank (1 = top row, largest ΔSJ)
    df_ranked['rank'] = range(1, len(df_ranked) + 1)

    return df_ranked


# Constants for coloring
BLACK = (0.0, 0.0, 0.0)
WHITE = (1.0, 1.0, 1.0)
GREY = (0.6, 0.6, 0.6)


def nucleotide_to_rgb(nucleotide):
    """Convert nucleotide to RGB color tuple (0-1 scale)."""
    color_map = {
        'T': (1.0, 0.0, 0.0),      # Red
        'G': (1.0, 0.65, 0.0),     # Orange
        'C': (0.0, 0.0, 1.0),      # Blue
        'A': (0.0, 0.8, 0.0),      # Green
        'N': (0.5, 0.5, 0.5),      # Grey for unknown
    }
    return color_map.get(nucleotide.upper(), (0.5, 0.5, 0.5))


def create_grid_visualization(df, output_path, seq_range):
    """
    Create simplified grid visualization.

    Layout per row:
    [GREY padding (seq_range bp)] [BLACK region (cryptic_length)] [GREY padding (seq_range bp)]

    - Y-axis: Ranked events (largest ΔSJ at top)
    - X-axis: Position in bp
    """
    print("  Creating simplified grid visualization...")

    n_events = len(df)
    if n_events == 0:
        print("  WARNING: No events to plot")
        return

    max_length = int(df['cryptic_length'].max())

    # Grid dimensions: padding + max_event_length + padding
    grid_width = seq_range + max_length + seq_range

    # Create RGB grid (white background)
    grid = np.ones((n_events, grid_width, 3))  # RGB array

    for row_idx, (_, event) in enumerate(df.iterrows()):
        cryptic_len = int(event['cryptic_length'])

        # Get sequences for coloring
        # canonical_sequence format: "upstream|downstream" - use upstream for left padding
        # cryptic_sequence format: "upstream|downstream" - use downstream for right padding
        canonical_seq = event.get('canonical_sequence', '')
        cryptic_seq = event.get('cryptic_sequence', '')

        if pd.notna(canonical_seq) and '|' in str(canonical_seq):
            left_seq = str(canonical_seq).split('|')[0]  # upstream of canonical
        else:
            left_seq = ''

        if pd.notna(cryptic_seq) and '|' in str(cryptic_seq):
            right_seq = str(cryptic_seq).split('|')[1]  # downstream of cryptic
        else:
            right_seq = ''

        col = 0

        # 1. Left padding - colored by canonical upstream sequence
        for i in range(seq_range):
            if col < grid_width:
                if i < len(left_seq):
                    grid[row_idx, col] = nucleotide_to_rgb(left_seq[i])
                else:
                    grid[row_idx, col] = GREY
            col += 1

        # 2. Black event region (cryptic_length bp)
        for _ in range(cryptic_len):
            if col < grid_width:
                grid[row_idx, col] = BLACK
            col += 1

        # 3. Right padding - colored by cryptic downstream sequence
        for i in range(seq_range):
            if col < grid_width:
                if i < len(right_seq):
                    grid[row_idx, col] = nucleotide_to_rgb(right_seq[i])
                else:
                    grid[row_idx, col] = GREY
            col += 1

        # Rest stays white (for events shorter than max)

    # Create figure - sized for visibility
    fig_height = max(6, min(20, n_events * 0.04))
    fig_width = max(12, min(24, grid_width * 0.08))
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Plot the grid
    extent = [0, grid_width, n_events, 0]
    ax.imshow(grid, aspect='auto', extent=extent, interpolation='nearest')

    # Add vertical lines to mark boundaries
    ax.axvline(x=seq_range, color='black', linestyle=':', linewidth=1, alpha=0.5)
    ax.axvline(x=seq_range + max_length, color='black', linestyle=':', linewidth=1, alpha=0.5)

    # Labels and formatting
    ax.set_xlabel('Position (bp)', fontsize=12)
    ax.set_ylabel('Ranked Events (largest ΔSJ at top)', fontsize=10)
    ax.set_title(f'Splice Site Spatial Distribution (Positive ΔSJ Only, n={n_events})', fontsize=14)

    # X-axis ticks
    tick_positions = [0, seq_range, seq_range + max_length // 2, seq_range + max_length, grid_width]
    tick_labels = [f'-{seq_range}', '0', f'+{max_length // 2}', f'+{max_length}', f'+{max_length + seq_range}']
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels)

    # Add legend
    legend_elements = [
        mpatches.Patch(facecolor=(0.0, 0.8, 0.0), label='A'),
        mpatches.Patch(facecolor=(0.0, 0.0, 1.0), label='C'),
        mpatches.Patch(facecolor=(1.0, 0.65, 0.0), label='G'),
        mpatches.Patch(facecolor=(1.0, 0.0, 0.0), label='T'),
        mpatches.Patch(facecolor=BLACK, label='Event region (ΔSJ)'),
        mpatches.Patch(facecolor=WHITE, edgecolor='black', label='Empty'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=8, ncol=2)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  Saved: {output_path}")


def print_summary_statistics(df):
    """Print summary statistics."""
    print("\n" + "=" * 60)
    print("SUMMARY STATISTICS (Positive ΔSJ / Exon Elongation Only)")
    print("=" * 60)

    print(f"\nTotal events: {len(df)}")

    if len(df) > 0:
        print(f"\nΔSJ (cryptic_length) statistics:")
        print(f"  Range: +{df['cryptic_length'].min():.0f} to +{df['cryptic_length'].max():.0f} bp")
        print(f"  Mean: +{df['cryptic_length'].mean():.1f} bp")
        print(f"  Median: +{df['cryptic_length'].median():.1f} bp")

        # Sequence extraction stats
        if 'cryptic_sequence' in df.columns:
            valid_cryptic = df['cryptic_sequence'].notna().sum()
            valid_canonical = df['canonical_sequence'].notna().sum()
            print(f"\nSequences extracted:")
            print(f"  Cryptic: {valid_cryptic}/{len(df)} ({100*valid_cryptic/len(df):.1f}%)")
            print(f"  Canonical: {valid_canonical}/{len(df)} ({100*valid_canonical/len(df):.1f}%)")


def main():
    args = parse_args()

    # Create output directory
    os.makedirs(args.output, exist_ok=True)

    print("=" * 70)
    print("SPLICE SITE GRID ANALYSIS")
    print("=" * 70)
    print(f"Input: {args.productive}")
    print(f"Reference: {args.reference}")
    print(f"Output: {args.output}")
    print(f"Sequence range: {args.sequence_range} bp")
    print(f"Junction type: {args.junction_type}")
    if args.IntronExon_range:
        print(f"IntronExon range filter: max {args.IntronExon_range} bp")
    print("=" * 70 + "\n")

    # Phase 1: Load and filter data
    print("Phase 1: Loading and processing data...")
    df = load_and_process_data(
        args.productive,
        junction_type=args.junction_type,
        chromosome=args.chromosome,
        intron_exon_range=args.IntronExon_range
    )

    if len(df) == 0:
        print("ERROR: No data remaining after filtering")
        sys.exit(1)

    # Phase 2: Rank events
    print("\nPhase 2: Ranking events by ΔSJ...")
    df = rank_events(df)
    print(f"  Events ranked: {len(df)}")

    # Phase 3: Extract sequences
    print("\nPhase 3: Extracting splice site sequences...")
    genome = pyfaidx.Fasta(args.reference)
    df = extract_all_sequences(df, genome, args.sequence_range)

    # Phase 4: Generate outputs
    print("\nPhase 4: Generating outputs...")

    # Tabular output
    output_tsv = os.path.join(args.output, 'splice_site_analysis.tsv')
    df.to_csv(output_tsv, sep='\t', index=False)
    print(f"  Saved tabular output: {output_tsv}")

    # Visualization
    output_png = os.path.join(args.output, 'splice_site_grid.png')
    create_grid_visualization(df, output_png, args.sequence_range)

    # Summary statistics
    print_summary_statistics(df)

    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)


if __name__ == '__main__':
    main()
