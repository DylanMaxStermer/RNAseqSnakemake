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
    parser.add_argument('--exon-change', default='longer',
                        choices=['longer', 'shorter'],
                        help='Filter by exon change type: longer (cryptic_length > 0) or shorter (cryptic_length < 0)')
    parser.add_argument('--merge', action='store_true', default=False,
                        help='Plot both longer and shorter on same plot (longer=right/+X, shorter=left/-X)')

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


def load_and_process_data(file_path, junction_type, chromosome, intron_exon_range, exon_change):
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

    # Filter based on exon change type
    before = len(df)
    if exon_change == 'longer':
        df = df[df['cryptic_length'] > 0].copy()
        print(f"  After filtering for positive ΔSJ (exon longer): {len(df)} (removed {before - len(df)} events)")
    else:  # shorter
        df = df[df['cryptic_length'] < 0].copy()
        print(f"  After filtering for negative ΔSJ (exon shorter): {len(df)} (removed {before - len(df)} events)")

    # Apply IntronExon-range filter
    if intron_exon_range is not None:
        before = len(df)
        df = df[df['cryptic_length_absolute'] <= intron_exon_range]
        print(f"  After IntronExon-range filter (max {intron_exon_range}): {len(df)} (removed {before - len(df)})")

    return df


def rank_events(df, exon_change):
    """
    Rank events for visualization:
    - Order: largest absolute value at top, smallest at bottom
    - For longer: sorted by cryptic_length descending
    - For shorter: sorted by cryptic_length ascending (most negative first)
    """
    if exon_change == 'longer':
        # Sort by cryptic_length descending (largest positive first)
        df_ranked = df.sort_values('cryptic_length', ascending=False).reset_index(drop=True)
    else:
        # Sort by cryptic_length ascending (most negative first = largest absolute value)
        df_ranked = df.sort_values('cryptic_length', ascending=True).reset_index(drop=True)

    # Assign rank (1 = top row, largest |ΔSJ|)
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


def create_grid_array(df, seq_range, exon_change):
    """
    Create grid array for visualization.

    Returns:
        grid: numpy array (n_events, grid_width, 3) - RGB grid
        max_length: int - maximum cryptic_length in dataset
        grid_width: int - total grid width
        junction_info: list of tuples (cryptic_len) for each row
    """
    n_events = len(df)
    if n_events == 0:
        return None, 0, 0, []

    # Use absolute value for sizing
    max_length = int(df['cryptic_length_absolute'].max())

    # Grid dimensions: padding + max_event_length + padding
    grid_width = seq_range + max_length + seq_range

    # Create RGB grid (white background)
    grid = np.ones((n_events, grid_width, 3))  # RGB array

    junction_info = []

    for row_idx, (_, event) in enumerate(df.iterrows()):
        # Use absolute value for event length
        cryptic_len = int(abs(event['cryptic_length']))
        junction_info.append(cryptic_len)

        # Get sequences for coloring
        canonical_seq = event.get('canonical_sequence', '')
        cryptic_seq = event.get('cryptic_sequence', '')

        if pd.notna(canonical_seq) and '|' in str(canonical_seq):
            canonical_upstream = str(canonical_seq).split('|')[0]
            canonical_downstream = str(canonical_seq).split('|')[1]
        else:
            canonical_upstream = ''
            canonical_downstream = ''

        if pd.notna(cryptic_seq) and '|' in str(cryptic_seq):
            cryptic_upstream = str(cryptic_seq).split('|')[0]
            cryptic_downstream = str(cryptic_seq).split('|')[1]
        else:
            cryptic_upstream = ''
            cryptic_downstream = ''

        if exon_change == 'longer':
            left_upstream = canonical_upstream
            left_downstream = canonical_downstream
            right_upstream = cryptic_upstream
            right_downstream = cryptic_downstream
        else:  # shorter
            left_upstream = cryptic_upstream
            left_downstream = cryptic_downstream
            right_upstream = canonical_upstream
            right_downstream = canonical_downstream

        col = 0

        # 1. Left padding
        if exon_change == 'longer':
            left_padding_size = seq_range
        else:  # shorter
            left_padding_size = seq_range + max_length - cryptic_len

        for i in range(left_padding_size):
            if col < grid_width:
                seq_idx = len(left_upstream) - left_padding_size + i
                if seq_idx >= 0 and seq_idx < len(left_upstream):
                    grid[row_idx, col] = nucleotide_to_rgb(left_upstream[seq_idx])
                else:
                    grid[row_idx, col] = GREY
            col += 1

        # 2. Event region
        if cryptic_len >= 2 * seq_range:
            for i in range(seq_range):
                if col < grid_width:
                    if i < len(left_downstream):
                        grid[row_idx, col] = nucleotide_to_rgb(left_downstream[i])
                    else:
                        grid[row_idx, col] = GREY
                col += 1
            gap_size = cryptic_len - 2 * seq_range
            for _ in range(gap_size):
                if col < grid_width:
                    grid[row_idx, col] = BLACK
                col += 1
            for i in range(seq_range):
                if col < grid_width:
                    if i < len(right_upstream):
                        grid[row_idx, col] = nucleotide_to_rgb(right_upstream[i])
                    else:
                        grid[row_idx, col] = GREY
                col += 1
        elif cryptic_len > seq_range:
            for i in range(seq_range):
                if col < grid_width:
                    if i < len(left_downstream):
                        grid[row_idx, col] = nucleotide_to_rgb(left_downstream[i])
                    else:
                        grid[row_idx, col] = GREY
                col += 1
            remainder = cryptic_len - seq_range
            start_idx = max(0, seq_range - remainder)
            for i in range(remainder):
                if col < grid_width:
                    idx = start_idx + i
                    if idx < len(right_upstream):
                        grid[row_idx, col] = nucleotide_to_rgb(right_upstream[idx])
                    else:
                        grid[row_idx, col] = GREY
                col += 1
        else:
            for i in range(cryptic_len):
                if col < grid_width:
                    if i < len(left_downstream):
                        grid[row_idx, col] = nucleotide_to_rgb(left_downstream[i])
                    else:
                        grid[row_idx, col] = GREY
                col += 1

        # 3. Right padding
        for i in range(seq_range):
            if col < grid_width:
                if i < len(right_downstream):
                    grid[row_idx, col] = nucleotide_to_rgb(right_downstream[i])
                else:
                    grid[row_idx, col] = GREY
            col += 1

    return grid, max_length, grid_width, junction_info


def create_grid_visualization(df, output_path, seq_range, exon_change):
    """
    Create grid visualization with sequence coloring.

    For longer (cryptic_length > 0):
    [canonical_upstream] [■] [event region] [■] [cryptic_downstream]

    For shorter (cryptic_length < 0):
    [cryptic_upstream] [■] [event region] [■] [canonical_downstream]

    - Y-axis: Ranked events (largest |ΔSJ| at top)
    - X-axis: Position in bp (negative for shorter)
    """
    print("  Creating grid visualization...")

    # Create grid array using helper function
    grid, max_length, grid_width, junction_info = create_grid_array(df, seq_range, exon_change)

    if grid is None:
        print("  WARNING: No events to plot")
        return

    n_events = len(df)

    # Create figure - sized for visibility
    fig_height = max(6, min(20, n_events * 0.04))
    fig_width = max(12, min(24, grid_width * 0.08))
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Plot the grid
    extent = [0, grid_width, n_events, 0]
    ax.imshow(grid, aspect='auto', extent=extent, interpolation='nearest')

    # Add junction markers
    if exon_change == 'longer':
        # For longer: canonical junction at seq_range (fixed, vertical line at X=0)
        ax.axvline(x=seq_range, color='black', linestyle='-', linewidth=2, label='Canonical Junction (0)')

        # Alternative junctions vary per-row (at positive X positions)
        # Also connect tops of adjacent lines to create a "staircase" outline
        prev_alt_junc_x = None
        for row_idx, (_, event) in enumerate(df.iterrows()):
            cryptic_len = int(abs(event['cryptic_length']))
            alt_junc_x = seq_range + cryptic_len

            # Draw vertical line segment for this row
            ax.plot([alt_junc_x, alt_junc_x], [row_idx, row_idx + 1],
                    color='black', linewidth=2, solid_capstyle='butt')

            # Connect to previous row's junction with a horizontal line
            if prev_alt_junc_x is not None:
                ax.plot([prev_alt_junc_x, alt_junc_x], [row_idx, row_idx],
                        color='black', linewidth=2, solid_capstyle='butt')

            prev_alt_junc_x = alt_junc_x
    else:  # shorter
        # For shorter: canonical junction at seq_range + max_length (fixed, vertical line at X=0)
        canonical_junc_x = seq_range + max_length
        ax.axvline(x=canonical_junc_x, color='black', linestyle='-', linewidth=2, label='Canonical Junction (0)')

        # Alternative junctions vary per-row (at negative X positions)
        # Also connect tops of adjacent lines to create a "staircase" outline
        prev_alt_junc_x = None
        for row_idx, (_, event) in enumerate(df.iterrows()):
            cryptic_len = int(abs(event['cryptic_length']))
            alt_junc_x = seq_range + max_length - cryptic_len

            # Draw vertical line segment for this row
            ax.plot([alt_junc_x, alt_junc_x], [row_idx, row_idx + 1],
                    color='black', linewidth=2, solid_capstyle='butt')

            # Connect to previous row's junction with a horizontal line
            if prev_alt_junc_x is not None:
                ax.plot([prev_alt_junc_x, alt_junc_x], [row_idx, row_idx],
                        color='black', linewidth=2, solid_capstyle='butt')

            prev_alt_junc_x = alt_junc_x

    # Labels and formatting
    ax.set_xlabel('Position (bp)', fontsize=12)
    ax.set_ylabel('Ranked Events (largest |ΔSJ| at top)', fontsize=10)

    if exon_change == 'longer':
        ax.set_title(f'Splice Site Spatial Distribution (Exon Longer, ΔSJ > 0, n={n_events})', fontsize=14)
        # X-axis ticks for positive direction
        tick_positions = [0, seq_range, seq_range + max_length // 2, seq_range + max_length, grid_width]
        tick_labels = [f'-{seq_range}', '0', f'+{max_length // 2}', f'+{max_length}', f'+{max_length + seq_range}']
    else:
        ax.set_title(f'Splice Site Spatial Distribution (Exon Shorter, ΔSJ < 0, n={n_events})', fontsize=14)
        # Canonical junction at position seq_range + max_length should be labeled '0'
        # Positions to the left are negative, positions to the right are positive
        canonical_pos = seq_range + max_length
        tick_positions = [0, seq_range, canonical_pos - max_length // 2, canonical_pos, grid_width]
        tick_labels = [
            f'-{canonical_pos}',           # Left edge: -(seq_range + max_length)
            f'-{max_length}',               # Cryptic junction: -max_length from canonical
            f'-{max_length // 2}',          # Midpoint
            '0',                            # Canonical junction: reference
            f'+{seq_range}'                 # Right edge: +seq_range from canonical
        ]

    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels)

    # Add legend
    legend_elements = [
        mpatches.Patch(facecolor=(0.0, 0.8, 0.0), label='A'),
        mpatches.Patch(facecolor=(0.0, 0.0, 1.0), label='C'),
        mpatches.Patch(facecolor=(1.0, 0.65, 0.0), label='G'),
        mpatches.Patch(facecolor=(1.0, 0.0, 0.0), label='T'),
        mpatches.Patch(facecolor=BLACK, label='Gap'),
        mpatches.Patch(facecolor=WHITE, edgecolor='black', label='Empty'),
        plt.Line2D([0], [0], color='black', linewidth=2, label='Junction (|)'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=8, ncol=2)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  Saved: {output_path}")


def create_merged_visualization(df_longer, df_shorter, output_path, seq_range):
    """
    Create merged visualization with both longer and shorter events side-by-side.

    Uses two separate subplots:
    - Left: Shorter events (ΔSJ < 0) with canonical junctions aligned
    - Right: Longer events (ΔSJ > 0) with canonical junctions aligned
    """
    print("  Creating merged grid visualization (side-by-side subplots)...")

    n_longer = len(df_longer)
    n_shorter = len(df_shorter)

    if n_longer == 0 and n_shorter == 0:
        print("  WARNING: No events to plot")
        return

    # Create grid arrays using helper function
    grid_shorter, max_shorter, width_shorter, junc_shorter = create_grid_array(
        df_shorter, seq_range, 'shorter'
    ) if n_shorter > 0 else (None, 0, 0, [])

    grid_longer, max_longer, width_longer, junc_longer = create_grid_array(
        df_longer, seq_range, 'longer'
    ) if n_longer > 0 else (None, 0, 0, [])

    # Create figure with two subplots side-by-side
    n_events = max(n_longer, n_shorter)
    fig_height = max(6, min(20, n_events * 0.04))
    fig_width = max(16, min(32, (width_shorter + width_longer) * 0.06))

    fig, (ax_shorter, ax_longer) = plt.subplots(1, 2, figsize=(fig_width, fig_height))

    # Plot shorter events (left subplot)
    if grid_shorter is not None:
        extent_shorter = [0, width_shorter, n_shorter, 0]
        ax_shorter.imshow(grid_shorter, aspect='auto', extent=extent_shorter, interpolation='nearest')

        # Add junction markers for shorter
        canonical_junc_x = seq_range + max_shorter
        ax_shorter.axvline(x=canonical_junc_x, color='black', linestyle='-', linewidth=2, label='Canonical Junction (0)')

        # Alternative junctions vary per-row
        prev_alt_junc_x = None
        for row_idx, cryptic_len in enumerate(junc_shorter):
            alt_junc_x = seq_range + max_shorter - cryptic_len
            ax_shorter.plot([alt_junc_x, alt_junc_x], [row_idx, row_idx + 1],
                           color='black', linewidth=2, solid_capstyle='butt')
            if prev_alt_junc_x is not None:
                ax_shorter.plot([prev_alt_junc_x, alt_junc_x], [row_idx, row_idx],
                               color='black', linewidth=2, solid_capstyle='butt')
            prev_alt_junc_x = alt_junc_x

        # Labels and ticks for shorter
        ax_shorter.set_xlabel('Position (bp)', fontsize=12)
        ax_shorter.set_ylabel('Ranked Events (largest |ΔSJ| at top)', fontsize=10)
        ax_shorter.set_title(f'Exon Shorter (ΔSJ < 0, n={n_shorter})', fontsize=14)

        canonical_pos = seq_range + max_shorter
        tick_positions = [0, seq_range, canonical_pos - max_shorter // 2, canonical_pos, width_shorter]
        tick_labels = [
            f'-{canonical_pos}',
            f'-{max_shorter}',
            f'-{max_shorter // 2}',
            '0',
            f'+{seq_range}'
        ]
        ax_shorter.set_xticks(tick_positions)
        ax_shorter.set_xticklabels(tick_labels)

    else:
        ax_shorter.text(0.5, 0.5, 'No shorter events', ha='center', va='center', transform=ax_shorter.transAxes)
        ax_shorter.set_title('Exon Shorter (ΔSJ < 0, n=0)', fontsize=14)

    # Plot longer events (right subplot)
    if grid_longer is not None:
        extent_longer = [0, width_longer, n_longer, 0]
        ax_longer.imshow(grid_longer, aspect='auto', extent=extent_longer, interpolation='nearest')

        # Add junction markers for longer
        ax_longer.axvline(x=seq_range, color='black', linestyle='-', linewidth=2, label='Canonical Junction (0)')

        # Alternative junctions vary per-row
        prev_alt_junc_x = None
        for row_idx, cryptic_len in enumerate(junc_longer):
            alt_junc_x = seq_range + cryptic_len
            ax_longer.plot([alt_junc_x, alt_junc_x], [row_idx, row_idx + 1],
                          color='black', linewidth=2, solid_capstyle='butt')
            if prev_alt_junc_x is not None:
                ax_longer.plot([prev_alt_junc_x, alt_junc_x], [row_idx, row_idx],
                              color='black', linewidth=2, solid_capstyle='butt')
            prev_alt_junc_x = alt_junc_x

        # Labels and ticks for longer
        ax_longer.set_xlabel('Position (bp)', fontsize=12)
        ax_longer.set_ylabel('Ranked Events (largest |ΔSJ| at top)', fontsize=10)
        ax_longer.set_title(f'Exon Longer (ΔSJ > 0, n={n_longer})', fontsize=14)

        tick_positions = [0, seq_range, seq_range + max_longer // 2, seq_range + max_longer, width_longer]
        tick_labels = [
            f'-{seq_range}',
            '0',
            f'+{max_longer // 2}',
            f'+{max_longer}',
            f'+{max_longer + seq_range}'
        ]
        ax_longer.set_xticks(tick_positions)
        ax_longer.set_xticklabels(tick_labels)

    else:
        ax_longer.text(0.5, 0.5, 'No longer events', ha='center', va='center', transform=ax_longer.transAxes)
        ax_longer.set_title('Exon Longer (ΔSJ > 0, n=0)', fontsize=14)

    # Add legend to the longer subplot
    legend_elements = [
        mpatches.Patch(facecolor=(0.0, 0.8, 0.0), label='A'),
        mpatches.Patch(facecolor=(0.0, 0.0, 1.0), label='C'),
        mpatches.Patch(facecolor=(1.0, 0.65, 0.0), label='G'),
        mpatches.Patch(facecolor=(1.0, 0.0, 0.0), label='T'),
        mpatches.Patch(facecolor=BLACK, label='Gap'),
        mpatches.Patch(facecolor=WHITE, edgecolor='black', label='Empty'),
        plt.Line2D([0], [0], color='black', linewidth=2, label='Junction (|)'),
    ]
    ax_longer.legend(handles=legend_elements, loc='upper right', fontsize=8, ncol=2)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  Saved: {output_path}")


def print_summary_statistics(df, exon_change):
    """Print summary statistics."""
    print("\n" + "=" * 60)
    if exon_change == 'longer':
        print("SUMMARY STATISTICS (Exon Longer / ΔSJ > 0)")
    else:
        print("SUMMARY STATISTICS (Exon Shorter / ΔSJ < 0)")
    print("=" * 60)

    print(f"\nTotal events: {len(df)}")

    if len(df) > 0:
        print(f"\nΔSJ (cryptic_length) statistics:")
        min_val = df['cryptic_length'].min()
        max_val = df['cryptic_length'].max()
        mean_val = df['cryptic_length'].mean()
        median_val = df['cryptic_length'].median()

        if exon_change == 'longer':
            print(f"  Range: +{min_val:.0f} to +{max_val:.0f} bp")
            print(f"  Mean: +{mean_val:.1f} bp")
            print(f"  Median: +{median_val:.1f} bp")
        else:
            print(f"  Range: {max_val:.0f} to {min_val:.0f} bp")
            print(f"  Mean: {mean_val:.1f} bp")
            print(f"  Median: {median_val:.1f} bp")

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
    if args.merge:
        print(f"Mode: MERGED (both longer and shorter)")
    else:
        print(f"Exon change: {args.exon_change}")
    if args.IntronExon_range:
        print(f"IntronExon range filter: max {args.IntronExon_range} bp")
    print("=" * 70 + "\n")

    # Load genome for sequence extraction
    genome = pyfaidx.Fasta(args.reference)

    if args.merge:
        # Merged mode: load both longer and shorter
        print("Phase 1: Loading and processing data (MERGED MODE)...")

        # Load longer events
        print("  Loading longer events...")
        df_longer = load_and_process_data(
            args.productive,
            junction_type=args.junction_type,
            chromosome=args.chromosome,
            intron_exon_range=args.IntronExon_range,
            exon_change='longer'
        )

        # Load shorter events
        print("  Loading shorter events...")
        df_shorter = load_and_process_data(
            args.productive,
            junction_type=args.junction_type,
            chromosome=args.chromosome,
            intron_exon_range=args.IntronExon_range,
            exon_change='shorter'
        )

        if len(df_longer) == 0 and len(df_shorter) == 0:
            print("ERROR: No data remaining after filtering")
            sys.exit(1)

        # Phase 2: Rank events
        print("\nPhase 2: Ranking events by |ΔSJ|...")
        if len(df_longer) > 0:
            df_longer = rank_events(df_longer, 'longer')
            print(f"  Longer events ranked: {len(df_longer)}")
        if len(df_shorter) > 0:
            df_shorter = rank_events(df_shorter, 'shorter')
            print(f"  Shorter events ranked: {len(df_shorter)}")

        # Phase 3: Extract sequences
        print("\nPhase 3: Extracting splice site sequences...")
        if len(df_longer) > 0:
            df_longer = extract_all_sequences(df_longer, genome, args.sequence_range)
        if len(df_shorter) > 0:
            df_shorter = extract_all_sequences(df_shorter, genome, args.sequence_range)

        # Phase 4: Generate outputs
        print("\nPhase 4: Generating outputs...")

        # Tabular outputs
        if len(df_longer) > 0:
            output_tsv_longer = os.path.join(args.output, 'splice_site_analysis_longer.tsv')
            df_longer.to_csv(output_tsv_longer, sep='\t', index=False)
            print(f"  Saved: {output_tsv_longer}")
        if len(df_shorter) > 0:
            output_tsv_shorter = os.path.join(args.output, 'splice_site_analysis_shorter.tsv')
            df_shorter.to_csv(output_tsv_shorter, sep='\t', index=False)
            print(f"  Saved: {output_tsv_shorter}")

        # Merged visualization
        output_png = os.path.join(args.output, 'splice_site_grid_merged.png')
        create_merged_visualization(df_longer, df_shorter, output_png, args.sequence_range)

        # Summary statistics
        if len(df_longer) > 0:
            print_summary_statistics(df_longer, 'longer')
        if len(df_shorter) > 0:
            print_summary_statistics(df_shorter, 'shorter')

    else:
        # Single mode: load only specified exon_change type
        print("Phase 1: Loading and processing data...")
        df = load_and_process_data(
            args.productive,
            junction_type=args.junction_type,
            chromosome=args.chromosome,
            intron_exon_range=args.IntronExon_range,
            exon_change=args.exon_change
        )

        if len(df) == 0:
            print("ERROR: No data remaining after filtering")
            sys.exit(1)

        # Phase 2: Rank events
        print("\nPhase 2: Ranking events by |ΔSJ|...")
        df = rank_events(df, args.exon_change)
        print(f"  Events ranked: {len(df)}")

        # Phase 3: Extract sequences
        print("\nPhase 3: Extracting splice site sequences...")
        df = extract_all_sequences(df, genome, args.sequence_range)

        # Phase 4: Generate outputs
        print("\nPhase 4: Generating outputs...")

        # Tabular output
        output_tsv = os.path.join(args.output, 'splice_site_analysis.tsv')
        df.to_csv(output_tsv, sep='\t', index=False)
        print(f"  Saved tabular output: {output_tsv}")

        # Visualization
        output_png = os.path.join(args.output, 'splice_site_grid.png')
        create_grid_visualization(df, output_png, args.sequence_range, args.exon_change)

        # Summary statistics
        print_summary_statistics(df, args.exon_change)

    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)


if __name__ == '__main__':
    main()
