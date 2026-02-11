#!/usr/bin/env python3
"""
Splice Site Distribution Analysis Tool

Analyzes the distribution of distances between alternative and canonical splice sites,
categorizing them as exon elongations or shortenings, and generates comparative
visualizations for Poison vs Productive events for both Donor and Acceptor splice sites.

Usage:
    python splice_site_distribution.py \
        --poison path/to/poison_da_flanking.tsv \
        --productive path/to/productive_da_flanking \
        --output output_directory
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Analyze splice site distance distributions between alternative and canonical splice sites.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
    python splice_site_distribution.py \\
        --poison poison_da_flanking.tsv \\
        --productive productive_da_flanking \\
        --output results/
        """
    )

    parser.add_argument('--poison', required=True,
                        help='Path to poison_da_flanking.tsv')
    parser.add_argument('--productive', required=True,
                        help='Path to productive_da_flanking file')
    parser.add_argument('--output', '-o', required=True,
                        help='Output directory for plots')
    parser.add_argument('--chromosome', default='All',
                        help='Filter to specific chromosome (default: All)')
    parser.add_argument('--intron-range', type=int, default=200,
                        help='X-axis limit for shortening/negative values (default: 200)')
    parser.add_argument('--exon-range', type=int, default=200,
                        help='X-axis limit for elongation/positive values (default: 200)')

    return parser.parse_args()


def parse_coordinate(coord_str):
    """
    Parse coordinate string 'chr:start-end' into (chr, start, end).

    Example: 'chr1:964167-964180' -> ('chr1', 964167, 964180)

    Returns None if coordinate is invalid or NA.
    """
    if not coord_str or coord_str == 'NA' or pd.isna(coord_str):
        return None
    try:
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


def calculate_splice_region_size(row, data_type):
    """
    Calculate splice region size with sign based on exon lengthening/shortening.

    For DONOR and ACCEPTOR events:
    - Compare intron_coord length vs reference junction length

    Sign determination:
    - if len(intron_coord) < len(reference) -> exon Lengthened (POSITIVE)
    - if len(intron_coord) > len(reference) -> exon Shortened (NEGATIVE)

    Parameters:
        row: DataFrame row
        data_type: 'poison' or 'productive'

    Returns:
        Signed integer representing splice region size in bp, or None
    """
    junction_type = row['junction_type']

    # Process both Donor and Acceptor events
    if junction_type not in ['poisonDonor', 'altDonor', 'poisonAcceptor', 'altAcceptor']:
        return None

    # Get splice region column based on data type
    if data_type == 'poison':
        splice_region_col = 'poison_event_splice_region'
        reference_col = 'productive_junctions'
    else:  # productive
        splice_region_col = 'productive_event_splice_region'
        reference_col = 'spanning_productive'

    # Get splice region value
    splice_region = row[splice_region_col]
    if not splice_region or splice_region == 'NA' or pd.isna(splice_region):
        return None

    # Handle comma-separated regions (take first one)
    if ',' in str(splice_region):
        splice_region = str(splice_region).split(',')[0]

    parsed_region = parse_coordinate(splice_region)
    if parsed_region is None:
        return None
    region_length = parsed_region[2] - parsed_region[1]  # stop - start

    # Calculate intron lengths for sign determination
    intron_coord = row['intron_coord']
    reference_junc = row[reference_col]

    intron_len = calculate_intron_length(intron_coord)

    # Handle comma-separated reference junctions (take first)
    if reference_junc and not pd.isna(reference_junc) and ',' in str(reference_junc):
        reference_junc = str(reference_junc).split(',')[0]
    reference_len = calculate_intron_length(reference_junc)

    if intron_len is None or reference_len is None:
        return None

    # Determine sign
    if intron_len < reference_len:
        # Exon Lengthened - positive
        sign = 1
    elif intron_len > reference_len:
        # Exon Shortened - negative
        sign = -1
    else:
        # No difference - return 0 (will be excluded)
        return 0

    return sign * region_length


def load_and_process_data(file_path, data_type, event_type, chromosome='All'):
    """
    Load TSV file and calculate splice region sizes.

    Parameters:
        file_path: Path to input file
        data_type: 'poison' or 'productive'
        event_type: 'donor' or 'acceptor'
        chromosome: Chromosome filter ('All' for no filter)

    Returns:
        pandas Series of signed splice region sizes (excluding zeros)
    """
    df = pd.read_csv(file_path, sep='\t')

    # Filter by chromosome if specified
    if chromosome != 'All':
        df = df[df['chrom'] == chromosome]

    # Filter by event type (donor or acceptor)
    if event_type == 'donor':
        if data_type == 'poison':
            df = df[df['junction_type'] == 'poisonDonor']
        else:
            df = df[df['junction_type'] == 'altDonor']
    else:  # acceptor
        if data_type == 'poison':
            df = df[df['junction_type'] == 'poisonAcceptor']
        else:
            df = df[df['junction_type'] == 'altAcceptor']

    # Calculate sizes
    sizes = df.apply(lambda row: calculate_splice_region_size(row, data_type), axis=1)

    # Remove None and 0 values
    sizes = sizes.dropna()
    sizes = sizes[sizes != 0]

    # Flip sign for acceptor events (elongation/shortening interpretation is reversed)
    if event_type == 'acceptor':
        sizes = -sizes

    return sizes


def create_histogram(sizes, title, output_path, intron_range, exon_range):
    """
    Create histogram with orange bars for multiples of 3, blue for non-multiples.

    Parameters:
        sizes: pandas Series of signed splice region sizes
        title: Plot title
        output_path: Where to save the figure
        intron_range: X-axis limit for negative values (shortenings)
        exon_range: X-axis limit for positive values (elongations)
    """
    if len(sizes) == 0:
        print(f"Warning: No data to plot for {title}")
        return

    fig, ax = plt.subplots(figsize=(12, 6))

    # Convert to integers for modulo operation
    sizes_int = sizes.astype(int)

    # Separate multiples of 3 vs non-multiples
    mult3 = sizes_int[sizes_int % 3 == 0]
    non_mult3 = sizes_int[sizes_int % 3 != 0]

    # Determine bin edges
    x_min = -intron_range
    x_max = exon_range
    bins = np.arange(x_min - 0.5, x_max + 1.5, 1)

    # Plot histograms
    ax.hist(non_mult3, bins=bins, color='steelblue', alpha=0.7,
            label=f'Non-multiple of 3', edgecolor='darkblue', linewidth=0.5)
    ax.hist(mult3, bins=bins, color='orange', alpha=0.7,
            label=f'Multiple of 3', edgecolor='darkorange', linewidth=0.5)

    # Add red dashed vertical line at 0
    ax.axvline(x=0, color='red', linestyle='--', linewidth=2, label='Zero')

    # Count positives and negatives
    n_positive = len(sizes_int[sizes_int > 0])

    # Formatting
    ax.set_xlabel('Region Size (bp)', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title(f'{title}\n(Orange = multiples of 3, Blue = non-multiples)', fontsize=12)
    ax.set_xlim(x_min, x_max)

    # Create custom legend
    legend_elements = [
        Patch(facecolor='orange', edgecolor='darkorange', alpha=0.7, label=f'Multiple of 3: {len(mult3)}'),
        Patch(facecolor='steelblue', edgecolor='darkblue', alpha=0.7, label=f'Non-multiple of 3: {len(non_mult3)}'),
        Line2D([0], [0], color='red', linestyle='--', linewidth=2, label='Zero'),
        Patch(facecolor='white', edgecolor='white', label=f'Pos: {n_positive}')
    ]
    ax.legend(handles=legend_elements, loc='upper right')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def create_proportion_histogram(sizes, title, output_path, intron_range, exon_range):
    """
    Create proportion-normalized histogram.
    """
    if len(sizes) == 0:
        print(f"Warning: No data to plot for {title}")
        return

    fig, ax = plt.subplots(figsize=(12, 6))

    # Convert to integers for modulo operation
    sizes_int = sizes.astype(int)

    mult3 = sizes_int[sizes_int % 3 == 0]
    non_mult3 = sizes_int[sizes_int % 3 != 0]

    x_min = -intron_range
    x_max = exon_range
    bins = np.arange(x_min - 0.5, x_max + 1.5, 1)

    # Use weights to normalize to proportions
    total = len(sizes_int)
    weights_mult3 = np.ones(len(mult3)) / total
    weights_non_mult3 = np.ones(len(non_mult3)) / total

    ax.hist(non_mult3, bins=bins, weights=weights_non_mult3, color='steelblue', alpha=0.7,
            label=f'Non-multiple of 3', edgecolor='darkblue', linewidth=0.5)
    ax.hist(mult3, bins=bins, weights=weights_mult3, color='orange', alpha=0.7,
            label=f'Multiple of 3', edgecolor='darkorange', linewidth=0.5)

    ax.axvline(x=0, color='red', linestyle='--', linewidth=2)

    ax.set_xlabel('Region Size (bp)', fontsize=12)
    ax.set_ylabel('Proportion', fontsize=12)
    ax.set_title(f'{title}\n(Orange = multiples of 3, Blue = non-multiples)', fontsize=12)
    ax.set_xlim(x_min, x_max)

    # Legend
    legend_elements = [
        Patch(facecolor='orange', edgecolor='darkorange', alpha=0.7,
              label=f'Multiple of 3: {len(mult3)/total*100:.1f}%'),
        Patch(facecolor='steelblue', edgecolor='darkblue', alpha=0.7,
              label=f'Non-multiple of 3: {len(non_mult3)/total*100:.1f}%'),
        Line2D([0], [0], color='red', linestyle='--', linewidth=2, label='Zero')
    ]
    ax.legend(handles=legend_elements, loc='upper right')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def create_overlay_histogram(poison_sizes, productive_sizes, title, output_path,
                              intron_range, exon_range, event_type, mode='counts'):
    """
    Create overlay histogram comparing poison and productive datasets.
    """
    fig, ax = plt.subplots(figsize=(12, 6))

    x_min = -intron_range
    x_max = exon_range
    bins = np.arange(x_min - 0.5, x_max + 1.5, 1)

    event_label = event_type.capitalize()

    if mode == 'proportion':
        weights_poison = np.ones(len(poison_sizes)) / len(poison_sizes) if len(poison_sizes) > 0 else None
        weights_prod = np.ones(len(productive_sizes)) / len(productive_sizes) if len(productive_sizes) > 0 else None
        ylabel = 'Proportion'
    else:
        weights_poison = None
        weights_prod = None
        ylabel = 'Count'

    if len(poison_sizes) > 0:
        ax.hist(poison_sizes, bins=bins, weights=weights_poison, color='crimson', alpha=0.5,
                label=f'Poison {event_label} (n={len(poison_sizes)})', edgecolor='darkred', linewidth=0.5)

    if len(productive_sizes) > 0:
        ax.hist(productive_sizes, bins=bins, weights=weights_prod, color='forestgreen', alpha=0.5,
                label=f'Productive {event_label} (n={len(productive_sizes)})', edgecolor='darkgreen', linewidth=0.5)

    ax.axvline(x=0, color='red', linestyle='--', linewidth=2, label='Zero')

    ax.set_xlabel('Region Size (bp)', fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=12)
    ax.set_xlim(x_min, x_max)
    ax.legend(loc='upper right')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def print_summary_statistics(data_dict):
    """
    Print comprehensive summary with symmetry and periodicity analysis.

    Parameters:
        data_dict: Dictionary with keys like 'Poison Donor', 'Productive Donor', etc.
    """
    print("=" * 70)
    print("SPLICE SITE DISTRIBUTION ANALYSIS SUMMARY")
    print("=" * 70)

    for name, sizes in data_dict.items():
        print(f"\n{name}:")
        print("-" * 40)

        if len(sizes) == 0:
            print("  No data available")
            continue

        sizes_int = sizes.astype(int)
        n_total = len(sizes_int)
        n_positive = len(sizes_int[sizes_int > 0])
        n_negative = len(sizes_int[sizes_int < 0])

        # Symmetry analysis
        elongation_ratio = n_positive / n_total * 100 if n_total > 0 else 0
        shortening_ratio = n_negative / n_total * 100 if n_total > 0 else 0

        print(f"  Total events: {n_total}")
        print(f"  Elongations (positive): {n_positive} ({elongation_ratio:.1f}%)")
        print(f"  Shortenings (negative): {n_negative} ({shortening_ratio:.1f}%)")
        if n_negative > 0:
            print(f"  Symmetry ratio (elong/short): {n_positive/n_negative:.2f}")
        else:
            print(f"  Symmetry ratio: N/A (no shortenings)")

        # Periodicity analysis (multiples of 3)
        mult3 = sizes_int[sizes_int % 3 == 0]
        non_mult3 = sizes_int[sizes_int % 3 != 0]
        mult3_pct = len(mult3) / n_total * 100 if n_total > 0 else 0

        print(f"\n  Periodicity Analysis:")
        print(f"    Multiples of 3: {len(mult3)} ({mult3_pct:.1f}%)")
        print(f"    Non-multiples of 3: {len(non_mult3)} ({100-mult3_pct:.1f}%)")
        print(f"    Expected by chance: 33.3%")

        # Basic statistics
        print(f"\n  Size Statistics:")
        print(f"    Mean: {sizes.mean():.1f} bp")
        print(f"    Median: {sizes.median():.1f} bp")
        print(f"    Std Dev: {sizes.std():.1f} bp")
        print(f"    Range: [{sizes.min():.0f}, {sizes.max():.0f}] bp")

    print("\n" + "=" * 70)


def main():
    args = parse_args()

    # Create output directories
    donor_counts_dir = os.path.join(args.output, 'donor', 'counts')
    donor_proportion_dir = os.path.join(args.output, 'donor', 'proportion')
    acceptor_counts_dir = os.path.join(args.output, 'acceptor', 'counts')
    acceptor_proportion_dir = os.path.join(args.output, 'acceptor', 'proportion')

    for d in [donor_counts_dir, donor_proportion_dir, acceptor_counts_dir, acceptor_proportion_dir]:
        os.makedirs(d, exist_ok=True)

    # Load and process DONOR data
    print("Loading poison donor data...")
    poison_donor_sizes = load_and_process_data(args.poison, 'poison', 'donor', args.chromosome)
    print(f"  Loaded {len(poison_donor_sizes)} poison donor events")

    print("Loading productive donor data...")
    productive_donor_sizes = load_and_process_data(args.productive, 'productive', 'donor', args.chromosome)
    print(f"  Loaded {len(productive_donor_sizes)} productive donor events")

    # Load and process ACCEPTOR data
    print("Loading poison acceptor data...")
    poison_acceptor_sizes = load_and_process_data(args.poison, 'poison', 'acceptor', args.chromosome)
    print(f"  Loaded {len(poison_acceptor_sizes)} poison acceptor events")

    print("Loading productive acceptor data...")
    productive_acceptor_sizes = load_and_process_data(args.productive, 'productive', 'acceptor', args.chromosome)
    print(f"  Loaded {len(productive_acceptor_sizes)} productive acceptor events")

    # Generate suffix for chromosome filtering
    suffix = f"_{args.chromosome}" if args.chromosome != 'All' else ""

    # ===== DONOR PLOTS =====
    print("\nGenerating donor count histograms...")
    create_histogram(poison_donor_sizes,
                    'Distribution of poison_event_splice_region Sizes (poisonDonor only)',
                    os.path.join(donor_counts_dir, f'poison_histogram{suffix}.png'),
                    args.intron_range, args.exon_range)

    create_histogram(productive_donor_sizes,
                    'Distribution of productive_event_splice_region Sizes (altDonor only)',
                    os.path.join(donor_counts_dir, f'productive_histogram{suffix}.png'),
                    args.intron_range, args.exon_range)

    create_overlay_histogram(poison_donor_sizes, productive_donor_sizes,
                            'Poison vs Productive Donor Comparison',
                            os.path.join(donor_counts_dir, f'overlay_comparison{suffix}.png'),
                            args.intron_range, args.exon_range, 'donor', mode='counts')

    print("Generating donor proportion histograms...")
    create_proportion_histogram(poison_donor_sizes,
                               'Distribution of poison_event_splice_region Sizes (poisonDonor only)',
                               os.path.join(donor_proportion_dir, f'poison_histogram{suffix}.png'),
                               args.intron_range, args.exon_range)

    create_proportion_histogram(productive_donor_sizes,
                               'Distribution of productive_event_splice_region Sizes (altDonor only)',
                               os.path.join(donor_proportion_dir, f'productive_histogram{suffix}.png'),
                               args.intron_range, args.exon_range)

    create_overlay_histogram(poison_donor_sizes, productive_donor_sizes,
                            'Poison vs Productive Donor Comparison (Proportion)',
                            os.path.join(donor_proportion_dir, f'overlay_comparison{suffix}.png'),
                            args.intron_range, args.exon_range, 'donor', mode='proportion')

    # ===== ACCEPTOR PLOTS =====
    print("\nGenerating acceptor count histograms...")
    create_histogram(poison_acceptor_sizes,
                    'Distribution of poison_event_splice_region Sizes (poisonAcceptor only)',
                    os.path.join(acceptor_counts_dir, f'poison_histogram{suffix}.png'),
                    args.intron_range, args.exon_range)

    create_histogram(productive_acceptor_sizes,
                    'Distribution of productive_event_splice_region Sizes (altAcceptor only)',
                    os.path.join(acceptor_counts_dir, f'productive_histogram{suffix}.png'),
                    args.intron_range, args.exon_range)

    create_overlay_histogram(poison_acceptor_sizes, productive_acceptor_sizes,
                            'Poison vs Productive Acceptor Comparison',
                            os.path.join(acceptor_counts_dir, f'overlay_comparison{suffix}.png'),
                            args.intron_range, args.exon_range, 'acceptor', mode='counts')

    print("Generating acceptor proportion histograms...")
    create_proportion_histogram(poison_acceptor_sizes,
                               'Distribution of poison_event_splice_region Sizes (poisonAcceptor only)',
                               os.path.join(acceptor_proportion_dir, f'poison_histogram{suffix}.png'),
                               args.intron_range, args.exon_range)

    create_proportion_histogram(productive_acceptor_sizes,
                               'Distribution of productive_event_splice_region Sizes (altAcceptor only)',
                               os.path.join(acceptor_proportion_dir, f'productive_histogram{suffix}.png'),
                               args.intron_range, args.exon_range)

    create_overlay_histogram(poison_acceptor_sizes, productive_acceptor_sizes,
                            'Poison vs Productive Acceptor Comparison (Proportion)',
                            os.path.join(acceptor_proportion_dir, f'overlay_comparison{suffix}.png'),
                            args.intron_range, args.exon_range, 'acceptor', mode='proportion')

    # Print summary statistics
    summary_data = {
        'Poison Donor': poison_donor_sizes,
        'Productive Donor': productive_donor_sizes,
        'Poison Acceptor': poison_acceptor_sizes,
        'Productive Acceptor': productive_acceptor_sizes
    }
    print_summary_statistics(summary_data)

    print("\nAnalysis complete!")
    print(f"  Donor plots saved to: {os.path.join(args.output, 'donor')}")
    print(f"  Acceptor plots saved to: {os.path.join(args.output, 'acceptor')}")


if __name__ == '__main__':
    main()
