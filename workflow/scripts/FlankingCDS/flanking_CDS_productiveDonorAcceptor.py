#!/usr/bin/env python3
"""
Identify flanking CDS features for productive alternative donor/acceptor junctions.

Finds CDS features at the donor and acceptor splice sites of spanning productive
junctions that are associated with alternative donor/acceptor events. Annotates
which CDS is shared (unchanged) vs. alternative between spanning productive and
alternative junctions.

Author: Dylan Stermer (with Claude Code)
Date: February 2026
"""

import pandas as pd
import argparse
from collections import defaultdict
import sys

def parse_gtf_attributes(attr_string):
    """Parse GTF attributes string into dictionary."""
    attrs = {}
    for item in attr_string.strip().split(';'):
        item = item.strip()
        if not item:
            continue
        parts = item.split(' ', 1)
        if len(parts) == 2:
            key = parts[0]
            value = parts[1].strip('"')
            attrs[key] = value
    return attrs

def load_gtf_cds(gtf_file, chromosome=None):
    """
    Load CDS features from GTF, indexed by chromosome.

    Args:
        gtf_file: Path to GTF file
        chromosome: Optional chromosome filter (e.g., 'chr1')

    Returns:
        cds_index: Dict {chrom: [(start, end, strand, gene_name, transcript_id)]}
    """

    print(f"Loading GTF: {gtf_file}")
    if chromosome:
        print(f"  Filtering for chromosome: {chromosome}")

    cds_features = defaultdict(list)

    line_count = 0
    cds_count = 0

    with open(gtf_file, 'r') as f:
        for line in f:
            line_count += 1

            if line_count % 1000000 == 0:
                print(f"  Processed {line_count:,} lines...")

            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chrom = fields[0]

            # Skip if chromosome filter is active and doesn't match
            if chromosome and chrom != chromosome:
                continue

            feature = fields[2]
            start = int(fields[3])  # GTF is 1-based inclusive
            end = int(fields[4])    # GTF is 1-based inclusive
            strand = fields[6]
            attributes = parse_gtf_attributes(fields[8])

            gene_name = attributes.get('gene_name', 'NA')
            transcript_id = attributes.get('transcript_id', 'NA')

            if feature == 'CDS':
                cds_features[chrom].append((start, end, strand, gene_name, transcript_id))
                cds_count += 1

    print(f"  Total lines: {line_count:,}")
    print(f"  CDS features: {cds_count:,}")

    return cds_features

def parse_coordinate(coord_str):
    """
    Parse coordinate string 'chr:start-end' into (chr, start, end).

    Returns coordinates as they appear in the string (BED format).
    """
    chrom, pos = coord_str.split(':')
    start, end = pos.split('-')
    return chrom, int(start), int(end)

def filter_spanning_productive_junctions(spanning_productive_str, poison_start, poison_end,
                                         junction_type, strand, filter_mode='closest'):
    """
    Filter spanning productive junctions to keep only the most relevant ones.

    Args:
        spanning_productive_str: Comma-separated spanning productive junction coordinates
        poison_start: Alternative junction start coordinate
        poison_end: Alternative junction end coordinate
        junction_type: 'altDonor' or 'altAcceptor'
        strand: '+' or '-'
        filter_mode: 'closest', 'threshold', or 'none'

    Returns:
        Filtered comma-separated spanning productive junction string
    """
    if pd.isna(spanning_productive_str) or spanning_productive_str == '':
        return ''

    productive_junctions = str(spanning_productive_str).split(',')

    if filter_mode == 'none' or len(productive_junctions) <= 1:
        return spanning_productive_str

    # Parse each productive junction and calculate distance
    junction_distances = []
    for prod_junc in productive_junctions:
        prod_junc = prod_junc.strip()
        _, prod_start, prod_end = parse_coordinate(prod_junc)

        # Strand-aware distance calculation
        if junction_type == 'altAcceptor':
            # Plus strand: acceptor is at end position
            # Minus strand: acceptor is at start position (biologically)
            if strand == '+':
                distance = abs(prod_end - poison_end)
            else:  # strand == '-'
                distance = abs(prod_start - poison_start)
        elif junction_type == 'altDonor':
            # Plus strand: donor is at start position
            # Minus strand: donor is at end position (biologically)
            if strand == '+':
                distance = abs(prod_start - poison_start)
            else:  # strand == '-'
                distance = abs(prod_end - poison_end)
        else:
            distance = float('inf')

        junction_distances.append((prod_junc, distance))

    if filter_mode == 'closest':
        # Keep only the closest junction(s) - may be multiple with same distance
        min_distance = min(d for _, d in junction_distances)
        filtered = [junc for junc, dist in junction_distances if dist == min_distance]
        return ','.join(filtered)

    elif filter_mode == 'threshold':
        # Keep junctions within threshold
        MAX_DISTANCE = 500
        filtered = [junc for junc, dist in junction_distances if dist <= MAX_DISTANCE]
        return ','.join(filtered) if filtered else ''

    return spanning_productive_str

def find_cds_at_junction_boundary(chrom, junction_start, junction_end, strand,
                                  boundary_type, cds_index, tolerance=0):
    """
    Find CDS features at a junction boundary.

    Args:
        chrom: Chromosome
        junction_start: Junction start (BED format, 0-based)
        junction_end: Junction end (BED format)
        strand: Strand (+ or -)
        boundary_type: 'donor' or 'acceptor'
        cds_index: Dict of CDS features
        tolerance: Allowed coordinate deviation in bp

    Returns:
        List of (start, end, gene_name, transcript_id) tuples
    """

    if chrom not in cds_index:
        return []

    matching_cds = []

    # Determine expected CDS boundary position based on strand and boundary type
    if strand == '+':
        if boundary_type == 'donor':
            # Donor CDS (upstream): should END at junction_start (GTF 1-based)
            target_pos = junction_start
            for cds_start, cds_end, cds_strand, gene_name, transcript_id in cds_index[chrom]:
                if cds_strand != strand:
                    continue
                if abs(cds_end - target_pos) <= tolerance:
                    matching_cds.append((cds_start, cds_end, gene_name, transcript_id))

        elif boundary_type == 'acceptor':
            # Acceptor CDS (downstream): should START at junction_end + 1 (GTF 1-based)
            target_pos = junction_end + 1
            for cds_start, cds_end, cds_strand, gene_name, transcript_id in cds_index[chrom]:
                if cds_strand != strand:
                    continue
                if abs(cds_start - target_pos) <= tolerance:
                    matching_cds.append((cds_start, cds_end, gene_name, transcript_id))

    else:  # Minus strand
        if boundary_type == 'donor':
            # Donor CDS (biologically upstream): should START at junction_end + 1
            target_pos = junction_end + 1
            for cds_start, cds_end, cds_strand, gene_name, transcript_id in cds_index[chrom]:
                if cds_strand != strand:
                    continue
                if abs(cds_start - target_pos) <= tolerance:
                    matching_cds.append((cds_start, cds_end, gene_name, transcript_id))

        elif boundary_type == 'acceptor':
            # Acceptor CDS (biologically downstream): should END at junction_start
            target_pos = junction_start
            for cds_start, cds_end, cds_strand, gene_name, transcript_id in cds_index[chrom]:
                if cds_strand != strand:
                    continue
                if abs(cds_end - target_pos) <= tolerance:
                    matching_cds.append((cds_start, cds_end, gene_name, transcript_id))

    return matching_cds

def process_productive_donor_acceptor(productive_file, cds_index, tolerance=0, chromosome=None, filter_mode='closest'):
    """
    Process productive alternative donor/acceptor junctions and find flanking CDS features.

    Args:
        productive_file: Productive alternative donor/acceptor BED file
        cds_index: Dict of CDS features
        tolerance: Coordinate matching tolerance
        chromosome: Optional chromosome filter
        filter_mode: Spanning productive junction filtering mode ('closest', 'threshold', or 'none')

    Returns:
        (productive_with_flanks, flanking_cds_list)
    """

    print(f"\nReading productive alternative donor/acceptor junctions: {productive_file}")
    pd_df = pd.read_csv(productive_file, sep='\t')

    # Handle column name with or without '#' prefix
    chr_col = '#chr' if '#chr' in pd_df.columns else 'chr'

    # Filter by chromosome if specified
    if chromosome:
        print(f"  Filtering for chromosome: {chromosome}")
        pd_df = pd_df[pd_df[chr_col] == chromosome].copy()

    print(f"  Total junctions to process: {len(pd_df)}")

    results = []
    flanking_cds_features = []

    found_count = 0
    not_found_count = 0

    for idx, row in pd_df.iterrows():
        if idx % 1000 == 0 and idx > 0:
            print(f"  Processed {idx:,}/{len(pd_df):,} junctions...")

        chrom = row[chr_col]
        start = row['start']
        end = row['end']
        gene_name = row['name']
        strand = row['strand']
        junction_type = row['junction_type']
        intron_coord = row['intron_coord']
        spanning_productive_str = row['spanning_productive']

        # Skip if no spanning productive junctions (NaN or empty)
        if pd.isna(spanning_productive_str) or spanning_productive_str == '':
            not_found_count += 1
            continue

        # Filter spanning productive junctions to keep only relevant ones
        spanning_productive_filtered = filter_spanning_productive_junctions(
            spanning_productive_str, start, end, junction_type, strand, filter_mode
        )

        if not spanning_productive_filtered:
            not_found_count += 1
            continue

        # Parse spanning productive junctions (may be comma-separated)
        spanning_productive_list = str(spanning_productive_filtered).split(',')

        all_donor_cds = []
        all_acceptor_cds = []

        # Process each spanning productive junction
        for prod_junc in spanning_productive_list:
            prod_junc = prod_junc.strip()

            # Parse productive junction coordinates
            _, prod_start, prod_end = parse_coordinate(prod_junc)

            # Find donor and acceptor CDS for this productive junction
            donor_cds = find_cds_at_junction_boundary(
                chrom, prod_start, prod_end, strand, 'donor', cds_index, tolerance
            )
            acceptor_cds = find_cds_at_junction_boundary(
                chrom, prod_start, prod_end, strand, 'acceptor', cds_index, tolerance
            )

            all_donor_cds.extend(donor_cds)
            all_acceptor_cds.extend(acceptor_cds)

        # Deduplicate CDS lists
        all_donor_cds = list(set(all_donor_cds))
        all_acceptor_cds = list(set(all_acceptor_cds))

        # Determine shared vs. alternative based on junction type
        if junction_type == 'altDonor':
            shared_type = 'acceptor_CDS'
            alternative_type = 'donor_CDS'
        elif junction_type == 'altAcceptor':
            shared_type = 'donor_CDS'
            alternative_type = 'acceptor_CDS'
        else:
            shared_type = 'NA'
            alternative_type = 'NA'

        # Format CDS coordinates as comma-separated lists (0-based BED format)
        # Use dict to deduplicate by coordinate while preserving one instance
        donor_coords_dict = {f"{chrom}:{cds[0]-1}-{cds[1]}": None for cds in all_donor_cds}
        acceptor_coords_dict = {f"{chrom}:{cds[0]-1}-{cds[1]}": None for cds in all_acceptor_cds}

        donor_cds_str = ",".join(sorted(donor_coords_dict.keys())) if donor_coords_dict else "NA"
        acceptor_cds_str = ",".join(sorted(acceptor_coords_dict.keys())) if acceptor_coords_dict else "NA"

        # Deduplicate spanning productive junctions as well
        spanning_productive_unique = list(dict.fromkeys([pj.strip() for pj in spanning_productive_list]))
        spanning_productive_str_unique = ",".join(spanning_productive_unique)

        # Alternative junction identifier
        pj_id = f"{chrom}:{start}-{end}_{gene_name}"

        # Calculate productive event exon regions
        # On + strand: donor is at start, acceptor is at end
        # On - strand: donor is at end (biologically), acceptor is at start (biologically)

        if junction_type == 'altAcceptor':
            splice_regions = []

            if strand == '+':
                # Alternative acceptor is at end position
                alt_acceptor_pos = end
                for prod_junc in spanning_productive_unique:
                    _, prod_start, prod_end = parse_coordinate(prod_junc.strip())
                    # Splice region: from spanning productive acceptor to alternative acceptor
                    if prod_end < alt_acceptor_pos:
                        splice_regions.append(f"{chrom}:{prod_end}-{alt_acceptor_pos}")
                    elif prod_end > alt_acceptor_pos:
                        splice_regions.append(f"{chrom}:{alt_acceptor_pos}-{prod_end}")
                    else:
                        splice_regions.append(f"{chrom}:{prod_end}-{alt_acceptor_pos}")

                # Full exon: from alternative acceptor to acceptor CDS end
                if acceptor_cds_str != "NA":
                    first_cds = acceptor_cds_str.split(',')[0].strip()
                    _, cds_start, cds_end = parse_coordinate(first_cds)
                    productive_full_exon = f"{chrom}:{alt_acceptor_pos}-{cds_end}"
                else:
                    productive_full_exon = "NA"

            else:  # strand == '-'
                # Alternative acceptor is at start position (biologically)
                alt_acceptor_pos = start
                for prod_junc in spanning_productive_unique:
                    _, prod_start, prod_end = parse_coordinate(prod_junc.strip())
                    # Splice region: from spanning productive acceptor to alternative acceptor
                    if prod_start < alt_acceptor_pos:
                        splice_regions.append(f"{chrom}:{prod_start}-{alt_acceptor_pos}")
                    elif prod_start > alt_acceptor_pos:
                        splice_regions.append(f"{chrom}:{alt_acceptor_pos}-{prod_start}")
                    else:
                        splice_regions.append(f"{chrom}:{prod_start}-{alt_acceptor_pos}")

                # Full exon: from acceptor CDS start to alternative acceptor
                if acceptor_cds_str != "NA":
                    first_cds = acceptor_cds_str.split(',')[0].strip()
                    _, cds_start, cds_end = parse_coordinate(first_cds)
                    productive_full_exon = f"{chrom}:{cds_start}-{alt_acceptor_pos}"
                else:
                    productive_full_exon = "NA"

            productive_splice_region = ",".join(splice_regions) if splice_regions else "NA"

        elif junction_type == 'altDonor':
            splice_regions = []

            if strand == '+':
                # Alternative donor is at start position
                alt_donor_pos = start
                for prod_junc in spanning_productive_unique:
                    _, prod_start, prod_end = parse_coordinate(prod_junc.strip())
                    # Splice region: from alternative donor to spanning productive donor
                    if alt_donor_pos < prod_start:
                        splice_regions.append(f"{chrom}:{alt_donor_pos}-{prod_start}")
                    elif alt_donor_pos > prod_start:
                        splice_regions.append(f"{chrom}:{prod_start}-{alt_donor_pos}")
                    else:
                        splice_regions.append(f"{chrom}:{alt_donor_pos}-{prod_start}")

                # Full exon: from donor CDS start to alternative donor
                if donor_cds_str != "NA":
                    first_cds = donor_cds_str.split(',')[0].strip()
                    _, cds_start, cds_end = parse_coordinate(first_cds)
                    productive_full_exon = f"{chrom}:{cds_start}-{alt_donor_pos}"
                else:
                    productive_full_exon = "NA"

            else:  # strand == '-'
                # Alternative donor is at end position (biologically)
                alt_donor_pos = end
                for prod_junc in spanning_productive_unique:
                    _, prod_start, prod_end = parse_coordinate(prod_junc.strip())
                    # Splice region: from alternative donor to spanning productive donor
                    if alt_donor_pos < prod_end:
                        splice_regions.append(f"{chrom}:{alt_donor_pos}-{prod_end}")
                    elif alt_donor_pos > prod_end:
                        splice_regions.append(f"{chrom}:{prod_end}-{alt_donor_pos}")
                    else:
                        splice_regions.append(f"{chrom}:{alt_donor_pos}-{prod_end}")

                # Full exon: from alternative donor to donor CDS end
                if donor_cds_str != "NA":
                    first_cds = donor_cds_str.split(',')[0].strip()
                    _, cds_start, cds_end = parse_coordinate(first_cds)
                    productive_full_exon = f"{chrom}:{alt_donor_pos}-{cds_end}"
                else:
                    productive_full_exon = "NA"

            productive_splice_region = ",".join(splice_regions) if splice_regions else "NA"

        else:
            productive_splice_region = "NA"
            productive_full_exon = "NA"

        results.append({
            'chrom': chrom,
            'start': start,
            'end': end,
            'name': gene_name,
            'length': end - start,
            'strand': strand,
            'junction_type': junction_type,
            'intron_coord': intron_coord,
            'spanning_productive': spanning_productive_str_unique,
            'donor_CDS': donor_cds_str,
            'acceptor_CDS': acceptor_cds_str,
            'shared_CDS': shared_type,
            'alternative_CDS': alternative_type,
            'productive_event_splice_region': productive_splice_region,
            'productive_event_CDS_region': productive_full_exon
        })

        # Add ALL flanking CDS to separate list
        for cds in all_donor_cds:
            for prod_junc in spanning_productive_list:
                flanking_cds_features.append({
                    'chrom': chrom,
                    'start': cds[0] - 1,  # Convert to 0-based BED
                    'end': cds[1],
                    'name': f"{chrom}:{cds[0]-1}-{cds[1]}_{cds[2]}",
                    'length': cds[1] - cds[0] + 1,
                    'strand': strand,
                    'position': 'donor',
                    'splice_site_status': 'shared' if shared_type == 'donor_CDS' else 'alternative',
                    'altJunction': pj_id,
                    'junction_type': junction_type,
                    'intron_coord': intron_coord,
                    'spanning_productive_junction': prod_junc.strip(),
                    'transcript_id': cds[3]
                })

        for cds in all_acceptor_cds:
            for prod_junc in spanning_productive_list:
                flanking_cds_features.append({
                    'chrom': chrom,
                    'start': cds[0] - 1,  # Convert to 0-based BED
                    'end': cds[1],
                    'name': f"{chrom}:{cds[0]-1}-{cds[1]}_{cds[2]}",
                    'length': cds[1] - cds[0] + 1,
                    'strand': strand,
                    'position': 'acceptor',
                    'splice_site_status': 'shared' if shared_type == 'acceptor_CDS' else 'alternative',
                    'altJunction': pj_id,
                    'junction_type': junction_type,
                    'intron_coord': intron_coord,
                    'spanning_productive_junction': prod_junc.strip(),
                    'transcript_id': cds[3]
                })

        if all_donor_cds or all_acceptor_cds:
            found_count += 1
        else:
            not_found_count += 1

    print(f"\n  Found flanking CDS for: {found_count:,} junctions")
    print(f"  No CDS found for: {not_found_count:,} junctions")

    return pd.DataFrame(results), pd.DataFrame(flanking_cds_features)

def main():
    parser = argparse.ArgumentParser(
        description='Identify flanking CDS features for productive alternative donor/acceptor junctions',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python flanking_CDS_productiveDonorAcceptor.py \\
      -g Reference.gtf \\
      -p _altDonorAcceptor.bed \\
      -o output/productive_da_flanking.tsv \\
      -f output/flanking_cds.tsv \\
      --tolerance 0 \\
      --chromosome chr1

Flanking CDS identification:
  - For each spanning productive junction associated with an alternative donor/acceptor event
  - Finds CDS at the donor (5' splice site) and acceptor (3' splice site) boundaries
  - Annotates which CDS is shared (unchanged) vs. alternative
  - Alternative Donor: acceptor is shared, donor is alternative
  - Alternative Acceptor: donor is shared, acceptor is alternative

Coordinate systems:
  - Input BED uses 0-based start
  - GTF uses 1-based coordinates
  - Output in 0-based BED format

Output format:
  - All coordinates in 0-based BED format (start position is -1 from 1-based)
        """
    )

    parser.add_argument('-g', '--gtf', required=True,
                       help='Reference GTF file')
    parser.add_argument('-p', '--productive-bed', required=True,
                       help='Productive alternative donor/acceptor BED file')
    parser.add_argument('-o', '--output', required=True,
                       help='Output: junctions with flanking CDS info')
    parser.add_argument('-f', '--flanking-output', required=True,
                       help='Output: flanking CDS features only')
    parser.add_argument('--tolerance', type=int, default=0,
                       help='Coordinate matching tolerance in bp (default: 0)')
    parser.add_argument('--chromosome', type=str, default=None,
                       help='Process only this chromosome (optional, for testing)')
    parser.add_argument('--filter-mode', type=str, default='closest',
                       choices=['closest', 'threshold', 'none'],
                       help='Spanning productive junction filtering mode (default: closest)')

    args = parser.parse_args()

    print("="*70)
    print("IDENTIFY FLANKING CDS FOR PRODUCTIVE ALTERNATIVE DONOR/ACCEPTOR JUNCTIONS")
    print("="*70)
    print(f"GTF: {args.gtf}")
    print(f"Productive donor/acceptor BED: {args.productive_bed}")
    print(f"Output 1: {args.output}")
    print(f"Output 2: {args.flanking_output}")
    print(f"Tolerance: {args.tolerance} bp")
    print(f"Filter mode: {args.filter_mode}")
    if args.chromosome:
        print(f"Chromosome filter: {args.chromosome}")
    print("="*70 + "\n")

    # Load GTF CDS features
    cds_index = load_gtf_cds(args.gtf, chromosome=args.chromosome)

    # Process productive alternative donor/acceptor junctions
    productive_with_flanks, flanking_cds = process_productive_donor_acceptor(
        args.productive_bed, cds_index, tolerance=args.tolerance, chromosome=args.chromosome,
        filter_mode=args.filter_mode
    )

    # Save outputs
    print(f"\nSaving outputs...")

    # Remove duplicates
    print(f"  Removing duplicate rows...")
    productive_before = len(productive_with_flanks)
    flanking_before = len(flanking_cds)

    productive_with_flanks = productive_with_flanks.drop_duplicates()
    flanking_cds = flanking_cds.drop_duplicates()

    productive_removed = productive_before - len(productive_with_flanks)
    flanking_removed = flanking_before - len(flanking_cds)

    if productive_removed > 0:
        print(f"  Removed {productive_removed:,} duplicate junction rows")
    else:
        print(f"  No duplicate junction rows found")

    if flanking_removed > 0:
        print(f"  Removed {flanking_removed:,} duplicate flanking CDS rows")
    else:
        print(f"  No duplicate flanking CDS rows found")

    # Output 1: Junctions with flanking CDS info
    productive_with_flanks.to_csv(args.output, sep='\t', index=False)
    print(f"  Saved {len(productive_with_flanks):,} unique junctions to: {args.output}")

    # Output 2: Flanking CDS only
    flanking_cds.to_csv(args.flanking_output, sep='\t', index=False)
    print(f"  Saved {len(flanking_cds):,} unique flanking CDS features to: {args.flanking_output}")

    print("\n" + "="*70)
    print("COMPLETE")
    print("="*70)

if __name__ == "__main__":
    main()
