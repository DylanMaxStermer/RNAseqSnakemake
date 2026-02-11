# #!/usr/bin/env python3
# """
# Identify flanking CDS features for poison exons using intron coordinates.

# Finds CDS features that are immediately adjacent to the upstream and downstream
# introns, which represent the flanking exons.

# Author: Dylan Stermer
# Date: December 2024
# """

# import pandas as pd
# import argparse
# from collections import defaultdict
# import sys

# def parse_gtf_attributes(attr_string):

#     """Parse GTF attributes string into dictionary."""
#     attrs = {}
#     for item in attr_string.strip().split(';'):
#         item = item.strip()
#         if not item:
#             continue
#         parts = item.split(' ', 1)
#         if len(parts) == 2:
#             key = parts[0]
#             value = parts[1].strip('"')
#             attrs[key] = value
#     return attrs

# def load_gtf_cds(gtf_file):
#     """
#     Load CDS features from GTF, indexed by chromosome.
    
#     Returns:
#         cds_index: Dict {chrom: [(start, end, strand, gene_name, transcript_id)]}
#     """
    
#     print(f"Loading GTF: {gtf_file}")
    
#     cds_features = defaultdict(list)
    
#     line_count = 0
#     cds_count = 0
    
#     with open(gtf_file, 'r') as f:
#         for line in f:
#             line_count += 1
            
#             if line_count % 1000000 == 0:
#                 print(f"  Processed {line_count:,} lines...")
            
#             if line.startswith('#'):
#                 continue
            
#             fields = line.strip().split('\t')
#             if len(fields) < 9:
#                 continue
            
#             chrom = fields[0]
#             feature = fields[2]
#             start = int(fields[3])  # GTF is 1-based inclusive
#             end = int(fields[4])    # GTF is 1-based inclusive
#             strand = fields[6]
#             attributes = parse_gtf_attributes(fields[8])
            
#             gene_name = attributes.get('gene_name', 'NA')
#             transcript_id = attributes.get('transcript_id', 'NA')
            
#             if feature == 'CDS':
#                 cds_features[chrom].append((start, end, strand, gene_name, transcript_id))
#                 cds_count += 1
    
#     print(f"  Total lines: {line_count:,}")
#     print(f"  CDS features: {cds_count:,}")
    
#     return cds_features

# def parse_coordinate(coord_str):
#     """
#     Parse coordinate string 'chr:start-end' into (chr, start, end).
    
#     Returns 1-based coordinates matching GTF.
#     """
#     chrom, pos = coord_str.split(':')
#     start, end = pos.split('-')
#     return chrom, int(start), int(end)

# def find_overlapping_cds(chrom, start, end, cds_index):
#     """
#     Find CDS features that overlap with given coordinates.
    
#     Args:
#         chrom: Chromosome
#         start, end: 1-based coordinates
#         cds_index: Dict of CDS features
    
#     Returns:
#         List of (start, end, strand, gene_name, transcript_id)
#     """
#     if chrom not in cds_index:
#         return []
    
#     overlapping = []
#     for cds_start, cds_end, strand, gene_name, transcript_id in cds_index[chrom]:
#         # Check for overlap
#         if not (end < cds_start or start > cds_end):
#             overlapping.append((cds_start, cds_end, strand, gene_name, transcript_id))
    
#     return overlapping

# def find_flanking_cds_by_introns(chrom, upstream_intron, downstream_intron, strand, cds_index):
#     """
#     Find CDS features that are adjacent to the introns (flanking the poison exon).
    
#     Returns ALL matching CDS features, not just one.
    
#     Returns:
#         (upstream_cds_list, downstream_cds_list): Each as list of (start, end, gene_name, transcript_id) tuples
#     """
    
#     if chrom not in cds_index:
#         return [], []
    
#     # Parse intron coordinates
#     _, up_intron_start, up_intron_end = parse_coordinate(upstream_intron)
#     _, down_intron_start, down_intron_end = parse_coordinate(downstream_intron)
    
#     upstream_cds_list = []
#     downstream_cds_list = []
    
#     if strand == '+':
#         # Upstream CDS: should end at up_intron_start - 1
#         for cds_start, cds_end, cds_strand, gene_name, transcript_id in cds_index[chrom]:
#             if cds_strand != strand:
#                 continue
            
#             if cds_end == up_intron_start - 1 or cds_end == up_intron_start or cds_end == up_intron_start - 2:
#                 upstream_cds_list.append((cds_start, cds_end, gene_name, transcript_id))
        
#         # Downstream CDS: should start at down_intron_end + 1
#         for cds_start, cds_end, cds_strand, gene_name, transcript_id in cds_index[chrom]:
#             if cds_strand != strand:
#                 continue
            
#             if cds_start == down_intron_end + 1 or cds_start == down_intron_end or cds_start == down_intron_end + 2:
#                 downstream_cds_list.append((cds_start, cds_end, gene_name, transcript_id))
    
#     else:  # Minus strand
#         # Upstream CDS (5', biologically): should start at down_intron_end + 1
#         for cds_start, cds_end, cds_strand, gene_name, transcript_id in cds_index[chrom]:
#             if cds_strand != strand:
#                 continue
            
#             if cds_start == down_intron_end + 1 or cds_start == down_intron_end or cds_start == down_intron_end + 2:
#                 upstream_cds_list.append((cds_start, cds_end, gene_name, transcript_id))
        
#         # Downstream CDS (3', biologically): should end at up_intron_start - 1
#         for cds_start, cds_end, cds_strand, gene_name, transcript_id in cds_index[chrom]:
#             if cds_strand != strand:
#                 continue
            
#             if cds_end == up_intron_start - 1 or cds_end == up_intron_start or cds_end == up_intron_start - 2:
#                 downstream_cds_list.append((cds_start, cds_end, gene_name, transcript_id))
    
#     return upstream_cds_list, downstream_cds_list

# def process_poison_exons(pe_file, cds_index):
#     """
#     Process poison exons and find flanking CDS features.
    
#     Returns:
#         (pe_with_flanks, flanking_cds_list)
#     """
    
#     print(f"\nReading poison exons: {pe_file}")
#     pe_df = pd.read_csv(pe_file, sep='\t')
    
#     print(f"  Total poison exons: {len(pe_df)}")
    
#     results = []
#     flanking_cds_features = []
    
#     found_count = 0
#     not_found_count = 0
    
#     for idx, row in pe_df.iterrows():
#         if idx % 1000 == 0 and idx > 0:
#             print(f"  Processed {idx:,}/{len(pe_df):,} poison exons...")
        
#         coding_junction = row['coding_junction']
#         upstream_intron = row['upstream_intron']
#         downstream_intron = row['downstream_intron']
        
#         # Parse intron coordinates to determine actual poison exon boundaries
#         up_chrom, up_start, up_end = parse_coordinate(upstream_intron)
#         down_chrom, down_start, down_end = parse_coordinate(downstream_intron)
        
#         chrom = up_chrom
        
#         # Actual poison exon coordinates (1-based, inclusive)
#         pe_actual_start = up_end + 1
#         pe_actual_end = down_start
        
#         # First, try to find flanking CDS to get strand/gene information
#         strand = None
#         gene_name = None
#         upstream_cds_list = []
#         downstream_cds_list = []
        
#         # Try plus strand first
#         for test_strand in ['+', '-']:
#             up_cds_list, down_cds_list = find_flanking_cds_by_introns(
#                 chrom, upstream_intron, downstream_intron, test_strand, cds_index
#             )
            
#             if up_cds_list or down_cds_list:
#                 # Found flanking CDS, get gene name from them
#                 strand = test_strand
#                 upstream_cds_list = up_cds_list
#                 downstream_cds_list = down_cds_list
                
#                 # Get gene name from whichever CDS we found
#                 if up_cds_list:
#                     gene_name = up_cds_list[0][2]
#                 elif down_cds_list:
#                     gene_name = down_cds_list[0][2]
                
#                 break
        
#         # If we didn't find flanking CDS, try overlapping the poison exon itself
#         if strand is None:
#             overlapping_cds = find_overlapping_cds(chrom, pe_actual_start, pe_actual_end, cds_index)
#             if overlapping_cds:
#                 cds_start, cds_end, strand, gene_name, transcript_id = overlapping_cds[0]
#             else:
#                 not_found_count += 1
#                 pe_id = f"{chrom}:{pe_actual_start}-{pe_actual_end}_NA"
#                 results.append({
#                     'chrom': chrom,
#                     'start': pe_actual_start,
#                     'end': pe_actual_end,
#                     'name': pe_id,
#                     'length': pe_actual_end - pe_actual_start + 1,
#                     'strand': 'NA',
#                     'upstream_intron': upstream_intron,
#                     'downstream_intron': downstream_intron,
#                     'upstream_exon': 'NA',
#                     'downstream_exon': 'NA'
#                 })
#                 continue
        
#         # Format CDS coordinates as comma-separated lists (deduplicated and sorted)
#         upstream_coords = list(set([f"{chrom}:{cds[0]}-{cds[1]}" for cds in upstream_cds_list]))
#         downstream_coords = list(set([f"{chrom}:{cds[0]}-{cds[1]}" for cds in downstream_cds_list]))
        
#         upstream_coords.sort()
#         downstream_coords.sort()
        
#         upstream_exon_str = ",".join(upstream_coords) if upstream_coords else "NA"
#         downstream_exon_str = ",".join(downstream_coords) if downstream_coords else "NA"
        
#         # Poison exon identifier
#         pe_id = f"{chrom}:{pe_actual_start}-{pe_actual_end}_{gene_name}"
        
#         results.append({
#             'chrom': chrom,
#             'start': pe_actual_start,
#             'end': pe_actual_end,
#             'name': pe_id,
#             'length': pe_actual_end - pe_actual_start + 1,
#             'strand': strand,
#             'upstream_intron': upstream_intron,
#             'downstream_intron': downstream_intron,
#             'upstream_exon': upstream_exon_str,
#             'downstream_exon': downstream_exon_str
#         })
        
#         # Add ALL flanking CDS to separate list
#         for cds in upstream_cds_list:
#             flanking_cds_features.append({
#                 'chrom': chrom,
#                 'start': cds[0],
#                 'end': cds[1],
#                 'name': f"{chrom}:{cds[0]}-{cds[1]}_{cds[2]}",
#                 'length': cds[1] - cds[0] + 1,
#                 'strand': strand,
#                 'position': 'upstream',
#                 'poisonExon': pe_id,
#                 'upstream_intron': upstream_intron,
#                 'downstream_intron': downstream_intron,
#                 'transcript_id': cds[3]
#             })
        
#         for cds in downstream_cds_list:
#             flanking_cds_features.append({
#                 'chrom': chrom,
#                 'start': cds[0],
#                 'end': cds[1],
#                 'name': f"{chrom}:{cds[0]}-{cds[1]}_{cds[2]}",
#                 'length': cds[1] - cds[0] + 1,
#                 'strand': strand,
#                 'position': 'downstream',
#                 'poisonExon': pe_id,
#                 'upstream_intron': upstream_intron,
#                 'downstream_intron': downstream_intron,
#                 'transcript_id': cds[3]
#             })
        
#         if upstream_cds_list or downstream_cds_list:
#             found_count += 1
    
#     print(f"\n  Found flanking CDS for: {found_count:,} poison exons")
#     print(f"  No CDS overlap for: {not_found_count:,} poison exons")
    
#     return pd.DataFrame(results), pd.DataFrame(flanking_cds_features)

# def main():
#     parser = argparse.ArgumentParser(
#         description='Identify flanking CDS features for poison exons using intron coordinates',
#         formatter_class=argparse.RawDescriptionHelpFormatter,
#         epilog="""
# Example:
#   python identify_flanking_cds.py \\
#       -g Reference.gtf \\
#       -p PE/_coding_noncoding_introns_all.tsv \\
#       -o output/poison_exons_with_flanks.bed \\
#       -f output/flanking_cds.bed

# Flanking CDS identification:
#   - Finds CDS features that END/START at intron boundaries
#   - Plus strand: upstream CDS ends at up_intron_start-1, downstream CDS starts at down_intron_end+1
#   - Minus strand: accounts for reversed biological direction
#   - Allows +/- 1bp tolerance for coordinate matching
#         """
#     )
    
#     parser.add_argument('-g', '--gtf', required=True,
#                        help='Reference GTF file')
#     parser.add_argument('-p', '--poison-exons', required=True,
#                        help='Poison exon TSV file all')
#     parser.add_argument('-o', '--output', required=True,
#                        help='Output: poison exons with flanking exons')
#     parser.add_argument('-f', '--flanking-output', required=True,
#                        help='Output: flanking CDS features only')
    
#     args = parser.parse_args()
    
#     print("="*70)
#     print("IDENTIFY FLANKING CDS FOR POISON EXONS")
#     print("="*70)
#     print(f"GTF: {args.gtf}")
#     print(f"Poison exons: {args.poison_exons}")
#     print(f"Output 1: {args.output}")
#     print(f"Output 2: {args.flanking_output}")
#     print("="*70 + "\n")
    
#     # Load GTF CDS features
#     cds_index = load_gtf_cds(args.gtf)
    
#     # Process poison exons
#     pe_with_flanks, flanking_cds = process_poison_exons(
#         args.poison_exons, cds_index
#     )
    
#     # Save outputs
#     print(f"\nSaving outputs...")
    
#     # Remove duplicates
#     print(f"  Removing duplicate rows...")
#     pe_before = len(pe_with_flanks)
#     flanking_before = len(flanking_cds)
    
#     pe_with_flanks = pe_with_flanks.drop_duplicates()
#     flanking_cds = flanking_cds.drop_duplicates()
    
#     pe_removed = pe_before - len(pe_with_flanks)
#     flanking_removed = flanking_before - len(flanking_cds)
    
#     if pe_removed > 0:
#         print(f"  Removed {pe_removed:,} duplicate poison exon rows")
#     else:
#         print(f"  No duplicate poison exon rows found")
    
#     if flanking_removed > 0:
#         print(f"  Removed {flanking_removed:,} duplicate flanking CDS rows")
#     else:
#         print(f"  No duplicate flanking CDS rows found")
    
#     # Output 1: Poison exons with flanking exons
#     pe_with_flanks.to_csv(args.output, sep='\t', index=False)
#     print(f"  Saved {len(pe_with_flanks):,} unique poison exons to: {args.output}")
    
#     # Output 2: Flanking CDS only
#     flanking_cds.to_csv(args.flanking_output, sep='\t', index=False)
#     print(f"  Saved {len(flanking_cds):,} unique flanking CDS features to: {args.flanking_output}")
    
#     print("\n" + "="*70)
#     print("COMPLETE")
#     print("="*70)

# if __name__ == "__main__":
#     main()

#!/usr/bin/env python3
"""
Identify flanking CDS features for poison exons using intron coordinates.

Finds CDS features that are immediately adjacent to the upstream and downstream
introns, which represent the flanking exons.

Author: Dylan Stermer
Date: December 2024
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

def load_gtf_cds(gtf_file):
    """
    Load CDS features from GTF, indexed by chromosome.
    
    Returns:
        cds_index: Dict {chrom: [(start, end, strand, gene_name, transcript_id)]}
    """
    
    print(f"Loading GTF: {gtf_file}")
    
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
    
    Returns 1-based coordinates matching GTF.
    """
    chrom, pos = coord_str.split(':')
    start, end = pos.split('-')
    return chrom, int(start), int(end)

def find_overlapping_cds(chrom, start, end, cds_index):
    """
    Find CDS features that overlap with given coordinates.
    
    Args:
        chrom: Chromosome
        start, end: 1-based coordinates
        cds_index: Dict of CDS features
    
    Returns:
        List of (start, end, strand, gene_name, transcript_id)
    """
    if chrom not in cds_index:
        return []
    
    overlapping = []
    for cds_start, cds_end, strand, gene_name, transcript_id in cds_index[chrom]:
        # Check for overlap
        if not (end < cds_start or start > cds_end):
            overlapping.append((cds_start, cds_end, strand, gene_name, transcript_id))
    
    return overlapping

def find_flanking_cds_by_introns(chrom, upstream_intron, downstream_intron, strand, cds_index):
    """
    Find CDS features that are adjacent to the introns (flanking the poison exon).
    
    Returns ALL matching CDS features, not just one.
    
    Returns:
        (upstream_cds_list, downstream_cds_list): Each as list of (start, end, gene_name, transcript_id) tuples
    """
    
    if chrom not in cds_index:
        return [], []
    
    # Parse intron coordinates
    _, up_intron_start, up_intron_end = parse_coordinate(upstream_intron)
    _, down_intron_start, down_intron_end = parse_coordinate(downstream_intron)
    
    upstream_cds_list = []
    downstream_cds_list = []
    
    if strand == '+':
        # Upstream CDS: should end at up_intron_start - 1
        for cds_start, cds_end, cds_strand, gene_name, transcript_id in cds_index[chrom]:
            if cds_strand != strand:
                continue
            
            if cds_end == up_intron_start - 1 or cds_end == up_intron_start or cds_end == up_intron_start - 2:
                upstream_cds_list.append((cds_start, cds_end, gene_name, transcript_id))
        
        # Downstream CDS: should start at down_intron_end + 1
        for cds_start, cds_end, cds_strand, gene_name, transcript_id in cds_index[chrom]:
            if cds_strand != strand:
                continue
            
            if cds_start == down_intron_end + 1 or cds_start == down_intron_end or cds_start == down_intron_end + 2:
                downstream_cds_list.append((cds_start, cds_end, gene_name, transcript_id))
    
    else:  # Minus strand
        # Upstream CDS (5', biologically): should start at down_intron_end + 1
        for cds_start, cds_end, cds_strand, gene_name, transcript_id in cds_index[chrom]:
            if cds_strand != strand:
                continue
            
            if cds_start == down_intron_end + 1 or cds_start == down_intron_end or cds_start == down_intron_end + 2:
                upstream_cds_list.append((cds_start, cds_end, gene_name, transcript_id))
        
        # Downstream CDS (3', biologically): should end at up_intron_start - 1
        for cds_start, cds_end, cds_strand, gene_name, transcript_id in cds_index[chrom]:
            if cds_strand != strand:
                continue
            
            if cds_end == up_intron_start - 1 or cds_end == up_intron_start or cds_end == up_intron_start - 2:
                downstream_cds_list.append((cds_start, cds_end, gene_name, transcript_id))
    
    return upstream_cds_list, downstream_cds_list

def process_poison_exons(pe_file, cds_index):
    """
    Process poison exons and find flanking CDS features.
    
    Returns:
        (pe_with_flanks, flanking_cds_list)
    """
    
    print(f"\nReading poison exons: {pe_file}")
    pe_df = pd.read_csv(pe_file, sep='\t')
    
    print(f"  Total poison exons: {len(pe_df)}")
    
    results = []
    flanking_cds_features = []
    
    found_count = 0
    not_found_count = 0
    
    for idx, row in pe_df.iterrows():
        if idx % 1000 == 0 and idx > 0:
            print(f"  Processed {idx:,}/{len(pe_df):,} poison exons...")
        
        coding_junction = row['coding_junction']
        upstream_intron = row['upstream_intron']
        downstream_intron = row['downstream_intron']
        
        # Parse intron coordinates to determine actual poison exon boundaries
        up_chrom, up_start, up_end = parse_coordinate(upstream_intron)
        down_chrom, down_start, down_end = parse_coordinate(downstream_intron)
        
        chrom = up_chrom
        
        # Actual poison exon coordinates (1-based, inclusive)
        pe_actual_start = up_end + 1
        pe_actual_end = down_start
        
        # First, try to find flanking CDS to get strand/gene information
        strand = None
        gene_name = None
        upstream_cds_list = []
        downstream_cds_list = []
        
        # Try plus strand first
        for test_strand in ['+', '-']:
            up_cds_list, down_cds_list = find_flanking_cds_by_introns(
                chrom, upstream_intron, downstream_intron, test_strand, cds_index
            )
            
            if up_cds_list or down_cds_list:
                # Found flanking CDS, get gene name from them
                strand = test_strand
                upstream_cds_list = up_cds_list
                downstream_cds_list = down_cds_list
                
                # Get gene name from whichever CDS we found
                if up_cds_list:
                    gene_name = up_cds_list[0][2]
                elif down_cds_list:
                    gene_name = down_cds_list[0][2]
                
                break
        
        # If we didn't find flanking CDS, try overlapping the poison exon itself
        if strand is None:
            overlapping_cds = find_overlapping_cds(chrom, pe_actual_start, pe_actual_end, cds_index)
            if overlapping_cds:
                cds_start, cds_end, strand, gene_name, transcript_id = overlapping_cds[0]
            else:
                not_found_count += 1
                pe_id = f"{chrom}:{pe_actual_start}-{pe_actual_end}_NA"
                results.append({
                    'chrom': chrom,
                    'start': pe_actual_start - 1,  # Convert to 0-based BED format
                    'end': pe_actual_end,
                    'name': pe_id,
                    'length': pe_actual_end - pe_actual_start + 1,
                    'strand': 'NA',
                    'upstream_intron': upstream_intron,
                    'downstream_intron': downstream_intron,
                    'upstream_exon': 'NA',
                    'downstream_exon': 'NA'
                })
                continue
        
        # Format CDS coordinates as comma-separated lists (deduplicated and sorted)
        # Convert to 0-based BED format (subtract 1 from start)
        upstream_coords = list(set([f"{chrom}:{cds[0]-1}-{cds[1]}" for cds in upstream_cds_list]))
        downstream_coords = list(set([f"{chrom}:{cds[0]-1}-{cds[1]}" for cds in downstream_cds_list]))
        
        upstream_coords.sort()
        downstream_coords.sort()
        
        upstream_exon_str = ",".join(upstream_coords) if upstream_coords else "NA"
        downstream_exon_str = ",".join(downstream_coords) if downstream_coords else "NA"
        
        # Poison exon identifier
        pe_id = f"{chrom}:{pe_actual_start}-{pe_actual_end}_{gene_name}"
        
        results.append({
            'chrom': chrom,
            'start': pe_actual_start - 1,  # Convert to 0-based BED format
            'end': pe_actual_end,
            'name': pe_id,
            'length': pe_actual_end - pe_actual_start + 1,
            'strand': strand,
            'upstream_intron': upstream_intron,
            'downstream_intron': downstream_intron,
            'upstream_exon': upstream_exon_str,
            'downstream_exon': downstream_exon_str
        })
        
        # Add ALL flanking CDS to separate list
        for cds in upstream_cds_list:
            flanking_cds_features.append({
                'chrom': chrom,
                'start': cds[0] - 1,  # Convert to 0-based BED format
                'end': cds[1],
                'name': f"{chrom}:{cds[0]-1}-{cds[1]}_{cds[2]}",
                'length': cds[1] - cds[0] + 1,
                'strand': strand,
                'position': 'upstream',
                'poisonExon': pe_id,
                'upstream_intron': upstream_intron,
                'downstream_intron': downstream_intron,
                'transcript_id': cds[3]
            })
        
        for cds in downstream_cds_list:
            flanking_cds_features.append({
                'chrom': chrom,
                'start': cds[0] - 1,  # Convert to 0-based BED format
                'end': cds[1],
                'name': f"{chrom}:{cds[0]-1}-{cds[1]}_{cds[2]}",
                'length': cds[1] - cds[0] + 1,
                'strand': strand,
                'position': 'downstream',
                'poisonExon': pe_id,
                'upstream_intron': upstream_intron,
                'downstream_intron': downstream_intron,
                'transcript_id': cds[3]
            })
        
        if upstream_cds_list or downstream_cds_list:
            found_count += 1
    
    print(f"\n  Found flanking CDS for: {found_count:,} poison exons")
    print(f"  No CDS overlap for: {not_found_count:,} poison exons")
    
    return pd.DataFrame(results), pd.DataFrame(flanking_cds_features)

def main():
    parser = argparse.ArgumentParser(
        description='Identify flanking CDS features for poison exons using intron coordinates',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python identify_flanking_cds.py \\
      -g Reference.gtf \\
      -p PE/_coding_noncoding_introns_all.tsv \\
      -o output/poison_exons_with_flanks.bed \\
      -f output/flanking_cds.bed

Flanking CDS identification:
  - Finds CDS features that END/START at intron boundaries
  - Plus strand: upstream CDS ends at up_intron_start-1, downstream CDS starts at down_intron_end+1
  - Minus strand: accounts for reversed biological direction
  - Allows +/- 1bp tolerance for coordinate matching
  
Output format:
  - All coordinates in 0-based BED format (start position is -1 from 1-based)
        """
    )
    
    parser.add_argument('-g', '--gtf', required=True,
                       help='Reference GTF file')
    parser.add_argument('-p', '--poison-exons', required=True,
                       help='Poison exon TSV file all')
    parser.add_argument('-o', '--output', required=True,
                       help='Output: poison exons with flanking exons')
    parser.add_argument('-f', '--flanking-output', required=True,
                       help='Output: flanking CDS features only')
    
    args = parser.parse_args()
    
    print("="*70)
    print("IDENTIFY FLANKING CDS FOR POISON EXONS")
    print("="*70)
    print(f"GTF: {args.gtf}")
    print(f"Poison exons: {args.poison_exons}")
    print(f"Output 1: {args.output}")
    print(f"Output 2: {args.flanking_output}")
    print("="*70 + "\n")
    
    # Load GTF CDS features
    cds_index = load_gtf_cds(args.gtf)
    
    # Process poison exons
    pe_with_flanks, flanking_cds = process_poison_exons(
        args.poison_exons, cds_index
    )
    
    # Save outputs
    print(f"\nSaving outputs...")
    
    # Remove duplicates
    print(f"  Removing duplicate rows...")
    pe_before = len(pe_with_flanks)
    flanking_before = len(flanking_cds)
    
    pe_with_flanks = pe_with_flanks.drop_duplicates()
    flanking_cds = flanking_cds.drop_duplicates()
    
    pe_removed = pe_before - len(pe_with_flanks)
    flanking_removed = flanking_before - len(flanking_cds)
    
    if pe_removed > 0:
        print(f"  Removed {pe_removed:,} duplicate poison exon rows")
    else:
        print(f"  No duplicate poison exon rows found")
    
    if flanking_removed > 0:
        print(f"  Removed {flanking_removed:,} duplicate flanking CDS rows")
    else:
        print(f"  No duplicate flanking CDS rows found")
    
    # Output 1: Poison exons with flanking exons
    pe_with_flanks.to_csv(args.output, sep='\t', index=False)
    print(f"  Saved {len(pe_with_flanks):,} unique poison exons to: {args.output}")
    
    # Output 2: Flanking CDS only
    flanking_cds.to_csv(args.flanking_output, sep='\t', index=False)
    print(f"  Saved {len(flanking_cds):,} unique flanking CDS features to: {args.flanking_output}")
    
    print("\n" + "="*70)
    print("COMPLETE")
    print("="*70)

if __name__ == "__main__":
    main()