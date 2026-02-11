# #!/usr/bin/env python3
# """
# Match LeafCutter cassette exons to GTF annotations.

# Prioritizes CDS matching, falls back to exon/UTR if no CDS match found.

# Author: Dylan Stermer
# Date: January 2026
# """

# import pandas as pd
# import argparse
# import sys
# from collections import defaultdict

# def parse_gtf_attributes(attr_string):
#     """Parse GTF attribute string into dictionary."""
#     attributes = {}
#     for item in attr_string.strip().split(';'):
#         item = item.strip()
#         if not item:
#             continue
#         try:
#             key, value = item.split(' ', 1)
#             value = value.strip('"')
#             attributes[key] = value
#         except ValueError:
#             continue
#     return attributes

# def load_gtf_features(gtf_file):
#     """
#     Load CDS, exon, and UTR features from GTF file.
    
#     Returns: 
#         cds_records: dict (chrom, start_0based, end_0based, strand) -> list of CDS records
#         exon_records: dict (chrom, start_0based, end_0based, strand) -> list of exon/UTR records
#         gene_transcripts: dict gene_name -> set of transcript_ids
#         gene_exons: dict gene_name -> dict of exon_key -> set of transcript_ids
#     """
#     print(f"Loading GTF file: {gtf_file}")
    
#     cds_records = defaultdict(list)
#     exon_records = defaultdict(list)
#     gene_transcripts = defaultdict(set)
#     gene_exons = defaultdict(lambda: defaultdict(set))
    
#     with open(gtf_file, 'r') as f:
#         for line_num, line in enumerate(f, 1):
#             if line.startswith('#'):
#                 continue
            
#             fields = line.strip().split('\t')
#             if len(fields) < 9:
#                 continue
            
#             feature_type = fields[2]
            
#             # Load CDS, exon, and UTR features
#             if feature_type not in ['CDS', 'exon', 'five_prime_utr', 'three_prime_utr', 'UTR']:
#                 continue
            
#             chrom = fields[0]
#             start_1based = int(fields[3])
#             end_1based = int(fields[4])
#             strand = fields[6]
            
#             # Convert to 0-based half-open (BED format)
#             start_0based = start_1based - 1
#             end_0based = end_1based
            
#             attributes = parse_gtf_attributes(fields[8])
            
#             gene_name = attributes.get('gene_name', 'NA')
#             gene_id = attributes.get('gene_id', 'NA')
#             transcript_id = attributes.get('transcript_id', 'NA')
            
#             exon_key = (chrom, start_0based, end_0based, strand)
            
#             record = {
#                 'chrom': chrom,
#                 'start': start_0based,
#                 'end': end_0based,
#                 'strand': strand,
#                 'feature_type': feature_type,
#                 'gene_id': gene_id,
#                 'gene_name': gene_name,
#                 'gene_type': attributes.get('gene_type', attributes.get('gene_biotype', 'NA')),
#                 'transcript_id': transcript_id,
#                 'transcript_type': attributes.get('transcript_type', attributes.get('transcript_biotype', 'NA')),
#                 'exon_number': attributes.get('exon_number', 'NA'),
#                 'exon_id': attributes.get('exon_id', 'NA')
#             }
            
#             # Separate CDS from exon/UTR
#             if feature_type == 'CDS':
#                 cds_records[exon_key].append(record)
#             else:
#                 exon_records[exon_key].append(record)
            
#             # Track for constitutive analysis
#             if gene_name != 'NA' and transcript_id != 'NA':
#                 gene_transcripts[gene_name].add(transcript_id)
#                 gene_exons[gene_name][exon_key].add(transcript_id)
            
#             if line_num % 100000 == 0:
#                 print(f"  Processed {line_num} GTF lines...")
    
#     print(f"Loaded {len(cds_records)} unique CDS coordinates")
#     print(f"Loaded {len(exon_records)} unique exon/UTR coordinates")
#     print(f"Tracked {len(gene_transcripts)} genes for constitutive analysis")
#     return cds_records, exon_records, gene_transcripts, gene_exons

# def extract_cds_from_gtf(gtf_file):
#     """Extract CDS regions per transcript for phase calculation."""
#     print("Extracting CDS coordinates for phase calculation...")
    
#     transcript_cds = defaultdict(list)
    
#     with open(gtf_file, 'r') as f:
#         for line in f:
#             if line.startswith('#'):
#                 continue
            
#             fields = line.strip().split('\t')
#             if len(fields) < 9:
#                 continue
            
#             if fields[2] != 'CDS':
#                 continue
            
#             chrom = fields[0]
#             start_1based = int(fields[3])
#             end_1based = int(fields[4])
#             strand = fields[6]
            
#             start_0based = start_1based - 1
#             end_0based = end_1based
            
#             attributes = parse_gtf_attributes(fields[8])
#             transcript_id = attributes.get('transcript_id', 'NA')
            
#             if transcript_id == 'NA':
#                 continue
            
#             transcript_cds[transcript_id].append({
#                 'chrom': chrom,
#                 'start': start_0based,
#                 'end': end_0based,
#                 'strand': strand
#             })
    
#     # Sort CDS exons by position
#     for tx_id in transcript_cds:
#         strand = transcript_cds[tx_id][0]['strand']
#         if strand == '+':
#             transcript_cds[tx_id].sort(key=lambda x: x['start'])
#         else:
#             transcript_cds[tx_id].sort(key=lambda x: x['start'], reverse=True)
    
#     print(f"Found CDS for {len(transcript_cds)} transcripts")
#     return transcript_cds

# def calculate_exon_phase(exon_coord, transcript_id, transcript_cds):
#     """Calculate phase and cumulative phase for a CDS exon."""
#     if transcript_id not in transcript_cds:
#         return None, None, None, None, None
    
#     cds_list = transcript_cds[transcript_id]
#     chrom, start, end, strand = exon_coord
    
#     # Find this exon in the CDS list
#     exon_idx = None
#     for i, cds in enumerate(cds_list):
#         if cds['chrom'] == chrom and cds['start'] == start and cds['end'] == end:
#             exon_idx = i
#             break
    
#     if exon_idx is None:
#         return None, None, None, None, None
    
#     # Calculate phase
#     if exon_idx == 0:
#         phase = 0
#     else:
#         cumulative_length = 0
#         for i in range(exon_idx):
#             upstream_cds = cds_list[i]
#             upstream_length = upstream_cds['end'] - upstream_cds['start']
#             cumulative_length += upstream_length
#         phase = cumulative_length % 3
    
#     # Calculate cumulative phase
#     current_length = end - start
#     cumulative_phase = (phase + current_length) % 3
    
#     # Determine exon position
#     total_cds_exons = len(cds_list)
#     exon_num = exon_idx + 1
    
#     if total_cds_exons == 1:
#         exon_position = 'single'
#     elif exon_num == 1:
#         exon_position = 'first'
#     elif exon_num == total_cds_exons:
#         exon_position = 'last'
#     else:
#         exon_position = 'internal'
    
#     return phase, cumulative_phase, exon_position, exon_num, total_cds_exons

# def parse_cassette_coord(coord_str):
#     """
#     Parse cassette coordinate string (BED format) to components.
#     Example: GL000194.1:62720-62949
#     Returns: (chrom, start_0based, end_0based)
#     """
#     chrom, pos = coord_str.split(':')
#     start, end = pos.split('-')
#     return chrom, int(start), int(end)

# def match_cassette_to_gtf(cassette, cds_records, exon_records, transcript_cds, 
#                          gene_transcripts, gene_exons):
#     """
#     Match a cassette exon to GTF. Prioritize CDS, fall back to exon/UTR.
    
#     Returns:
#         matched_dict or None
#     """
#     chrom, start_0based, end_0based = parse_cassette_coord(cassette['cassette_coord'])
#     strand = cassette['strand']
    
#     exon_key = (chrom, start_0based, end_0based, strand)
    
#     # Try CDS first
#     gtf_matches = None
#     matched_feature = None
    
#     if exon_key in cds_records:
#         gtf_matches = cds_records[exon_key]
#         matched_feature = 'CDS'
#     elif exon_key in exon_records:
#         gtf_matches = exon_records[exon_key]
#         matched_feature = 'exon/UTR'
#     else:
#         return None
    
#     # Aggregate information
#     transcript_ids = []
#     gtf_gene_names = set()
#     gene_ids = set()
#     gene_types = set()
#     tx_types = set()
#     feature_types = set()
    
#     phase_results = []
    
#     for match in gtf_matches:
#         transcript_ids.append(match['transcript_id'])
#         gtf_gene_names.add(match['gene_name'])
#         gene_ids.add(match['gene_id'])
#         gene_types.add(match['gene_type'])
#         tx_types.add(match['transcript_type'])
#         feature_types.add(match['feature_type'])
        
#         # Only calculate phase for CDS features
#         if match['feature_type'] == 'CDS':
#             phase_info = calculate_exon_phase(exon_key, match['transcript_id'], transcript_cds)
#             if phase_info[0] is not None:
#                 phase_results.append(phase_info)
    
#     # Determine constitutive status
#     constitutive_status = []
#     for gname in gtf_gene_names:
#         if gname in gene_transcripts and gname in gene_exons:
#             total_transcripts = len(gene_transcripts[gname])
#             transcripts_with_exon = len(gene_exons[gname][exon_key])
#             is_constitutive = (transcripts_with_exon == total_transcripts)
#             constitutive_status.append(is_constitutive)
#         else:
#             constitutive_status.append(False)
    
#     constitutive = all(constitutive_status) if constitutive_status else False
    
#     # Aggregate phase information
#     if phase_results:
#         phases = sorted(set([p[0] for p in phase_results]))
#         cumulative_phases = sorted(set([p[1] for p in phase_results]))
#         positions = sorted(set([p[2] for p in phase_results]))
#         exon_nums = sorted(set([p[3] for p in phase_results]))
#         total_exons = sorted(set([p[4] for p in phase_results]))
        
#         phase = ','.join(map(str, phases))
#         cumulative_phase = ','.join(map(str, cumulative_phases))
#         exon_position = ','.join(positions)
#         exon_num = ','.join(map(str, exon_nums))
#         total_exon_count = ','.join(map(str, total_exons))
#         all_phases = ','.join(map(str, phases))
#         all_cumulative_phases = ','.join(map(str, cumulative_phases))
#         phase_consistent = len(phases) == 1
#         cumulative_phase_consistent = len(cumulative_phases) == 1
        
#         exon_length = end_0based - start_0based
#         symmetry = exon_length % 3
#     else:
#         # Non-CDS feature or couldn't calculate phase
#         phase = 'NA'
#         cumulative_phase = 'NA'
#         exon_position = 'NA'
#         exon_num = 'NA'
#         total_exon_count = 'NA'
#         all_phases = 'NA'
#         all_cumulative_phases = 'NA'
#         phase_consistent = 'NA'
#         cumulative_phase_consistent = 'NA'
#         exon_length = end_0based - start_0based
#         symmetry = 'NA'
    
#     return {
#         'cassette_id': cassette['cassette_id'],
#         'input_gene_name': cassette['gene_name'],
#         'cassette_coord': cassette['cassette_coord'],
#         'chrom': chrom,
#         'start': start_0based,
#         'end': end_0based,
#         'strand': strand,
#         'cassette_length': cassette['cassette_length'],
#         'cluster': cassette['cluster'],
#         'skip_detected': cassette['skip_detected'],
#         'length': exon_length,
#         'symmetry': symmetry,
#         'phase': phase,
#         'all_phases': all_phases,
#         'phase_consistent': phase_consistent,
#         'cumulative_phase': cumulative_phase,
#         'all_cumulative_phases': all_cumulative_phases,
#         'cumulative_phase_consistent': cumulative_phase_consistent,
#         'exon_position': exon_position,
#         'exon_num': exon_num,
#         'total_exons': total_exon_count,
#         'constitutive': constitutive,
#         'feature_types': ','.join(sorted(feature_types)),
#         'tx_names': ','.join(transcript_ids),
#         'n_transcripts': len(transcript_ids),
#         'gtf_gene_name': ','.join(sorted(gtf_gene_names)),
#         'gene_id': ','.join(sorted(gene_ids)),
#         'gene_type': ','.join(sorted(gene_types)),
#         'tx_type': ','.join(sorted(tx_types))
#     }

# def process_cassette_exons(cassette_file, cds_records, exon_records, transcript_cds,
#                           gene_transcripts, gene_exons):
#     """
#     Process cassette exons and match to GTF.
    
#     Returns:
#         matched: list of matched cassettes
#         unmatched: list of unmatched cassettes
#     """
#     print(f"\nReading cassette exons: {cassette_file}")
#     df = pd.read_csv(cassette_file, sep='\t')
    
#     print(f"  Total cassette exons: {len(df)}")
    
#     matched = []
#     unmatched = []
    
#     for idx, row in df.iterrows():
#         if (idx + 1) % 1000 == 0:
#             print(f"  Processed {idx+1}/{len(df)} cassettes...")
        
#         cassette = {
#             'cassette_id': row['cassette_id'],
#             'gene_name': row['gene_name'],
#             'cassette_coord': row['cassette_coord'],
#             'strand': row['strand'],
#             'cassette_length': row['cassette_length'],
#             'cluster': row['cluster'],
#             'skip_detected': row['skip_detected']
#         }
        
#         result = match_cassette_to_gtf(
#             cassette, cds_records, exon_records, transcript_cds,
#             gene_transcripts, gene_exons
#         )
        
#         if result:
#             matched.append(result)
#         else:
#             unmatched.append(cassette)
    
#     print(f"\nMatching complete:")
#     print(f"  Matched: {len(matched)}")
#     print(f"  Unmatched: {len(unmatched)}")
    
#     return matched, unmatched

# def write_matched(matched_list, output_file):
#     """Write matched cassettes to file."""
#     df = pd.DataFrame(matched_list)
    
#     if len(df) == 0:
#         print(f"WARNING: No matched cassettes")
#         with open(output_file, 'w') as f:
#             f.write("# No matched cassettes\n")
#         return
    
#     # Sort
#     def chrom_sort_key(chrom):
#         if chrom.startswith('chr'):
#             chrom = chrom[3:]
#         if chrom.isdigit():
#             return (0, int(chrom))
#         elif chrom == 'X':
#             return (1, 0)
#         elif chrom == 'Y':
#             return (1, 1)
#         elif chrom == 'M' or chrom == 'MT':
#             return (1, 2)
#         else:
#             return (2, chrom)
    
#     df['sort_key'] = df['chrom'].apply(chrom_sort_key)
#     df = df.sort_values(['sort_key', 'start'])
#     df = df.drop('sort_key', axis=1)
    
#     with open(output_file, 'w') as f:
#         f.write("# LeafCutter cassette exons matched to GTF\n")
#         f.write("# Columns: cassette_id, input_gene_name, cassette_coord, chrom, start, end, strand, ")
#         f.write("cassette_length, cluster, skip_detected, length, symmetry, phase, all_phases, ")
#         f.write("phase_consistent, cumulative_phase, all_cumulative_phases, cumulative_phase_consistent, ")
#         f.write("exon_position, exon_num, total_exons, constitutive, feature_types, tx_names, ")
#         f.write("n_transcripts, gtf_gene_name, gene_id, gene_type, tx_type\n")
#         f.write("# Coordinates: 0-based half-open (BED format)\n")
#         f.write("# Prioritizes CDS matching, falls back to exon/UTR\n")
        
#         for _, row in df.iterrows():
#             f.write(f"{row['cassette_id']}\t{row['input_gene_name']}\t{row['cassette_coord']}\t")
#             f.write(f"{row['chrom']}\t{row['start']}\t{row['end']}\t{row['strand']}\t")
#             f.write(f"{row['cassette_length']}\t{row['cluster']}\t{row['skip_detected']}\t")
#             f.write(f"{row['length']}\t{row['symmetry']}\t{row['phase']}\t{row['all_phases']}\t")
#             f.write(f"{row['phase_consistent']}\t{row['cumulative_phase']}\t")
#             f.write(f"{row['all_cumulative_phases']}\t{row['cumulative_phase_consistent']}\t")
#             f.write(f"{row['exon_position']}\t{row['exon_num']}\t{row['total_exons']}\t")
#             f.write(f"{row['constitutive']}\t{row['feature_types']}\t{row['tx_names']}\t")
#             f.write(f"{row['n_transcripts']}\t{row['gtf_gene_name']}\t{row['gene_id']}\t")
#             f.write(f"{row['gene_type']}\t{row['tx_type']}\n")
    
#     print(f"\nWrote {len(df)} matched cassettes to: {output_file}")
    
#     # Summary stats
#     print(f"\nMATCHED CASSETTE SUMMARY:")
#     print(f"  Total matched: {len(df)}")
    
#     # Feature type breakdown
#     print(f"\n  Feature type breakdown:")
#     feature_counts = defaultdict(int)
#     for ft in df['feature_types']:
#         for feature in ft.split(','):
#             feature_counts[feature] += 1
#     for feature, count in sorted(feature_counts.items()):
#         print(f"    {feature}: {count}")
    
#     # Constitutive breakdown
#     constitutive_count = sum(df['constitutive'])
#     alternative_count = sum(~df['constitutive'])
#     print(f"\n  Constitutive vs Alternative:")
#     print(f"    Constitutive: {constitutive_count} ({constitutive_count/len(df)*100:.1f}%)")
#     print(f"    Alternative: {alternative_count} ({alternative_count/len(df)*100:.1f}%)")

# def write_unmatched(unmatched_list, output_file):
#     """Write unmatched cassettes to file."""
#     df = pd.DataFrame(unmatched_list)
    
#     if len(df) == 0:
#         print(f"INFO: All cassettes matched!")
#         with open(output_file, 'w') as f:
#             f.write("# All cassettes matched\n")
#         return
    
#     # Sort
#     def chrom_sort_key_from_coord(coord):
#         chrom = coord.split(':')[0]
#         if chrom.startswith('chr'):
#             chrom = chrom[3:]
#         if chrom.isdigit():
#             return (0, int(chrom))
#         elif chrom == 'X':
#             return (1, 0)
#         elif chrom == 'Y':
#             return (1, 1)
#         elif chrom == 'M' or chrom == 'MT':
#             return (1, 2)
#         else:
#             return (2, chrom)
    
#     df['sort_key'] = df['cassette_coord'].apply(chrom_sort_key_from_coord)
#     df = df.sort_values('sort_key')
#     df = df.drop('sort_key', axis=1)
    
#     with open(output_file, 'w') as f:
#         f.write("# LeafCutter cassette exons with no GTF match\n")
#         f.write("# Columns: cassette_id, cassette_coord, strand, cassette_length, skip_detected, cluster\n")
        
#         for _, row in df.iterrows():
#             f.write(f"{row['cassette_id']}\t{row['cassette_coord']}\t{row['strand']}\t")
#             f.write(f"{row['cassette_length']}\t{row['skip_detected']}\t{row['cluster']}\n")
    
#     print(f"\nWrote {len(df)} unmatched cassettes to: {output_file}")

# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(
#         description='Match LeafCutter cassette exons to GTF annotations',
#         formatter_class=argparse.RawDescriptionHelpFormatter,
#         epilog="""
# Matches LeafCutter cassette exons to GTF features.
# Prioritizes CDS matching, falls back to exon/UTR if no CDS match.

# Input:
#   - Cassette exon TSV with columns: cassette_id, gene_name, cassette_coord, 
#     strand, cassette_length, cluster, skip_detected
#   - Cassette coordinates in BED format (0-based, half-open)

# Output:
#   1. *_cassette_gtf_match.bed: Matched cassettes with full annotation
#   2. *_cassette_gtf_no_match.bed: Unmatched cassettes (minimal format)

# Example:
#   python match_cassette_to_gtf.py \\
#       -c cassetteExons.tsv \\
#       -g reference.gtf \\
#       -o output/leafcutter
#         """
#     )
    
#     parser.add_argument('-c', '--cassettes', required=True,
#                        help='Input cassette exon TSV file')
#     parser.add_argument('-g', '--gtf', required=True,
#                        help='Reference GTF file')
#     parser.add_argument('-o', '--output', required=True,
#                        help='Output prefix for result files')
    
#     args = parser.parse_args()
    
#     print("="*70)
#     print("MATCH LEAFCUTTER CASSETTES TO GTF")
#     print("="*70)
#     print(f"Cassettes: {args.cassettes}")
#     print(f"GTF: {args.gtf}")
#     print(f"Output prefix: {args.output}")
#     print("="*70 + "\n")
    
#     # Load GTF
#     cds_records, exon_records, gene_transcripts, gene_exons = load_gtf_features(args.gtf)
#     transcript_cds = extract_cds_from_gtf(args.gtf)
    
#     # Process cassettes
#     matched, unmatched = process_cassette_exons(
#         args.cassettes, cds_records, exon_records, transcript_cds,
#         gene_transcripts, gene_exons
#     )
    
#     # Write outputs
#     print("\n" + "="*70)
#     print("WRITING OUTPUT FILES")
#     print("="*70)
    
#     match_file = f"{args.output}_cassette_gtf_match.bed"
#     write_matched(matched, match_file)
    
#     print("\n" + "-"*70)
    
#     nomatch_file = f"{args.output}_cassette_gtf_no_match.bed"
#     write_unmatched(unmatched, nomatch_file)
    
#     print("\n" + "="*70)
#     print("MATCHING COMPLETE")
#     print("="*70)
#     print(f"\nOutput files:")
#     print(f"  1. {match_file}")
#     print(f"  2. {nomatch_file}")

## Monday 

#!/usr/bin/env python3
"""
Match LeafCutter cassette exons to GTF annotations.

Prioritizes CDS matching, falls back to exon/UTR if no CDS match found.

STANDARDIZED OUTPUT COLUMNS:
chrom, start, end, name, strand, length, symmetry, phase, all_phases, phase_consistent,
cumulative_phase, all_cumulative_phases, cumulative_phase_consistent, exon_position,
exon_num, total_exons, constitutive, feature_types, tx_names, n_transcripts, gene_name,
gene_id, gene_type, tx_type

Author: Dylan Stermer
Date: January 2026
"""

import pandas as pd
import argparse
import sys
from collections import defaultdict

def parse_gtf_attributes(attr_string):
    """Parse GTF attribute string into dictionary."""
    attributes = {}
    for item in attr_string.strip().split(';'):
        item = item.strip()
        if not item:
            continue
        try:
            key, value = item.split(' ', 1)
            value = value.strip('"')
            attributes[key] = value
        except ValueError:
            continue
    return attributes

def load_gtf_features(gtf_file):
    """
    Load CDS, exon, and UTR features from GTF file.
    
    Returns: 
        cds_records: dict (chrom, start_0based, end_0based, strand) -> list of CDS records
        exon_records: dict (chrom, start_0based, end_0based, strand) -> list of exon/UTR records
        gene_transcripts: dict gene_name -> set of transcript_ids
        gene_exons: dict gene_name -> dict of exon_key -> set of transcript_ids
    """
    print(f"Loading GTF file: {gtf_file}")
    
    cds_records = defaultdict(list)
    exon_records = defaultdict(list)
    gene_transcripts = defaultdict(set)
    gene_exons = defaultdict(lambda: defaultdict(set))
    
    with open(gtf_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            feature_type = fields[2]
            
            # Load CDS, exon, and UTR features
            if feature_type not in ['CDS', 'exon', 'five_prime_utr', 'three_prime_utr', 'UTR']:
                continue
            
            chrom = fields[0]
            start_1based = int(fields[3])
            end_1based = int(fields[4])
            strand = fields[6]
            
            # Convert to 0-based half-open (BED format)
            start_0based = start_1based - 1
            end_0based = end_1based
            
            attributes = parse_gtf_attributes(fields[8])
            
            gene_name = attributes.get('gene_name', 'NA')
            gene_id = attributes.get('gene_id', 'NA')
            transcript_id = attributes.get('transcript_id', 'NA')
            
            exon_key = (chrom, start_0based, end_0based, strand)
            
            record = {
                'chrom': chrom,
                'start': start_0based,
                'end': end_0based,
                'strand': strand,
                'feature_type': feature_type,
                'gene_id': gene_id,
                'gene_name': gene_name,
                'gene_type': attributes.get('gene_type', attributes.get('gene_biotype', 'NA')),
                'transcript_id': transcript_id,
                'transcript_type': attributes.get('transcript_type', attributes.get('transcript_biotype', 'NA')),
                'exon_number': attributes.get('exon_number', 'NA'),
                'exon_id': attributes.get('exon_id', 'NA')
            }
            
            # Separate CDS from exon/UTR
            if feature_type == 'CDS':
                cds_records[exon_key].append(record)
            else:
                exon_records[exon_key].append(record)
            
            # Track for constitutive analysis
            if gene_name != 'NA' and transcript_id != 'NA':
                gene_transcripts[gene_name].add(transcript_id)
                gene_exons[gene_name][exon_key].add(transcript_id)
            
            if line_num % 100000 == 0:
                print(f"  Processed {line_num} GTF lines...")
    
    print(f"Loaded {len(cds_records)} unique CDS coordinates")
    print(f"Loaded {len(exon_records)} unique exon/UTR coordinates")
    print(f"Tracked {len(gene_transcripts)} genes for constitutive analysis")
    return cds_records, exon_records, gene_transcripts, gene_exons

def extract_cds_from_gtf(gtf_file):
    """Extract CDS regions per transcript for phase calculation."""
    print("Extracting CDS coordinates for phase calculation...")
    
    transcript_cds = defaultdict(list)
    
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            if fields[2] != 'CDS':
                continue
            
            chrom = fields[0]
            start_1based = int(fields[3])
            end_1based = int(fields[4])
            strand = fields[6]
            
            start_0based = start_1based - 1
            end_0based = end_1based
            
            attributes = parse_gtf_attributes(fields[8])
            transcript_id = attributes.get('transcript_id', 'NA')
            
            if transcript_id == 'NA':
                continue
            
            transcript_cds[transcript_id].append({
                'chrom': chrom,
                'start': start_0based,
                'end': end_0based,
                'strand': strand
            })
    
    # Sort CDS exons by position
    for tx_id in transcript_cds:
        strand = transcript_cds[tx_id][0]['strand']
        if strand == '+':
            transcript_cds[tx_id].sort(key=lambda x: x['start'])
        else:
            transcript_cds[tx_id].sort(key=lambda x: x['start'], reverse=True)
    
    print(f"Found CDS for {len(transcript_cds)} transcripts")
    return transcript_cds

def calculate_exon_phase(exon_coord, transcript_id, transcript_cds):
    """Calculate phase and cumulative phase for a CDS exon."""
    if transcript_id not in transcript_cds:
        return None, None, None, None, None
    
    cds_list = transcript_cds[transcript_id]
    chrom, start, end, strand = exon_coord
    
    # Find this exon in the CDS list
    exon_idx = None
    for i, cds in enumerate(cds_list):
        if cds['chrom'] == chrom and cds['start'] == start and cds['end'] == end:
            exon_idx = i
            break
    
    if exon_idx is None:
        return None, None, None, None, None
    
    # Calculate phase
    if exon_idx == 0:
        phase = 0
    else:
        cumulative_length = 0
        for i in range(exon_idx):
            upstream_cds = cds_list[i]
            upstream_length = upstream_cds['end'] - upstream_cds['start']
            cumulative_length += upstream_length
        phase = cumulative_length % 3
    
    # Calculate cumulative phase
    current_length = end - start
    cumulative_phase = (phase + current_length) % 3
    
    # Determine exon position
    total_cds_exons = len(cds_list)
    exon_num = exon_idx + 1
    
    if total_cds_exons == 1:
        exon_position = 'single'
    elif exon_num == 1:
        exon_position = 'first'
    elif exon_num == total_cds_exons:
        exon_position = 'last'
    else:
        exon_position = 'internal'
    
    return phase, cumulative_phase, exon_position, exon_num, total_cds_exons

def parse_cassette_coord(coord_str):
    """
    Parse cassette coordinate string (BED format) to components.
    Example: GL000194.1:62720-62949
    Returns: (chrom, start_0based, end_0based)
    """
    chrom, pos = coord_str.split(':')
    start, end = pos.split('-')
    return chrom, int(start), int(end)

def match_cassette_to_gtf(cassette, cds_records, exon_records, transcript_cds, 
                         gene_transcripts, gene_exons):
    """
    Match a cassette exon to GTF. Prioritize CDS, fall back to exon/UTR.
    
    Returns:
        matched_dict or None
    """
    chrom, start_0based, end_0based = parse_cassette_coord(cassette['cassette_coord'])
    strand = cassette['strand']
    
    exon_key = (chrom, start_0based, end_0based, strand)
    
    # Try CDS first
    gtf_matches = None
    matched_feature = None
    
    if exon_key in cds_records:
        gtf_matches = cds_records[exon_key]
        matched_feature = 'CDS'
    elif exon_key in exon_records:
        gtf_matches = exon_records[exon_key]
        matched_feature = 'exon/UTR'
    else:
        return None
    
    # Aggregate information
    transcript_ids = []
    gtf_gene_names = set()
    gene_ids = set()
    gene_types = set()
    tx_types = set()
    feature_types = set()
    
    phase_results = []
    
    for match in gtf_matches:
        transcript_ids.append(match['transcript_id'])
        gtf_gene_names.add(match['gene_name'])
        gene_ids.add(match['gene_id'])
        gene_types.add(match['gene_type'])
        tx_types.add(match['transcript_type'])
        feature_types.add(match['feature_type'])
        
        # Only calculate phase for CDS features
        if match['feature_type'] == 'CDS':
            phase_info = calculate_exon_phase(exon_key, match['transcript_id'], transcript_cds)
            if phase_info[0] is not None:
                phase_results.append(phase_info)
    
    # Determine constitutive status
    constitutive_status = []
    for gname in gtf_gene_names:
        if gname in gene_transcripts and gname in gene_exons:
            total_transcripts = len(gene_transcripts[gname])
            transcripts_with_exon = len(gene_exons[gname][exon_key])
            is_constitutive = (transcripts_with_exon == total_transcripts)
            constitutive_status.append(is_constitutive)
        else:
            constitutive_status.append(False)
    
    constitutive = all(constitutive_status) if constitutive_status else False
    
    exon_length = end_0based - start_0based
    
    # Aggregate phase information
    if phase_results:
        phases = sorted(set([p[0] for p in phase_results]))
        cumulative_phases = sorted(set([p[1] for p in phase_results]))
        positions = sorted(set([p[2] for p in phase_results]))
        exon_nums = sorted(set([p[3] for p in phase_results]))
        total_exons = sorted(set([p[4] for p in phase_results]))
        
        # STANDARDIZED: Always use comma-separated strings for phase
        phase = ','.join(map(str, phases))
        cumulative_phase = ','.join(map(str, cumulative_phases))
        exon_position = ','.join(positions)
        exon_num = ','.join(map(str, exon_nums))
        total_exon_count = ','.join(map(str, total_exons))
        all_phases = ','.join(map(str, phases))
        all_cumulative_phases = ','.join(map(str, cumulative_phases))
        phase_consistent = len(phases) == 1
        cumulative_phase_consistent = len(cumulative_phases) == 1
        
        symmetry = exon_length % 3
    else:
        # Non-CDS feature or couldn't calculate phase
        phase = 'NA'
        cumulative_phase = 'NA'
        exon_position = 'NA'
        exon_num = 'NA'
        total_exon_count = 'NA'
        all_phases = 'NA'
        all_cumulative_phases = 'NA'
        phase_consistent = 'NA'
        cumulative_phase_consistent = 'NA'
        symmetry = 'NA'
    
    # Create standardized name
    name = f"{cassette['cassette_id']}_{cassette['gene_name']}"
    
    # STANDARDIZED OUTPUT - matching the expected column order
    return {
        'chrom': chrom,
        'start': start_0based,
        'end': end_0based,
        'name': name,
        'strand': strand,
        'length': exon_length,
        'symmetry': symmetry,
        'phase': phase,
        'all_phases': all_phases,
        'phase_consistent': phase_consistent,
        'cumulative_phase': cumulative_phase,
        'all_cumulative_phases': all_cumulative_phases,
        'cumulative_phase_consistent': cumulative_phase_consistent,
        'exon_position': exon_position,
        'exon_num': exon_num,
        'total_exons': total_exon_count,
        'constitutive': constitutive,
        'feature_types': ','.join(sorted(feature_types)),
        'tx_names': ','.join(transcript_ids),
        'n_transcripts': len(transcript_ids),
        'gene_name': ','.join(sorted(gtf_gene_names)),
        'gene_id': ','.join(sorted(gene_ids)),
        'gene_type': ','.join(sorted(gene_types)),
        'tx_type': ','.join(sorted(tx_types)),
        # Additional cassette-specific columns (appended after standard)
        'cassette_id': cassette['cassette_id'],
        'input_gene_name': cassette['gene_name'],
        'cassette_coord': cassette['cassette_coord'],
        'cassette_length': cassette['cassette_length'],
        'cluster': cassette['cluster'],
        'skip_detected': cassette['skip_detected']
    }

def process_cassette_exons(cassette_file, cds_records, exon_records, transcript_cds,
                          gene_transcripts, gene_exons):
    """
    Process cassette exons and match to GTF.
    
    Returns:
        matched: list of matched cassettes
        unmatched: list of unmatched cassettes
    """
    print(f"\nReading cassette exons: {cassette_file}")
    df = pd.read_csv(cassette_file, sep='\t')
    
    print(f"  Total cassette exons: {len(df)}")
    
    matched = []
    unmatched = []
    
    for idx, row in df.iterrows():
        if (idx + 1) % 1000 == 0:
            print(f"  Processed {idx+1}/{len(df)} cassettes...")
        
        cassette = {
            'cassette_id': row['cassette_id'],
            'gene_name': row['gene_name'],
            'cassette_coord': row['cassette_coord'],
            'strand': row['strand'],
            'cassette_length': row['cassette_length'],
            'cluster': row['cluster'],
            'skip_detected': row['skip_detected']
        }
        
        result = match_cassette_to_gtf(
            cassette, cds_records, exon_records, transcript_cds,
            gene_transcripts, gene_exons
        )
        
        if result:
            matched.append(result)
        else:
            unmatched.append(cassette)
    
    print(f"\nMatching complete:")
    print(f"  Matched: {len(matched)}")
    print(f"  Unmatched: {len(unmatched)}")
    
    return matched, unmatched

def write_matched(matched_list, output_file):
    """Write matched cassettes to file with STANDARDIZED columns."""
    df = pd.DataFrame(matched_list)
    
    if len(df) == 0:
        print(f"WARNING: No matched cassettes")
        with open(output_file, 'w') as f:
            f.write("# No matched cassettes\n")
        return
    
    # Sort
    def chrom_sort_key(chrom):
        if chrom.startswith('chr'):
            chrom = chrom[3:]
        if chrom.isdigit():
            return (0, int(chrom))
        elif chrom == 'X':
            return (1, 0)
        elif chrom == 'Y':
            return (1, 1)
        elif chrom == 'M' or chrom == 'MT':
            return (1, 2)
        else:
            return (2, chrom)
    
    df['sort_key'] = df['chrom'].apply(chrom_sort_key)
    df = df.sort_values(['sort_key', 'start'])
    df = df.drop('sort_key', axis=1)
    
    # STANDARDIZED column order
    standard_cols = [
        'chrom', 'start', 'end', 'name', 'strand', 'length', 'symmetry',
        'phase', 'all_phases', 'phase_consistent',
        'cumulative_phase', 'all_cumulative_phases', 'cumulative_phase_consistent',
        'exon_position', 'exon_num', 'total_exons', 'constitutive',
        'feature_types', 'tx_names', 'n_transcripts',
        'gene_name', 'gene_id', 'gene_type', 'tx_type',
        # Cassette-specific columns after standard
        'cassette_id', 'input_gene_name', 'cassette_coord', 'cassette_length', 
        'cluster', 'skip_detected'
    ]
    
    # Reorder columns
    df = df[standard_cols]
    
    with open(output_file, 'w') as f:
        f.write("# LeafCutter cassette exons matched to GTF\n")
        f.write("# STANDARDIZED COLUMNS: " + ", ".join(standard_cols) + "\n")
        f.write("# Coordinates: 0-based half-open (BED format)\n")
        f.write("# Phase: Reading frame at START of exon (comma-separated if multiple isoforms)\n")
        f.write("# Prioritizes CDS matching, falls back to exon/UTR\n")
        
        # Write header
        f.write('\t'.join(standard_cols) + '\n')
        
        for _, row in df.iterrows():
            values = [str(row[col]) for col in standard_cols]
            f.write('\t'.join(values) + '\n')
    
    print(f"\nWrote {len(df)} matched cassettes to: {output_file}")
    
    # Summary stats
    print(f"\nMATCHED CASSETTE SUMMARY:")
    print(f"  Total matched: {len(df)}")
    
    # Feature type breakdown
    print(f"\n  Feature type breakdown:")
    feature_counts = defaultdict(int)
    for ft in df['feature_types']:
        for feature in ft.split(','):
            feature_counts[feature] += 1
    for feature, count in sorted(feature_counts.items()):
        print(f"    {feature}: {count}")
    
    # Constitutive breakdown
    constitutive_count = sum(df['constitutive'])
    alternative_count = sum(~df['constitutive'])
    print(f"\n  Constitutive vs Alternative:")
    print(f"    Constitutive: {constitutive_count} ({constitutive_count/len(df)*100:.1f}%)")
    print(f"    Alternative: {alternative_count} ({alternative_count/len(df)*100:.1f}%)")

def write_unmatched(unmatched_list, output_file):
    """Write unmatched cassettes to file."""
    df = pd.DataFrame(unmatched_list)
    
    if len(df) == 0:
        print(f"INFO: All cassettes matched!")
        with open(output_file, 'w') as f:
            f.write("# All cassettes matched\n")
        return
    
    # Sort
    def chrom_sort_key_from_coord(coord):
        chrom = coord.split(':')[0]
        if chrom.startswith('chr'):
            chrom = chrom[3:]
        if chrom.isdigit():
            return (0, int(chrom))
        elif chrom == 'X':
            return (1, 0)
        elif chrom == 'Y':
            return (1, 1)
        elif chrom == 'M' or chrom == 'MT':
            return (1, 2)
        else:
            return (2, chrom)
    
    df['sort_key'] = df['cassette_coord'].apply(chrom_sort_key_from_coord)
    df = df.sort_values('sort_key')
    df = df.drop('sort_key', axis=1)
    
    with open(output_file, 'w') as f:
        f.write("# LeafCutter cassette exons with no GTF match\n")
        f.write("# Columns: cassette_id, cassette_coord, strand, cassette_length, skip_detected, cluster\n")
        
        for _, row in df.iterrows():
            f.write(f"{row['cassette_id']}\t{row['cassette_coord']}\t{row['strand']}\t")
            f.write(f"{row['cassette_length']}\t{row['skip_detected']}\t{row['cluster']}\n")
    
    print(f"\nWrote {len(df)} unmatched cassettes to: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Match LeafCutter cassette exons to GTF annotations',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Matches LeafCutter cassette exons to GTF features.
Prioritizes CDS matching, falls back to exon/UTR if no CDS match.

STANDARDIZED OUTPUT COLUMNS:
chrom, start, end, name, strand, length, symmetry, phase, all_phases, phase_consistent,
cumulative_phase, all_cumulative_phases, cumulative_phase_consistent, exon_position,
exon_num, total_exons, constitutive, feature_types, tx_names, n_transcripts, gene_name,
gene_id, gene_type, tx_type, [cassette-specific columns]

Input:
  - Cassette exon TSV with columns: cassette_id, gene_name, cassette_coord, 
    strand, cassette_length, cluster, skip_detected
  - Cassette coordinates in BED format (0-based, half-open)

Output:
  1. *_cassette_gtf_match.bed: Matched cassettes with full annotation
  2. *_cassette_gtf_no_match.bed: Unmatched cassettes (minimal format)

Example:
  python match_cassette_to_gtf.py \\
      -c cassetteExons.tsv \\
      -g reference.gtf \\
      -o output/leafcutter
        """
    )
    
    parser.add_argument('-c', '--cassettes', required=True,
                       help='Input cassette exon TSV file')
    parser.add_argument('-g', '--gtf', required=True,
                       help='Reference GTF file')
    parser.add_argument('-o', '--output', required=True,
                       help='Output prefix for result files')
    
    args = parser.parse_args()
    
    print("="*70)
    print("MATCH LEAFCUTTER CASSETTES TO GTF")
    print("="*70)
    print(f"Cassettes: {args.cassettes}")
    print(f"GTF: {args.gtf}")
    print(f"Output prefix: {args.output}")
    print("="*70 + "\n")
    
    # Load GTF
    cds_records, exon_records, gene_transcripts, gene_exons = load_gtf_features(args.gtf)
    transcript_cds = extract_cds_from_gtf(args.gtf)
    
    # Process cassettes
    matched, unmatched = process_cassette_exons(
        args.cassettes, cds_records, exon_records, transcript_cds,
        gene_transcripts, gene_exons
    )
    
    # Write outputs
    print("\n" + "="*70)
    print("WRITING OUTPUT FILES")
    print("="*70)
    
    match_file = f"{args.output}_cassette_gtf_match.bed"
    write_matched(matched, match_file)
    
    print("\n" + "-"*70)
    
    nomatch_file = f"{args.output}_cassette_gtf_no_match.bed"
    write_unmatched(unmatched, nomatch_file)
    
    print("\n" + "="*70)
    print("MATCHING COMPLETE")
    print("="*70)
    print(f"\nOutput files:")
    print(f"  1. {match_file}")
    print(f"  2. {nomatch_file}")