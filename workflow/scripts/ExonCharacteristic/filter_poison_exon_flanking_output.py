# #!/usr/bin/env python3
# """
# Filter poison exon flanking output and match flanking exons to GTF annotations.

# Creates three output files:
# 1. Flanking exons matched to GTF with full annotation
# 2. Flanking exons not matched to GTF
# 3. Poison exons with phase information from flanking exons

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

# def load_gtf_exons(gtf_file):
#     """
#     Load CDS coordinates from GTF file for flanking exon matching.
    
#     Returns: 
#         exon_records: dict (chrom, start_0based, end_0based, strand) -> list of CDS records
#         gene_transcripts: dict gene_name -> set of transcript_ids
#         gene_exons: dict gene_name -> dict of exon_key -> set of transcript_ids
    
#     Note: Only loads CDS features. Converts GTF 1-based closed coordinates to 0-based half-open.
#     """
#     print(f"Loading GTF file (CDS only): {gtf_file}")
    
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
            
#             # ONLY load CDS features
#             if feature_type != 'CDS':
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
            
#             exon_records[exon_key].append({
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
#             })
            
#             if gene_name != 'NA' and transcript_id != 'NA':
#                 gene_transcripts[gene_name].add(transcript_id)
#                 gene_exons[gene_name][exon_key].add(transcript_id)
            
#             if line_num % 100000 == 0:
#                 print(f"  Processed {line_num} GTF lines...")
    
#     print(f"Loaded {len(exon_records)} unique CDS coordinates")
#     print(f"Tracked {len(gene_transcripts)} genes for constitutive analysis")
#     return exon_records, gene_transcripts, gene_exons

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

# def parse_coordinate_string(coord_str):
#     """Parse coordinate string chr:start-end into (chr, start, end)."""
#     chrom, pos = coord_str.split(':')
#     start, end = pos.split('-')
#     return chrom, int(start), int(end)

# def match_flanking_exon(coord_str, pe_name, position_type, pe_strand, gtf_exons, 
#                         transcript_cds, gene_transcripts, gene_exons):
#     """
#     Match a single flanking exon coordinate to GTF.
    
#     Returns:
#         matched: list of dicts (one per GTF match)
#         unmatched: dict or None
#     """
#     chrom, start, end = parse_coordinate_string(coord_str)
    
#     matched = []
    
#     # Try to match with the given strand first, then try opposite
#     for test_strand in [pe_strand, '+' if pe_strand == '-' else '-']:
#         exon_key = (chrom, start, end, test_strand)
        
#         if exon_key in gtf_exons:
#             gtf_matches = gtf_exons[exon_key]
            
#             # Aggregate across all matching transcripts
#             transcript_ids = []
#             gene_names = set()
#             gene_ids = set()
#             gene_types = set()
#             tx_types = set()
#             feature_types = set()
            
#             phase_results = []
            
#             for match in gtf_matches:
#                 transcript_ids.append(match['transcript_id'])
#                 gene_names.add(match['gene_name'])
#                 gene_ids.add(match['gene_id'])
#                 gene_types.add(match['gene_type'])
#                 tx_types.add(match['transcript_type'])
#                 feature_types.add(match['feature_type'])
                
#                 # All matches are CDS, so always calculate phase
#                 phase_info = calculate_exon_phase(exon_key, match['transcript_id'], transcript_cds)
#                 if phase_info[0] is not None:
#                     phase_results.append(phase_info)
            
#             # Determine constitutive status
#             constitutive_status = []
#             for gene_name in gene_names:
#                 if gene_name in gene_transcripts and gene_name in gene_exons:
#                     total_transcripts = len(gene_transcripts[gene_name])
#                     transcripts_with_exon = len(gene_exons[gene_name][exon_key])
#                     is_constitutive = (transcripts_with_exon == total_transcripts)
#                     constitutive_status.append(is_constitutive)
#                 else:
#                     constitutive_status.append(False)
            
#             constitutive = all(constitutive_status) if constitutive_status else False
            
#             # Aggregate phase information (should always have phases since all are CDS)
#             if phase_results:
#                 phases = sorted(set([p[0] for p in phase_results]))
#                 cumulative_phases = sorted(set([p[1] for p in phase_results]))
#                 positions = sorted(set([p[2] for p in phase_results]))
#                 exon_nums = sorted(set([p[3] for p in phase_results]))
#                 total_exons = sorted(set([p[4] for p in phase_results]))
                
#                 # Always use comma-separated list, even if only one value
#                 phase = ','.join(map(str, phases))
#                 cumulative_phase = ','.join(map(str, cumulative_phases))
#                 exon_position = ','.join(positions)
#                 exon_num = ','.join(map(str, exon_nums))
#                 total_exon_count = ','.join(map(str, total_exons))
#                 all_phases = ','.join(map(str, phases))
#                 all_cumulative_phases = ','.join(map(str, cumulative_phases))
#                 phase_consistent = len(phases) == 1
#                 cumulative_phase_consistent = len(cumulative_phases) == 1
                
#                 exon_length = end - start
#                 symmetry = exon_length % 3
#             else:
#                 # Shouldn't happen for CDS, but handle gracefully
#                 phase = 'NA'
#                 cumulative_phase = 'NA'
#                 exon_position = 'NA'
#                 exon_num = 'NA'
#                 total_exon_count = 'NA'
#                 all_phases = 'NA'
#                 all_cumulative_phases = 'NA'
#                 phase_consistent = 'NA'
#                 cumulative_phase_consistent = 'NA'
#                 exon_length = end - start
#                 symmetry = 'NA'
            
#             matched.append({
#                 'chrom': chrom,
#                 'start': start,
#                 'end': end,
#                 'pe_name': pe_name,
#                 'position_relative_pe': position_type,
#                 'score': len(transcript_ids),
#                 'strand': test_strand,
#                 'length': exon_length,
#                 'symmetry': symmetry,
#                 'phase': phase,
#                 'all_phases': all_phases,
#                 'phase_consistent': phase_consistent,
#                 'cumulative_phase': cumulative_phase,
#                 'all_cumulative_phases': all_cumulative_phases,
#                 'cumulative_phase_consistent': cumulative_phase_consistent,
#                 'exon_position': exon_position,
#                 'exon_num': exon_num,
#                 'total_exons': total_exon_count,
#                 'constitutive': constitutive,
#                 'feature_types': ','.join(sorted(feature_types)),
#                 'tx_names': ','.join(transcript_ids),
#                 'n_transcripts': len(transcript_ids),
#                 'gene_name': ','.join(sorted(gene_names)),
#                 'gene_id': ','.join(sorted(gene_ids)),
#                 'gene_type': ','.join(sorted(gene_types)),
#                 'tx_type': ','.join(sorted(tx_types))
#             })
            
#             return matched, None
    
#     # No match found
#     unmatched = {
#         'chrom': chrom,
#         'start': start,
#         'end': end,
#         'pe_name': pe_name,
#         'position_relative_pe': position_type,
#         'strand': pe_strand
#     }
    
#     return [], unmatched

# def process_poison_exons(pe_file, gtf_exons, transcript_cds, gene_transcripts, gene_exons):
#     """
#     Process poison exons and match flanking exons to GTF.
    
#     Returns:
#         flanking_matched: list of matched flanking exons
#         flanking_unmatched: list of unmatched flanking exons
#         poison_exons_annotated: list of poison exons with phase info
#     """
#     print(f"\nReading poison exons: {pe_file}")
#     pe_df = pd.read_csv(pe_file, sep='\t')
    
#     print(f"  Total poison exons: {len(pe_df)}")
    
#     flanking_matched = []
#     flanking_unmatched = []
#     poison_exons_annotated = []
    
#     for idx, row in pe_df.iterrows():
#         if (idx + 1) % 1000 == 0:
#             print(f"  Processed {idx+1}/{len(pe_df)} poison exons...")
        
#         pe_name = row['name']
#         pe_chrom = row['chrom']
#         pe_start = row['start']
#         pe_end = row['end']
#         pe_strand = row['strand']
#         pe_length = row['length']
        
#         upstream_exon = row['upstream_exon']
#         downstream_exon = row['downstream_exon']
        
#         # Extract gene name from poison exon name (format: chr:start-end_GENE)
#         if '_' in pe_name:
#             pe_gene_name = pe_name.split('_')[-1]
#         else:
#             pe_gene_name = 'NA'
        
#         # Process upstream flanking exons
#         upstream_cumulative_phases = []
#         upstream_gene_names = set()
#         upstream_gene_ids = set()
        
#         if upstream_exon != 'NA':
#             for coord in upstream_exon.split(','):
#                 coord = coord.strip()
#                 matched, unmatched = match_flanking_exon(
#                     coord, pe_name, 'upstream_flank', pe_strand,
#                     gtf_exons, transcript_cds, gene_transcripts, gene_exons
#                 )
                
#                 if matched:
#                     flanking_matched.extend(matched)
#                     for m in matched:
#                         if m['cumulative_phase'] != 'NA':
#                             # Parse comma-separated phases
#                             phases = m['cumulative_phase'].split(',')
#                             for phase_str in phases:
#                                 upstream_cumulative_phases.append(int(phase_str))
#                         if m['gene_name'] != 'NA':
#                             upstream_gene_names.update(m['gene_name'].split(','))
#                         if m['gene_id'] != 'NA':
#                             upstream_gene_ids.update(m['gene_id'].split(','))
#                 else:
#                     flanking_unmatched.append(unmatched)
        
#         # Process downstream flanking exons
#         downstream_cumulative_phases = []
#         downstream_gene_names = set()
#         downstream_gene_ids = set()
        
#         if downstream_exon != 'NA':
#             for coord in downstream_exon.split(','):
#                 coord = coord.strip()
#                 matched, unmatched = match_flanking_exon(
#                     coord, pe_name, 'downstream_flank', pe_strand,
#                     gtf_exons, transcript_cds, gene_transcripts, gene_exons
#                 )
                
#                 if matched:
#                     flanking_matched.extend(matched)
#                     for m in matched:
#                         if m['cumulative_phase'] != 'NA':
#                             # Parse comma-separated phases
#                             phases = m['cumulative_phase'].split(',')
#                             for phase_str in phases:
#                                 downstream_cumulative_phases.append(int(phase_str))
#                         if m['gene_name'] != 'NA':
#                             downstream_gene_names.update(m['gene_name'].split(','))
#                         if m['gene_id'] != 'NA':
#                             downstream_gene_ids.update(m['gene_id'].split(','))
#                 else:
#                     flanking_unmatched.append(unmatched)
        
#         # Calculate poison exon phase (strand-aware)
#         if pe_strand == '+':
#             # Plus strand: upstream exon comes first
#             relevant_cumulative_phases = upstream_cumulative_phases
#         else:
#             # Minus strand: downstream exon comes first (biologically)
#             relevant_cumulative_phases = downstream_cumulative_phases
        
#         # Determine poison exon phase
#         if relevant_cumulative_phases:
#             # Use cumulative phases from preceding flanking exon(s)
#             pe_phases = sorted(set(relevant_cumulative_phases))
#             pe_phase = ','.join(map(str, pe_phases))
#         else:
#             pe_phase = 'No Flanking GTF match'
        
#         # Calculate cumulative phase for poison exon
#         pe_symmetry = pe_length % 3
        
#         if pe_phase != 'No Flanking GTF match':
#             # Calculate cumulative phase for each starting phase
#             pe_cumulative_phases = []
#             for p in pe_phases:
#                 cum_p = (p + pe_length) % 3
#                 pe_cumulative_phases.append(cum_p)
#             pe_cumulative_phase = ','.join(map(str, sorted(set(pe_cumulative_phases))))
#         else:
#             pe_cumulative_phase = 'No Flanking GTF match'
        
#         # Collect all gene names
#         all_gene_names = set([pe_gene_name]) | upstream_gene_names | downstream_gene_names
#         all_gene_names.discard('NA')
#         gene_name_str = ','.join(sorted(all_gene_names)) if all_gene_names else pe_gene_name
        
#         all_gene_ids = upstream_gene_ids | downstream_gene_ids
#         all_gene_ids.discard('NA')
#         gene_id_str = ','.join(sorted(all_gene_ids)) if all_gene_ids else 'NA'
        
#         poison_exons_annotated.append({
#             'chrom': pe_chrom,
#             'start': pe_start,
#             'end': pe_end,
#             'name': pe_name,
#             'strand': pe_strand,
#             'length': pe_length,
#             'symmetry': pe_symmetry,
#             'phase': pe_phase,
#             'cumulative_phase': pe_cumulative_phase,
#             'gene_name': gene_name_str,
#             'gene_id': gene_id_str,
#             'feature_types': 'poisonExon',
#             'tx_type': 'poisonExon'
#         })
    
#     print(f"\nProcessing complete:")
#     print(f"  Flanking exons matched: {len(flanking_matched)}")
#     print(f"  Flanking exons unmatched: {len(flanking_unmatched)}")
#     print(f"  Poison exons annotated: {len(poison_exons_annotated)}")
    
#     return flanking_matched, flanking_unmatched, poison_exons_annotated

# def write_flanking_matched(matched_list, output_file):
#     """Write matched flanking exons to file."""
#     df = pd.DataFrame(matched_list)
    
#     if len(df) == 0:
#         print(f"WARNING: No matched flanking exons")
#         with open(output_file, 'w') as f:
#             f.write("# No matched flanking exons\n")
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
#         f.write("# Flanking exons matched to GTF\n")
#         f.write("# Columns: chrom, start, end, pe_name, position_relative_pe, score, strand, length, symmetry, ")
#         f.write("phase, all_phases, phase_consistent, cumulative_phase, all_cumulative_phases, ")
#         f.write("cumulative_phase_consistent, exon_position, exon_num, total_exons, constitutive, ")
#         f.write("feature_types, tx_names, n_transcripts, gene_name, gene_id, gene_type, tx_type\n")
#         f.write("# Coordinates: 0-based half-open (BED format)\n")
        
#         for _, row in df.iterrows():
#             f.write(f"{row['chrom']}\t{row['start']}\t{row['end']}\t{row['pe_name']}\t")
#             f.write(f"{row['position_relative_pe']}\t{row['score']}\t{row['strand']}\t")
#             f.write(f"{row['length']}\t{row['symmetry']}\t{row['phase']}\t{row['all_phases']}\t")
#             f.write(f"{row['phase_consistent']}\t{row['cumulative_phase']}\t")
#             f.write(f"{row['all_cumulative_phases']}\t{row['cumulative_phase_consistent']}\t")
#             f.write(f"{row['exon_position']}\t{row['exon_num']}\t{row['total_exons']}\t")
#             f.write(f"{row['constitutive']}\t{row['feature_types']}\t{row['tx_names']}\t")
#             f.write(f"{row['n_transcripts']}\t{row['gene_name']}\t{row['gene_id']}\t")
#             f.write(f"{row['gene_type']}\t{row['tx_type']}\n")
    
#     print(f"\nWrote {len(df)} matched flanking exons to: {output_file}")

# def write_flanking_unmatched(unmatched_list, output_file):
#     """Write unmatched flanking exons to file."""
#     df = pd.DataFrame(unmatched_list)
    
#     if len(df) == 0:
#         print(f"INFO: All flanking exons matched!")
#         with open(output_file, 'w') as f:
#             f.write("# All flanking exons matched\n")
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
#         f.write("# Flanking exons with no GTF match\n")
#         f.write("# Columns: chrom, start, end, pe_name, position_relative_pe, strand\n")
        
#         for _, row in df.iterrows():
#             f.write(f"{row['chrom']}\t{row['start']}\t{row['end']}\t")
#             f.write(f"{row['pe_name']}\t{row['position_relative_pe']}\t{row['strand']}\n")
    
#     print(f"\nWrote {len(df)} unmatched flanking exons to: {output_file}")

# def write_poison_exons(pe_list, output_file):
#     """Write annotated poison exons to file."""
#     df = pd.DataFrame(pe_list)
    
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
#         f.write("# Poison exons with phase annotations from flanking exons\n")
#         f.write("# Columns: chrom, start, end, name, strand, length, symmetry, phase, ")
#         f.write("cumulative_phase, gene_name, gene_id, feature_types, tx_type\n")
#         f.write("# Coordinates: 0-based half-open (BED format)\n")
#         f.write("# Phase: Inherited from preceding flanking exon (strand-aware)\n")
#         f.write("# Cumulative_phase: (phase + length) % 3\n")
        
#         for _, row in df.iterrows():
#             f.write(f"{row['chrom']}\t{row['start']}\t{row['end']}\t{row['name']}\t")
#             f.write(f"{row['strand']}\t{row['length']}\t{row['symmetry']}\t")
#             f.write(f"{row['phase']}\t{row['cumulative_phase']}\t{row['gene_name']}\t")
#             f.write(f"{row['gene_id']}\t{row['feature_types']}\t{row['tx_type']}\n")
    
#     print(f"\nWrote {len(df)} poison exons to: {output_file}")

# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(
#         description='Filter poison exon flanking output and match to GTF',
#         formatter_class=argparse.RawDescriptionHelpFormatter,
#         epilog="""
# Processes poison exon flanking output to:
# 1. Match flanking exons to GTF annotations
# 2. Annotate poison exons with phase information from flanking exons

# Output Files:
#   1. *_flanking_gtfMatch.bed: Flanking exons matched to GTF
#   2. *_flanking_gtfNoMatch.bed: Flanking exons not in GTF
#   3. *_poisonExons_annotatedFlanks.bed: Poison exons with phase info

# Example:
#   python filter_poison_exon_flanking_output.py \\
#       -p poison_exons_with_flanks.tsv \\
#       -g reference.gtf \\
#       -o output/filtered
#         """
#     )
    
#     parser.add_argument('-p', '--poison-exons', required=True,
#                        help='Input poison exon TSV file with flanking exons')
#     parser.add_argument('-g', '--gtf', required=True,
#                        help='Reference GTF file')
#     parser.add_argument('-o', '--output', required=True,
#                        help='Output prefix for result files')
    
#     args = parser.parse_args()
    
#     print("="*70)
#     print("FILTER POISON EXON FLANKING OUTPUT")
#     print("="*70)
#     print(f"Poison exons: {args.poison_exons}")
#     print(f"GTF: {args.gtf}")
#     print(f"Output prefix: {args.output}")
#     print("="*70 + "\n")
    
#     # Load GTF
#     gtf_exons, gene_transcripts, gene_exons = load_gtf_exons(args.gtf)
#     transcript_cds = extract_cds_from_gtf(args.gtf)
    
#     # Process poison exons and flanking exons
#     flanking_matched, flanking_unmatched, poison_exons = process_poison_exons(
#         args.poison_exons, gtf_exons, transcript_cds, gene_transcripts, gene_exons
#     )
    
#     # Write outputs
#     print("\n" + "="*70)
#     print("WRITING OUTPUT FILES")
#     print("="*70)
    
#     flanking_match_file = f"{args.output}_flanking_gtfMatch.bed"
#     write_flanking_matched(flanking_matched, flanking_match_file)
    
#     print("\n" + "-"*70)
    
#     flanking_nomatch_file = f"{args.output}_flanking_gtfNoMatch.bed"
#     write_flanking_unmatched(flanking_unmatched, flanking_nomatch_file)
    
#     print("\n" + "-"*70)
    
#     pe_file = f"{args.output}_poisonExons_annotatedFlanks.bed"
#     write_poison_exons(poison_exons, pe_file)
    
#     print("\n" + "="*70)
#     print("PROCESSING COMPLETE")
#     print("="*70)
#     print(f"\nOutput files:")
#     print(f"  1. {flanking_match_file}")
#     print(f"  2. {flanking_nomatch_file}")
#     print(f"  3. {pe_file}")



# #
## Monday 

#!/usr/bin/env python3
"""
Filter poison exon flanking output and match flanking exons to GTF annotations.

Creates three output files:
1. Flanking exons matched to GTF with full annotation
2. Flanking exons not matched to GTF
3. Poison exons with phase information from flanking exons

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

def load_gtf_exons(gtf_file):
    """
    Load CDS coordinates from GTF file for flanking exon matching.
    """
    print(f"Loading GTF file (CDS only): {gtf_file}")
    
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
            
            if feature_type != 'CDS':
                continue
            
            chrom = fields[0]
            start_1based = int(fields[3])
            end_1based = int(fields[4])
            strand = fields[6]
            
            start_0based = start_1based - 1
            end_0based = end_1based
            
            attributes = parse_gtf_attributes(fields[8])
            
            gene_name = attributes.get('gene_name', 'NA')
            gene_id = attributes.get('gene_id', 'NA')
            transcript_id = attributes.get('transcript_id', 'NA')
            
            exon_key = (chrom, start_0based, end_0based, strand)
            
            exon_records[exon_key].append({
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
            })
            
            if gene_name != 'NA' and transcript_id != 'NA':
                gene_transcripts[gene_name].add(transcript_id)
                gene_exons[gene_name][exon_key].add(transcript_id)
            
            if line_num % 100000 == 0:
                print(f"  Processed {line_num} GTF lines...")
    
    print(f"Loaded {len(exon_records)} unique CDS coordinates")
    print(f"Tracked {len(gene_transcripts)} genes for constitutive analysis")
    return exon_records, gene_transcripts, gene_exons

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
    
    exon_idx = None
    for i, cds in enumerate(cds_list):
        if cds['chrom'] == chrom and cds['start'] == start and cds['end'] == end:
            exon_idx = i
            break
    
    if exon_idx is None:
        return None, None, None, None, None
    
    if exon_idx == 0:
        phase = 0
    else:
        cumulative_length = 0
        for i in range(exon_idx):
            upstream_cds = cds_list[i]
            upstream_length = upstream_cds['end'] - upstream_cds['start']
            cumulative_length += upstream_length
        phase = cumulative_length % 3
    
    current_length = end - start
    cumulative_phase = (phase + current_length) % 3
    
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

def parse_coordinate_string(coord_str):
    """Parse coordinate string chr:start-end into (chr, start, end)."""
    chrom, pos = coord_str.split(':')
    start, end = pos.split('-')
    return chrom, int(start), int(end)

def match_flanking_exon(coord_str, pe_name, position_type, pe_strand, gtf_exons, 
                        transcript_cds, gene_transcripts, gene_exons):
    """Match a single flanking exon coordinate to GTF."""
    chrom, start, end = parse_coordinate_string(coord_str)
    
    matched = []
    
    for test_strand in [pe_strand, '+' if pe_strand == '-' else '-']:
        exon_key = (chrom, start, end, test_strand)
        
        if exon_key in gtf_exons:
            gtf_matches = gtf_exons[exon_key]
            
            transcript_ids = []
            gene_names = set()
            gene_ids = set()
            gene_types = set()
            tx_types = set()
            feature_types = set()
            
            phase_results = []
            
            for match in gtf_matches:
                transcript_ids.append(match['transcript_id'])
                gene_names.add(match['gene_name'])
                gene_ids.add(match['gene_id'])
                gene_types.add(match['gene_type'])
                tx_types.add(match['transcript_type'])
                feature_types.add(match['feature_type'])
                
                phase_info = calculate_exon_phase(exon_key, match['transcript_id'], transcript_cds)
                if phase_info[0] is not None:
                    phase_results.append(phase_info)
            
            constitutive_status = []
            for gene_name in gene_names:
                if gene_name in gene_transcripts and gene_name in gene_exons:
                    total_transcripts = len(gene_transcripts[gene_name])
                    transcripts_with_exon = len(gene_exons[gene_name][exon_key])
                    is_constitutive = (transcripts_with_exon == total_transcripts)
                    constitutive_status.append(is_constitutive)
                else:
                    constitutive_status.append(False)
            
            constitutive = all(constitutive_status) if constitutive_status else False
            
            exon_length = end - start
            
            if phase_results:
                phases = sorted(set([p[0] for p in phase_results]))
                cumulative_phases = sorted(set([p[1] for p in phase_results]))
                positions = sorted(set([p[2] for p in phase_results]))
                exon_nums = sorted(set([p[3] for p in phase_results]))
                total_exons = sorted(set([p[4] for p in phase_results]))
                
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
            
            matched.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'name': f"{pe_name}_{position_type}",
                'strand': test_strand,
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
                'gene_name': ','.join(sorted(gene_names)),
                'gene_id': ','.join(sorted(gene_ids)),
                'gene_type': ','.join(sorted(gene_types)),
                'tx_type': ','.join(sorted(tx_types)),
                'pe_name': pe_name,
                'position_relative_pe': position_type
            })
            
            return matched, None
    
    unmatched = {
        'chrom': chrom,
        'start': start,
        'end': end,
        'pe_name': pe_name,
        'position_relative_pe': position_type,
        'strand': pe_strand
    }
    
    return [], unmatched

def process_poison_exons(pe_file, gtf_exons, transcript_cds, gene_transcripts, gene_exons):
    """Process poison exons and match flanking exons to GTF."""
    print(f"\nReading poison exons: {pe_file}")
    pe_df = pd.read_csv(pe_file, sep='\t')
    
    print(f"  Total poison exons: {len(pe_df)}")
    
    flanking_matched = []
    flanking_unmatched = []
    poison_exons_annotated = []
    
    for idx, row in pe_df.iterrows():
        if (idx + 1) % 1000 == 0:
            print(f"  Processed {idx+1}/{len(pe_df)} poison exons...")
        
        pe_name = row['name']
        pe_chrom = row['chrom']
        pe_start = row['start']
        pe_end = row['end']
        pe_strand = row['strand']
        pe_length = row['length']
        
        upstream_exon = row['upstream_exon']
        downstream_exon = row['downstream_exon']
        
        if '_' in pe_name:
            pe_gene_name = pe_name.split('_')[-1]
        else:
            pe_gene_name = 'NA'
        
        upstream_cumulative_phases = []
        upstream_gene_names = set()
        upstream_gene_ids = set()
        
        if upstream_exon != 'NA':
            for coord in upstream_exon.split(','):
                coord = coord.strip()
                matched, unmatched = match_flanking_exon(
                    coord, pe_name, 'upstream_flank', pe_strand,
                    gtf_exons, transcript_cds, gene_transcripts, gene_exons
                )
                
                if matched:
                    flanking_matched.extend(matched)
                    for m in matched:
                        if m['cumulative_phase'] != 'NA':
                            phases = m['cumulative_phase'].split(',')
                            for phase_str in phases:
                                upstream_cumulative_phases.append(int(phase_str))
                        if m['gene_name'] != 'NA':
                            upstream_gene_names.update(m['gene_name'].split(','))
                        if m['gene_id'] != 'NA':
                            upstream_gene_ids.update(m['gene_id'].split(','))
                else:
                    flanking_unmatched.append(unmatched)
        
        downstream_cumulative_phases = []
        downstream_gene_names = set()
        downstream_gene_ids = set()
        
        if downstream_exon != 'NA':
            for coord in downstream_exon.split(','):
                coord = coord.strip()
                matched, unmatched = match_flanking_exon(
                    coord, pe_name, 'downstream_flank', pe_strand,
                    gtf_exons, transcript_cds, gene_transcripts, gene_exons
                )
                
                if matched:
                    flanking_matched.extend(matched)
                    for m in matched:
                        if m['cumulative_phase'] != 'NA':
                            phases = m['cumulative_phase'].split(',')
                            for phase_str in phases:
                                downstream_cumulative_phases.append(int(phase_str))
                        if m['gene_name'] != 'NA':
                            downstream_gene_names.update(m['gene_name'].split(','))
                        if m['gene_id'] != 'NA':
                            downstream_gene_ids.update(m['gene_id'].split(','))
                else:
                    flanking_unmatched.append(unmatched)
        
        if pe_strand == '+':
            relevant_cumulative_phases = upstream_cumulative_phases
        else:
            relevant_cumulative_phases = downstream_cumulative_phases
        
        if relevant_cumulative_phases:
            pe_phases = sorted(set(relevant_cumulative_phases))
            pe_phase = ','.join(map(str, pe_phases))
            phase_consistent = len(pe_phases) == 1
        else:
            pe_phase = 'NA'
            pe_phases = []
            phase_consistent = 'NA'
        
        pe_symmetry = pe_length % 3
        
        if pe_phase != 'NA':
            pe_cumulative_phases = []
            for p in pe_phases:
                cum_p = (p + pe_length) % 3
                pe_cumulative_phases.append(cum_p)
            pe_cumulative_phase = ','.join(map(str, sorted(set(pe_cumulative_phases))))
            cumulative_phase_consistent = len(set(pe_cumulative_phases)) == 1
        else:
            pe_cumulative_phase = 'NA'
            cumulative_phase_consistent = 'NA'
        
        all_gene_names = set([pe_gene_name]) | upstream_gene_names | downstream_gene_names
        all_gene_names.discard('NA')
        gene_name_str = ','.join(sorted(all_gene_names)) if all_gene_names else pe_gene_name
        
        all_gene_ids = upstream_gene_ids | downstream_gene_ids
        all_gene_ids.discard('NA')
        gene_id_str = ','.join(sorted(all_gene_ids)) if all_gene_ids else 'NA'
        
        # STANDARDIZED OUTPUT for poison exons
        # FIXED: gene_type = 'poisonExon', exon_position = 'internal'
        poison_exons_annotated.append({
            'chrom': pe_chrom,
            'start': pe_start,
            'end': pe_end,
            'name': pe_name,
            'strand': pe_strand,
            'length': pe_length,
            'symmetry': pe_symmetry,
            'phase': pe_phase,
            'all_phases': pe_phase,
            'phase_consistent': phase_consistent,
            'cumulative_phase': pe_cumulative_phase,
            'all_cumulative_phases': pe_cumulative_phase,
            'cumulative_phase_consistent': cumulative_phase_consistent,
            'exon_position': 'internal',  # FIXED: Poison exons are always internal
            'exon_num': 'NA',
            'total_exons': 'NA',
            'constitutive': False,
            'feature_types': 'poisonExon',
            'tx_names': 'NA',
            'n_transcripts': 0,
            'gene_name': gene_name_str,
            'gene_id': gene_id_str,
            'gene_type': 'poisonExon',  # FIXED: gene_type = 'poisonExon'
            'tx_type': 'poisonExon'
        })
    
    print(f"\nProcessing complete:")
    print(f"  Flanking exons matched: {len(flanking_matched)}")
    print(f"  Flanking exons unmatched: {len(flanking_unmatched)}")
    print(f"  Poison exons annotated: {len(poison_exons_annotated)}")
    
    return flanking_matched, flanking_unmatched, poison_exons_annotated

def write_flanking_matched(matched_list, output_file):
    """Write matched flanking exons to file with STANDARDIZED columns."""
    df = pd.DataFrame(matched_list)
    
    if len(df) == 0:
        print(f"WARNING: No matched flanking exons")
        with open(output_file, 'w') as f:
            f.write("# No matched flanking exons\n")
        return
    
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
    
    standard_cols = [
        'chrom', 'start', 'end', 'name', 'strand', 'length', 'symmetry',
        'phase', 'all_phases', 'phase_consistent',
        'cumulative_phase', 'all_cumulative_phases', 'cumulative_phase_consistent',
        'exon_position', 'exon_num', 'total_exons', 'constitutive',
        'feature_types', 'tx_names', 'n_transcripts',
        'gene_name', 'gene_id', 'gene_type', 'tx_type',
        'pe_name', 'position_relative_pe'
    ]
    
    df = df[standard_cols]
    
    with open(output_file, 'w') as f:
        f.write("# Flanking exons matched to GTF\n")
        f.write("# STANDARDIZED COLUMNS: " + ", ".join(standard_cols) + "\n")
        f.write("# Coordinates: 0-based half-open (BED format)\n")
        
        f.write('\t'.join(standard_cols) + '\n')
        
        for _, row in df.iterrows():
            values = [str(row[col]) for col in standard_cols]
            f.write('\t'.join(values) + '\n')
    
    print(f"\nWrote {len(df)} matched flanking exons to: {output_file}")

def write_flanking_unmatched(unmatched_list, output_file):
    """Write unmatched flanking exons to file."""
    df = pd.DataFrame(unmatched_list)
    
    if len(df) == 0:
        print(f"INFO: All flanking exons matched!")
        with open(output_file, 'w') as f:
            f.write("# All flanking exons matched\n")
        return
    
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
    
    with open(output_file, 'w') as f:
        f.write("# Flanking exons with no GTF match\n")
        f.write("# Columns: chrom, start, end, pe_name, position_relative_pe, strand\n")
        
        for _, row in df.iterrows():
            f.write(f"{row['chrom']}\t{row['start']}\t{row['end']}\t")
            f.write(f"{row['pe_name']}\t{row['position_relative_pe']}\t{row['strand']}\n")
    
    print(f"\nWrote {len(df)} unmatched flanking exons to: {output_file}")

def write_poison_exons(pe_list, output_file):
    """Write annotated poison exons to file with STANDARDIZED columns."""
    df = pd.DataFrame(pe_list)
    
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
    
    standard_cols = [
        'chrom', 'start', 'end', 'name', 'strand', 'length', 'symmetry',
        'phase', 'all_phases', 'phase_consistent',
        'cumulative_phase', 'all_cumulative_phases', 'cumulative_phase_consistent',
        'exon_position', 'exon_num', 'total_exons', 'constitutive',
        'feature_types', 'tx_names', 'n_transcripts',
        'gene_name', 'gene_id', 'gene_type', 'tx_type'
    ]
    
    df = df[standard_cols]
    
    with open(output_file, 'w') as f:
        f.write("# Poison exons with phase annotations from flanking exons\n")
        f.write("# STANDARDIZED COLUMNS: " + ", ".join(standard_cols) + "\n")
        f.write("# Coordinates: 0-based half-open (BED format)\n")
        f.write("# Phase: Inherited from preceding flanking exon (strand-aware)\n")
        f.write("# Cumulative_phase: (phase + length) % 3\n")
        f.write("# exon_position: 'internal' (poison exons are always internal by definition)\n")
        f.write("# gene_type: 'poisonExon'\n")
        
        f.write('\t'.join(standard_cols) + '\n')
        
        for _, row in df.iterrows():
            values = [str(row[col]) for col in standard_cols]
            f.write('\t'.join(values) + '\n')
    
    print(f"\nWrote {len(df)} poison exons to: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Filter poison exon flanking output and match to GTF',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Processes poison exon flanking output to:
1. Match flanking exons to GTF annotations
2. Annotate poison exons with phase information from flanking exons

STANDARDIZED OUTPUT COLUMNS:
chrom, start, end, name, strand, length, symmetry, phase, all_phases, phase_consistent,
cumulative_phase, all_cumulative_phases, cumulative_phase_consistent, exon_position,
exon_num, total_exons, constitutive, feature_types, tx_names, n_transcripts, gene_name,
gene_id, gene_type, tx_type

FIXES in this version:
- gene_type = 'poisonExon' (not NA)
- exon_position = 'internal' (not NA)

Output Files:
  1. *_flanking_gtfMatch.bed: Flanking exons matched to GTF
  2. *_flanking_gtfNoMatch.bed: Flanking exons not in GTF
  3. *_poisonExons_annotatedFlanks.bed: Poison exons with phase info

Example:
  python filter_poison_exon_flanking_output.py \\
      -p poison_exons_with_flanks.tsv \\
      -g reference.gtf \\
      -o output/filtered
        """
    )
    
    parser.add_argument('-p', '--poison-exons', required=True,
                       help='Input poison exon TSV file with flanking exons')
    parser.add_argument('-g', '--gtf', required=True,
                       help='Reference GTF file')
    parser.add_argument('-o', '--output', required=True,
                       help='Output prefix for result files')
    
    args = parser.parse_args()
    
    print("="*70)
    print("FILTER POISON EXON FLANKING OUTPUT")
    print("="*70)
    print(f"Poison exons: {args.poison_exons}")
    print(f"GTF: {args.gtf}")
    print(f"Output prefix: {args.output}")
    print("="*70 + "\n")
    
    gtf_exons, gene_transcripts, gene_exons = load_gtf_exons(args.gtf)
    transcript_cds = extract_cds_from_gtf(args.gtf)
    
    flanking_matched, flanking_unmatched, poison_exons = process_poison_exons(
        args.poison_exons, gtf_exons, transcript_cds, gene_transcripts, gene_exons
    )
    
    print("\n" + "="*70)
    print("WRITING OUTPUT FILES")
    print("="*70)
    
    flanking_match_file = f"{args.output}_flanking_gtfMatch.bed"
    write_flanking_matched(flanking_matched, flanking_match_file)
    
    print("\n" + "-"*70)
    
    flanking_nomatch_file = f"{args.output}_flanking_gtfNoMatch.bed"
    write_flanking_unmatched(flanking_unmatched, flanking_nomatch_file)
    
    print("\n" + "-"*70)
    
    pe_file = f"{args.output}_poisonExons_annotatedFlanks.bed"
    write_poison_exons(poison_exons, pe_file)
    
    print("\n" + "="*70)
    print("PROCESSING COMPLETE")
    print("="*70)
    print(f"\nOutput files:")
    print(f"  1. {flanking_match_file}")
    print(f"  2. {flanking_nomatch_file}")
    print(f"  3. {pe_file}")