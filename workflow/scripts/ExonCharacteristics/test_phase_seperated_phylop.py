# Author: Dylan Stermer
# Date: December 2024
# Updated: January 2026 - Added phase-specific output
# """

import pyBigWig
import numpy as np
import pandas as pd
from scipy import stats
import sys
import os
import argparse

def get_full_exon_conservation(bw, chrom, start, end):
    """
    Get mean conservation score across the ENTIRE exon.
    Completely independent of edge calculations.
    
    Args:
        bw: pyBigWig object
        chrom: Chromosome
        start: Exon start coordinate
        end: Exon end coordinate
    
    Returns:
        float: Mean conservation score across entire exon, or np.nan if failed
    """
    try:
        scores = np.array(bw.values(chrom, start, end))
        scores = np.nan_to_num(scores, nan=0.0)
        return np.mean(scores)
    except (RuntimeError, Exception):
        return np.nan

def fix_bigwig_path(path):
    """
    pyBigWig has issues with absolute paths in some environments.
    Convert to relative if possible.
    """
    if os.path.isabs(path):
        try:
            rel_path = os.path.relpath(path)
            if len(rel_path) < len(path):
                return rel_path
        except ValueError:
            pass
    return path

def parse_bed12_to_internal_exons(bed_file, chrom_filter=None, min_length=50):
    """
    Extract exons from BED file.
    Handles both standard BED12 and custom BED format from extract_constitutive_alternative.py
    Filter by minimum exon length and protein_coding, poisonExon, or lncRNA type.
    
    Args:
        bed_file: Path to BED file
        chrom_filter: Chromosome to analyze (None = all chromosomes)
        min_length: Minimum exon length in nt (default: 50)
    
    Returns:
        DataFrame with exons
    """
    
    # Try to read as TSV first (for standardized files with headers)
    try:
        df = pd.read_csv(bed_file, sep='\t', comment='#')
        
        # Check if it has the standard columns
        required_cols = ['chrom', 'start', 'end', 'strand', 'length', 'gene_name', 'gene_type', 'tx_type']
        if all(col in df.columns for col in required_cols):
            print(f"Detected standardized TSV format with {len(df.columns)} columns")
            
            # Apply filters
            if chrom_filter:
                df = df[df['chrom'] == chrom_filter]
            
            # Filter by gene/tx type
            df = df[
                (df['gene_type'].isin(['protein_coding', 'poisonExon', 'lncRNA'])) &
                (df['tx_type'].isin(['protein_coding', 'poisonExon', 'lncRNA']))
            ]
            
            # Filter by minimum length
            df = df[df['length'] >= min_length]
            
            # Add required columns if missing
            if 'tx_name' not in df.columns:
                df['tx_name'] = df.get('name', 'NA')
            if 'exon_num' not in df.columns:
                df['exon_num'] = 0
            if 'total_exons' not in df.columns:
                df['total_exons'] = 0
            
            print(f"Loaded {len(df)} exons from standardized TSV")
            return df
    except Exception as e:
        print(f"Not a TSV with headers, trying BED format parsing: {e}")
    
    # Fall back to manual parsing for BED files
    exons = []
    skipped_non_target = 0
    
    with open(bed_file, 'r') as f:
        for line in f:
            # Skip comments
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            
            # Skip header row (detect if first field is "chrom" or similar text)
            if len(fields) > 0 and fields[0].lower() in ['chrom', 'chromosome', 'chr']:
                continue
            
            if len(fields) < 12:
                print(f"WARNING: Line has only {len(fields)} columns, skipping", file=sys.stderr)
                continue
            
            chrom = fields[0]
            
            if chrom_filter is not None and chrom != chrom_filter:
                continue    
            
            # Check if this is custom format by testing column 10
            # Standard BED12: column 10 is blockSizes (comma-separated numbers)
            # Custom format: column 10 is exon_position (text: "first", "internal", "last")
            try:
                # Try to parse column 10 as integer(s)
                test_val = fields[10].rstrip(',').split(',')[0]
                int(test_val)
                is_bed12 = True
            except ValueError:
                # Column 10 is text, must be custom format
                is_bed12 = False
            
            if is_bed12:
                # Standard BED12 format
                tx_start = int(fields[1])
                tx_name = fields[3]
                strand = fields[5]
                exon_count = int(fields[9])
                exon_sizes = [int(x) for x in fields[10].rstrip(',').split(',')]
                exon_starts = [int(x) for x in fields[11].rstrip(',').split(',')]
                
                # Extract metadata (columns 12+)
                gene_name = fields[12] if len(fields) > 12 else 'NA'
                gene_type = fields[13] if len(fields) > 13 else 'NA'
                tx_type = fields[14] if len(fields) > 14 else 'NA'
                
                # Filter for protein_coding, poisonExon, or lncRNA
                if gene_type not in ['protein_coding', 'poisonExon', 'lncRNA'] or tx_type not in ['protein_coding', 'poisonExon', 'lncRNA']:
                    skipped_non_target += 1
                    continue
                
                # Skip transcripts with only 1 or 2 exons (no internal exons)
                if exon_count <= 2:
                    continue
                
                # Extract internal exons only (skip first i=0 and last i=exon_count-1)
                for i in range(1, exon_count - 1):
                    exon_start = tx_start + exon_starts[i]
                    exon_end = exon_start + exon_sizes[i]
                    exon_length = exon_sizes[i]
                    
                    # Filter by minimum length
                    if exon_length < min_length:
                        continue
                    
                    exons.append({
                        'chrom': chrom,
                        'start': exon_start,
                        'end': exon_end,
                        'length': exon_length,
                        'tx_name': tx_name,
                        'gene_name': gene_name,
                        'gene_type': gene_type,
                        'tx_type': tx_type,
                        'strand': strand,
                        'exon_num': i + 1,
                        'total_exons': exon_count
                    })
            else:
                # Custom format (standardized files or extract_constitutive_alternative.py)
                # Standardized format: chrom, start, end, name, strand, length, gene_name, gene_type, tx_type, exon_position
                exon_start = int(fields[1])
                exon_end = int(fields[2])
                exon_name = fields[3]
                strand = fields[4]
                exon_length = int(fields[5])
                gene_name = fields[6] if len(fields) > 6 else 'NA'
                gene_type = fields[7] if len(fields) > 7 else 'NA'
                tx_type = fields[8] if len(fields) > 8 else 'NA'
                exon_position = fields[9] if len(fields) > 9 else 'NA'
                
                # Try to extract phase if available (might be in different columns depending on format)
                phase = 'NA'
                for idx in [14, 15, 16]:  # Common positions for phase column
                    if len(fields) > idx:
                        try:
                            # Check if this looks like a phase value (number or comma-separated numbers)
                            test_phase = str(fields[idx])
                            if test_phase in ['0', '1', '2', '0,1', '0,2', '1,2', '0,1,2'] or test_phase.replace(',', '').isdigit():
                                phase = test_phase
                                break
                        except:
                            pass
                
                # Filter for protein_coding, poisonExon, or lncRNA
                if gene_type not in ['protein_coding', 'poisonExon', 'lncRNA'] or tx_type not in ['protein_coding', 'poisonExon', 'lncRNA']:
                    skipped_non_target += 1
                    continue
                
                # Filter by minimum length
                if exon_length < min_length:
                    continue
                
                exons.append({
                    'chrom': chrom,
                    'start': exon_start,
                    'end': exon_end,
                    'length': exon_length,
                    'tx_name': exon_name,
                    'gene_name': gene_name,
                    'gene_type': gene_type,
                    'tx_type': tx_type,
                    'strand': strand,
                    'exon_num': 0,
                    'total_exons': 0,
                    'phase': phase  # Add phase to the dictionary
                })
    
    if skipped_non_target > 0:
        print(f"Filtered out {skipped_non_target} entries (not protein_coding, poisonExon, or lncRNA)")
    
    return pd.DataFrame(exons)

def get_exon_and_flanking_conservation(bw, chrom, start, end, strand, 
                                       flank_size=50, exon_edge_bp=39, exon_phase=None):
    """
    Get conservation scores for flanking introns AND exonic edges.
    
    Returns conservation for:
    1. 5' flanking intron (50nt)
    2. First edge of exon (phase-adjusted to end at codon boundary)
    3. Last edge of exon (phase-adjusted to start at codon boundary)
    4. 3' flanking intron (50nt)
    
    Strand-aware coordinates.
    Phase-aware edge sizes for proper codon alignment.
    
    Args:
        bw: pyBigWig object
        chrom: Chromosome
        start: Exon start coordinate
        end: Exon end coordinate
        strand: '+' or '-'
        flank_size: Size of flanking intronic regions in nt (default: 50)
        exon_edge_bp: Base exon edge size in bp (default: 39)
        exon_phase: Reading frame phase at start of exon (0, 1, 2, or None)
    
    Returns:
        tuple: (fiveprime_intron, first_exon_edge, last_exon_edge, 
                threeprime_intron, first_edge_size, last_edge_size, valid)
    """
    
    exon_length = end - start
    
    # Determine first edge size (phase-adjusted to end at codon boundary)
    if exon_length >= (exon_edge_bp * 4):
        # Use phase-specific edge sizes for first edge
        if exon_phase == 0:
            first_edge_size = 39  # 13 codons
        elif exon_phase == 1:
            first_edge_size = 38  # (1 + 38) % 3 = 0
        elif exon_phase == 2:
            first_edge_size = 37  # (2 + 37) % 3 = 0
        else:
            # Phase unknown, use default
            first_edge_size = exon_edge_bp
        
        # For last edge, we want it to start at a codon boundary (phase 0)
        # The cumulative phase after the entire exon determines this
        if exon_phase is not None:
            exon_cumulative_phase = (exon_phase + exon_length) % 3
            # We want last edge to start at phase 0
            # So we need to back up to the nearest codon boundary
            if exon_cumulative_phase == 0:
                last_edge_size = 39  # Ends at phase 0, so 39bp back also starts at phase 0
            elif exon_cumulative_phase == 1:
                last_edge_size = 37  # Back 37bp: (1 - 37) % 3 = 0
            elif exon_cumulative_phase == 2:
                last_edge_size = 38  # Back 38bp: (2 - 38) % 3 = 0
        else:
            last_edge_size = exon_edge_bp
    else:
        # Short exon: use 25% rule, rounded to multiple of 3
        first_edge_size = max(3, int(exon_length * 0.25 / 3) * 3)
        last_edge_size = first_edge_size
    
    # Genomic coordinates
    upstream_intron_start = start - flank_size
    upstream_intron_end = start
    downstream_intron_start = end
    downstream_intron_end = end + flank_size
    
    first_exon_start = start
    first_exon_end = start + first_edge_size
    last_exon_start = end - last_edge_size
    last_exon_end = end
    
    # Check validity
    if upstream_intron_start < 0:
        return None, None, None, None, 0, 0, False
    
    # Check for overlap between first and last edges
    if first_exon_end > last_exon_start:
        return None, None, None, None, 0, 0, False
    
    try:
        # Get scores
        upstream_intron = np.array(bw.values(chrom, upstream_intron_start, upstream_intron_end))
        downstream_intron = np.array(bw.values(chrom, downstream_intron_start, downstream_intron_end))
        first_exon_edge = np.array(bw.values(chrom, first_exon_start, first_exon_end))
        last_exon_edge = np.array(bw.values(chrom, last_exon_start, last_exon_end))
        
        # Handle NaN values
        upstream_intron = np.nan_to_num(upstream_intron, nan=0.0)
        downstream_intron = np.nan_to_num(downstream_intron, nan=0.0)
        first_exon_edge = np.nan_to_num(first_exon_edge, nan=0.0)
        last_exon_edge = np.nan_to_num(last_exon_edge, nan=0.0)
        
        # Check lengths
        if (len(upstream_intron) != flank_size or 
            len(downstream_intron) != flank_size or
            len(first_exon_edge) != first_edge_size or 
            len(last_exon_edge) != last_edge_size):
            return None, None, None, None, 0, 0, False
        
        # Strand correction
        if strand == '+':
            fiveprime_intron = upstream_intron
            threeprime_intron = downstream_intron
            first_edge = first_exon_edge
            last_edge = last_exon_edge
        else:
            # Minus strand: flip everything
            fiveprime_intron = downstream_intron[::-1]
            threeprime_intron = upstream_intron[::-1]
            first_edge = last_exon_edge[::-1]
            last_edge = first_exon_edge[::-1]
        
        return fiveprime_intron, first_edge, last_edge, threeprime_intron, first_edge_size, last_edge_size, True
        
    except (RuntimeError, Exception) as e:
        print(f"Error at {chrom}:{start}-{end}: {e}", file=sys.stderr)
        return None, None, None, None, 0, 0, False

def bootstrap_confidence_interval(data, n_bootstrap=1000, ci=95):
    """Calculate bootstrap confidence intervals for mean at each position."""
    n_exons, n_positions = data.shape
    lower_bounds = np.zeros(n_positions)
    upper_bounds = np.zeros(n_positions)
    
    alpha = (100 - ci) / 2
    
    for pos in range(n_positions):
        bootstrap_means = []
        for _ in range(n_bootstrap):
            sample = np.random.choice(data[:, pos], size=n_exons, replace=True)
            bootstrap_means.append(np.mean(sample))
        
        lower_bounds[pos] = np.percentile(bootstrap_means, alpha)
        upper_bounds[pos] = np.percentile(bootstrap_means, 100 - alpha)
    
    return lower_bounds, upper_bounds


def calculate_conservation_profile(bed_file, phastcons_file, output_prefix, 
                                   chrom=None, flank_size=50, exon_edge_bp=40, 
                                   min_length=50, phylop_file=None, split_by_phase=False):
    """
    Calculate conservation profiles including exonic and intronic regions.
    
    Args:
        bed_file: Path to BED12 file
        phastcons_file: Path to PhastCons bigWig file
        output_prefix: Prefix for output files
        chrom: Chromosome to analyze (None = all chromosomes)
        flank_size: Size of flanking intronic regions in nt
        exon_edge_bp: Size of exonic edge regions in bp
        min_length: Minimum exon length to include
        phylop_file: Optional path to PhyloP bigWig file
        split_by_phase: If True, create separate output files for each phase (0, 1, 2)
    
    Returns:
        (summary_df, per_exon_df): DataFrames with results
    """
    
    chrom_msg = chrom if chrom else "ALL"
    print(f"Parsing BED12 for internal exons (min length: {min_length} nt)...")
    exons_df = parse_bed12_to_internal_exons(bed_file, chrom_filter=chrom, min_length=min_length)
    print(f"Found {len(exons_df)} internal exons on {chrom_msg}")
    
    if len(exons_df) == 0:
        print("ERROR: No internal exons found. Check chromosome name and min_length filter.")
        sys.exit(1)
    
    # Exon length distribution
    print(f"\nExon length distribution:")
    print(f"  Min: {exons_df['length'].min()} nt")
    print(f"  Median: {exons_df['length'].median():.0f} nt")
    print(f"  Max: {exons_df['length'].max()} nt")
    
    # Open bigWig files
    bw_phastcons = None
    if phastcons_file:
        print(f"\nOpening PhastCons bigWig: {phastcons_file}")
        phastcons_file = fix_bigwig_path(phastcons_file)
        bw_phastcons = pyBigWig.open(phastcons_file)
    
    bw_phylop = None
    if phylop_file:
        print(f"Opening PhyloP bigWig: {phylop_file}")
        phylop_file = fix_bigwig_path(phylop_file)
        bw_phylop = pyBigWig.open(phylop_file)
    
    # Store profiles
    profiles_data = []
    valid_exons = []
    
    print("\nExtracting conservation scores...")
    for idx, row in exons_df.iterrows():
        if idx % 5000 == 0 and idx > 0:
            print(f"  Processed {idx}/{len(exons_df)} exons...")
        
        # Get exon phase if available
        exon_phase = None
        if 'phase' in row and pd.notna(row['phase']) and row['phase'] != 'NA':
            try:
                phase_str = str(row['phase']).split(',')[0]  # Take first if comma-separated
                exon_phase = int(phase_str)
            except (ValueError, TypeError):
                pass
        
        profile_dict = {}
        valid = False
        first_edge_size = 0
        last_edge_size = 0
        
        # Get PhastCons scores if provided
        if bw_phastcons:
            five_intron_pc, first_edge_pc, last_edge_pc, three_intron_pc, first_sz, last_sz, valid = \
                get_exon_and_flanking_conservation(
                    bw_phastcons, row['chrom'], row['start'], row['end'], row['strand'], 
                    flank_size, exon_edge_bp, exon_phase
                )
            
            if valid:
                profile_dict['fiveprime_intron_phastcons'] = five_intron_pc
                profile_dict['first_exon_edge_phastcons'] = first_edge_pc
                profile_dict['last_exon_edge_phastcons'] = last_edge_pc
                profile_dict['threeprime_intron_phastcons'] = three_intron_pc
                profile_dict['first_edge_size'] = first_sz
                profile_dict['last_edge_size'] = last_sz
                first_edge_size = first_sz
                last_edge_size = last_sz
        
        # Get PhyloP scores if provided
        if bw_phylop:
            five_intron_pp, first_edge_pp, last_edge_pp, three_intron_pp, first_sz_pp, last_sz_pp, valid_pp = \
                get_exon_and_flanking_conservation(
                    bw_phylop, row['chrom'], row['start'], row['end'], row['strand'], 
                    flank_size, exon_edge_bp, exon_phase
                )
            
            if valid_pp:
                profile_dict['fiveprime_intron_phylop'] = five_intron_pp
                profile_dict['first_exon_edge_phylop'] = first_edge_pp
                profile_dict['last_exon_edge_phylop'] = last_edge_pp
                profile_dict['threeprime_intron_phylop'] = three_intron_pp
                if first_edge_size == 0:  # Use PhyloP edge sizes if PhastCons wasn't available
                    profile_dict['first_edge_size'] = first_sz_pp
                    profile_dict['last_edge_size'] = last_sz_pp
                    first_edge_size = first_sz_pp
                    last_edge_size = last_sz_pp
                valid = valid or valid_pp
        
        # Only include if at least one metric was valid
        if not valid or 'first_edge_size' not in profile_dict:
            continue
        
        profiles_data.append(profile_dict)
        valid_exons.append(row.to_dict())
    
    if bw_phastcons:
        bw_phastcons.close()
    if bw_phylop:
        bw_phylop.close()
    
    print(f"\nSuccessfully analyzed {len(valid_exons)} exons")
    print(f"Excluded {len(exons_df) - len(valid_exons)} exons due to coordinate issues")
    
    if len(valid_exons) == 0:
        print("ERROR: No valid exons to analyze.")
        sys.exit(1)
    
    valid_exons_df = pd.DataFrame(valid_exons)
    
    # Separate by exon edge size
    first_edge_sizes = [p['first_edge_size'] for p in profiles_data]
    last_edge_sizes = [p['last_edge_size'] for p in profiles_data]
    
    # For overall summary, only use exons where both edges equal the standard size
    standard_edge_mask = np.array([(f == exon_edge_bp and l == exon_edge_bp) 
                                    for f, l in zip(first_edge_sizes, last_edge_sizes)])
    
    print(f"\nExon edge size distribution:")
    print(f"  Standard (both edges {exon_edge_bp}bp): {standard_edge_mask.sum()} exons")
    print(f"  Variable (25% rule or phase-adjusted): {(~standard_edge_mask).sum()} exons")
    
    # Process standard-size edges for summary statistics
    summary_df = None
    if standard_edge_mask.sum() > 0:
        print("\nCalculating summary statistics (standard-size exons only)...")
        
        standard_profiles = [profiles_data[i] for i in range(len(profiles_data)) 
                            if standard_edge_mask[i]]
        
        # Initialize summary dataframe with positions and regions
        positions = (
            list(range(-flank_size, 0)) +
            [f'E{i+1}' for i in range(exon_edge_bp)] +
            [f'E-{exon_edge_bp-i}' for i in range(exon_edge_bp)] +
            list(range(1, flank_size + 1))
        )
        
        regions = (
            ["5'_intron"] * flank_size +
            ["first_exon_edge"] * exon_edge_bp +
            ["last_exon_edge"] * exon_edge_bp +
            ["3'_intron"] * flank_size
        )
        
        summary_df = pd.DataFrame({
            'position': positions,
            'region': regions
        })
        
        # PhastCons statistics if available
        has_phastcons = 'fiveprime_intron_phastcons' in standard_profiles[0]
        if has_phastcons:
            print("Calculating PhastCons statistics...")
            five_intron_pc_arr = np.array([p['fiveprime_intron_phastcons'] for p in standard_profiles])
            first_edge_pc_arr = np.array([p['first_exon_edge_phastcons'] for p in standard_profiles])
            last_edge_pc_arr = np.array([p['last_exon_edge_phastcons'] for p in standard_profiles])
            three_intron_pc_arr = np.array([p['threeprime_intron_phastcons'] for p in standard_profiles])
            
            combined_phastcons = np.concatenate([
                five_intron_pc_arr, first_edge_pc_arr, 
                last_edge_pc_arr, three_intron_pc_arr
            ], axis=1)
            
            mean_pc = np.mean(combined_phastcons, axis=0)
            median_pc = np.median(combined_phastcons, axis=0)
            sd_pc = np.std(combined_phastcons, axis=0, ddof=1)
            
            summary_df['phastcons_mean'] = mean_pc
            summary_df['phastcons_median'] = median_pc
            summary_df['phastcons_sd'] = sd_pc
        
        # PhyloP statistics if available
        has_phylop = 'fiveprime_intron_phylop' in standard_profiles[0]
        if has_phylop:
            print("Calculating PhyloP statistics...")
            five_intron_pp_arr = np.array([p['fiveprime_intron_phylop'] for p in standard_profiles])
            first_edge_pp_arr = np.array([p['first_exon_edge_phylop'] for p in standard_profiles])
            last_edge_pp_arr = np.array([p['last_exon_edge_phylop'] for p in standard_profiles])
            three_intron_pp_arr = np.array([p['threeprime_intron_phylop'] for p in standard_profiles])
            
            combined_phylop = np.concatenate([
                five_intron_pp_arr, first_edge_pp_arr,
                last_edge_pp_arr, three_intron_pp_arr
            ], axis=1)
            
            mean_pp = np.mean(combined_phylop, axis=0)
            median_pp = np.median(combined_phylop, axis=0)
            sd_pp = np.std(combined_phylop, axis=0, ddof=1)
            
            summary_df['phylop_mean'] = mean_pp
            summary_df['phylop_median'] = median_pp
            summary_df['phylop_sd'] = sd_pp
        
        summary_file = f"{output_prefix}_summary.tsv"
        summary_df.to_csv(summary_file, sep='\t', index=False)
        print(f"\nSaved summary to: {summary_file}")
    
    # Per-exon data file
    print("Creating per-exon output file...")
    
    for idx, profile_data in enumerate(profiles_data):
        # PhastCons means (if available)
        if 'fiveprime_intron_phastcons' in profile_data:
            valid_exons_df.loc[idx, 'phastcons_mean_5prime_intron'] = np.mean(profile_data['fiveprime_intron_phastcons'])
            valid_exons_df.loc[idx, 'phastcons_mean_first_exon_edge'] = np.mean(profile_data['first_exon_edge_phastcons'])
            valid_exons_df.loc[idx, 'phastcons_mean_last_exon_edge'] = np.mean(profile_data['last_exon_edge_phastcons'])
            valid_exons_df.loc[idx, 'phastcons_mean_3prime_intron'] = np.mean(profile_data['threeprime_intron_phastcons'])
        
        # PhyloP means (if available)
        if 'fiveprime_intron_phylop' in profile_data:
            valid_exons_df.loc[idx, 'phylop_mean_5prime_intron'] = np.mean(profile_data['fiveprime_intron_phylop'])
            valid_exons_df.loc[idx, 'phylop_mean_first_exon_edge'] = np.mean(profile_data['first_exon_edge_phylop'])
            valid_exons_df.loc[idx, 'phylop_mean_last_exon_edge'] = np.mean(profile_data['last_exon_edge_phylop'])
            valid_exons_df.loc[idx, 'phylop_mean_3prime_intron'] = np.mean(profile_data['threeprime_intron_phylop'])
        
        valid_exons_df.loc[idx, 'first_exon_edge_size'] = profile_data['first_edge_size']
        valid_exons_df.loc[idx, 'last_exon_edge_size'] = profile_data['last_edge_size']
    
    # =========================================================================
    # FULL EXON CONSERVATION (COMPLETELY SEPARATE FROM EDGE CALCULATIONS)
    # =========================================================================
    print("\nCalculating full exon conservation scores (separate from edge analysis)...")
    
    # Re-open bigWig files for full exon calculation
    bw_phastcons_full = None
    bw_phylop_full = None
    
    if phastcons_file:
        bw_phastcons_full = pyBigWig.open(fix_bigwig_path(phastcons_file))
    if phylop_file:
        bw_phylop_full = pyBigWig.open(fix_bigwig_path(phylop_file))
    
    # Initialize columns
    valid_exons_df['mean_exon_phastcons'] = np.nan
    valid_exons_df['mean_exon_phylop'] = np.nan
    
    # Iterate through valid exons and calculate full exon conservation
    for idx in range(len(valid_exons_df)):
        if idx % 5000 == 0 and idx > 0:
            print(f"  Processed {idx}/{len(valid_exons_df)} exons for full exon conservation...")
        
        row = valid_exons_df.iloc[idx]
        chrom_val = row['chrom']
        start_val = int(row['start'])
        end_val = int(row['end'])
        
        if bw_phastcons_full:
            mean_pc = get_full_exon_conservation(bw_phastcons_full, chrom_val, start_val, end_val)
            valid_exons_df.loc[valid_exons_df.index[idx], 'mean_exon_phastcons'] = mean_pc
        
        if bw_phylop_full:
            mean_pp = get_full_exon_conservation(bw_phylop_full, chrom_val, start_val, end_val)
            valid_exons_df.loc[valid_exons_df.index[idx], 'mean_exon_phylop'] = mean_pp
    
    # Close bigWig files
    if bw_phastcons_full:
        bw_phastcons_full.close()
    if bw_phylop_full:
        bw_phylop_full.close()
    
    print(f"  Added mean_exon_phastcons and mean_exon_phylop to {len(valid_exons_df)} exons")
    # =========================================================================
    # END FULL EXON CONSERVATION BLOCK
    # =========================================================================
    
    # =========================================================================
    # CLEAN UP PER-EXON DATAFRAME - REMOVE UNUSED/REDUNDANT COLUMNS
    # =========================================================================
    columns_to_remove = [
        'name',
        'all_phases',
        'phase_consistent', 
        'all_cumulative_phases',
        'cumulative_phase_consistent',
        'exon_num',
        'total_exons',
        'feature_types',
        'tx_names',
        'input_gene_name',
        'cassette_length',
        'skip_detected'
    ]
    
    # Only drop columns that actually exist in the dataframe
    cols_to_drop = [col for col in columns_to_remove if col in valid_exons_df.columns]
    if cols_to_drop:
        valid_exons_df = valid_exons_df.drop(columns=cols_to_drop)
        print(f"\nRemoved {len(cols_to_drop)} redundant columns: {', '.join(cols_to_drop)}")
    
    # Reorder columns to place mean_exon_phastcons and mean_exon_phylop after cassette_id
    current_cols = list(valid_exons_df.columns)
    
    # Find position of cassette_id (if it exists)
    if 'cassette_id' in current_cols:
        cassette_idx = current_cols.index('cassette_id')
        
        # Remove mean_exon columns from their current position
        cols_to_move = ['mean_exon_phastcons', 'mean_exon_phylop']
        for col in cols_to_move:
            if col in current_cols:
                current_cols.remove(col)
        
        # Find new cassette_id position after removal
        cassette_idx = current_cols.index('cassette_id')
        
        # Insert mean_exon columns right after cassette_id
        for i, col in enumerate(cols_to_move):
            if col in valid_exons_df.columns:
                current_cols.insert(cassette_idx + 1 + i, col)
        
        # Reorder dataframe
        valid_exons_df = valid_exons_df[current_cols]
        print(f"Reordered columns: mean_exon_phastcons and mean_exon_phylop now after cassette_id")

    # =========================================================================
    # END CLEANUP BLOCK
    # =========================================================================



    per_exon_file = f"{output_prefix}_per_exon.tsv"
    valid_exons_df.to_csv(per_exon_file, sep='\t', index=False)
    print(f"Saved per-exon data to: {per_exon_file}")
    
    # Phase-specific output
    if split_by_phase:
        print("\n" + "="*60)
        print("CREATING PHASE-SPECIFIC OUTPUT FILES")
        print("="*60)
        
        # Check if phase column exists
        if 'phase' not in valid_exons_df.columns:
            print("WARNING: No 'phase' column found, skipping phase split")
        else:
            # Parse phase values (handle comma-separated and NA)
            def get_primary_phase(phase_val):
                """Extract primary phase from phase column (handles '0,1' format)"""
                if pd.isna(phase_val) or phase_val == 'NA':
                    return None
                phase_str = str(phase_val)
                # Take first phase if comma-separated
                if ',' in phase_str:
                    phase_str = phase_str.split(',')[0]
                try:
                    return int(phase_str)
                except (ValueError, TypeError):
                    return None
            
            valid_exons_df['primary_phase'] = valid_exons_df['phase'].apply(get_primary_phase)
            
            # Split by phase
            for phase in [0, 1, 2]:
                phase_df = valid_exons_df[valid_exons_df['primary_phase'] == phase].copy()
                
                if len(phase_df) == 0:
                    print(f"\nPhase {phase}: No exons found")
                    continue
                
                # Remove the helper column
                phase_df = phase_df.drop('primary_phase', axis=1)
                
                phase_file = f"{output_prefix}_per_exon_phase{phase}.tsv"
                phase_df.to_csv(phase_file, sep='\t', index=False)
                print(f"\nPhase {phase}:")
                print(f"  Exons: {len(phase_df)}")
                print(f"  Saved to: {phase_file}")
                
                # Calculate phase-specific summary if we have standard-size exons
                # For each phase, determine the expected edge sizes
                if phase == 0:
                    expected_first_edge = 39
                    expected_last_edge = 39
                elif phase == 1:
                    expected_first_edge = 38
                    expected_last_edge = 37
                elif phase == 2:
                    expected_first_edge = 37
                    expected_last_edge = 38
                
                phase_indices = phase_df.index.tolist()
                # Filter for exons with the expected edge sizes for this phase
                # Need to check the actual edge sizes in the profiles_data
                phase_profiles_filtered = []
                for i in phase_indices:
                    if i < len(profiles_data):
                        prof = profiles_data[i]
                        if (prof.get('first_edge_size') == expected_first_edge and 
                            prof.get('last_edge_size') == expected_last_edge):
                            phase_profiles_filtered.append(prof)
                
                if len(phase_profiles_filtered) > 0:
                    print(f"  Calculating summary statistics for phase {phase} ({len(phase_profiles_filtered)} phase-standard exons)...")
                    
                    phase_profiles = phase_profiles_filtered
                    
                    # Check actual edge sizes present
                    actual_first_sizes = [p.get('first_edge_size', 0) for p in phase_profiles]
                    actual_last_sizes = [p.get('last_edge_size', 0) for p in phase_profiles]
                    unique_first = set(actual_first_sizes)
                    unique_last = set(actual_last_sizes)
                    
                    print(f"    Actual first edge sizes: {unique_first}")
                    print(f"    Actual last edge sizes: {unique_last}")
                    
                    # Use the maximum edge size for alignment (pad shorter ones with NaN)
                    max_first_edge = max(unique_first)
                    max_last_edge = max(unique_last)
                    
                    # Create phase-specific position labels using max sizes
                    phase_positions = (
                        list(range(-flank_size, 0)) +
                        [f'E{i+1}' for i in range(max_first_edge)] +
                        [f'E-{max_last_edge-i}' for i in range(max_last_edge)] +
                        list(range(1, flank_size + 1))
                    )
                    
                    phase_regions = (
                        ["5'_intron"] * flank_size +
                        ["first_exon_edge"] * max_first_edge +
                        ["last_exon_edge"] * max_last_edge +
                        ["3'_intron"] * flank_size
                    )
                    
                    # Initialize phase summary with positions and regions
                    phase_summary_df = pd.DataFrame({
                        'position': phase_positions,
                        'region': phase_regions
                    })
                    
                    # PhastCons statistics if available
                    has_phastcons = 'fiveprime_intron_phastcons' in phase_profiles[0]
                    if has_phastcons:
                        # Collect arrays with padding
                        five_intron_pc_list = []
                        first_edge_pc_list = []
                        last_edge_pc_list = []
                        three_intron_pc_list = []
                        
                        for p in phase_profiles:
                            five_intron_pc_list.append(p['fiveprime_intron_phastcons'])
                            
                            # Pad first edge to max_first_edge - always check actual length
                            first_edge = np.array(p['first_exon_edge_phastcons'])
                            actual_first_len = len(first_edge)
                            if actual_first_len < max_first_edge:
                                # Pad at the end (so E1 aligns)
                                padding = np.full(max_first_edge - actual_first_len, np.nan)
                                first_edge = np.concatenate([first_edge, padding])
                            elif actual_first_len > max_first_edge:
                                # Truncate if somehow longer
                                first_edge = first_edge[:max_first_edge]
                            first_edge_pc_list.append(first_edge)
                            
                            # Pad last edge to max_last_edge - always check actual length
                            last_edge = np.array(p['last_exon_edge_phastcons'])
                            actual_last_len = len(last_edge)
                            if actual_last_len < max_last_edge:
                                # Pad at the beginning (so E-1 aligns)
                                padding = np.full(max_last_edge - actual_last_len, np.nan)
                                last_edge = np.concatenate([padding, last_edge])
                            elif actual_last_len > max_last_edge:
                                # Truncate if somehow longer
                                last_edge = last_edge[-max_last_edge:]
                            last_edge_pc_list.append(last_edge)
                            
                            three_intron_pc_list.append(p['threeprime_intron_phastcons'])
                        
                        # Verify all arrays are the correct length
                        first_edge_lens = [len(x) for x in first_edge_pc_list]
                        last_edge_lens = [len(x) for x in last_edge_pc_list]
                        if len(set(first_edge_lens)) > 1 or len(set(last_edge_lens)) > 1:
                            print(f"    ERROR: Padding failed!")
                            print(f"      First edge lengths after padding: {set(first_edge_lens)}")
                            print(f"      Last edge lengths after padding: {set(last_edge_lens)}")
                            print(f"      Expected: first={max_first_edge}, last={max_last_edge}")
                            continue
                        
                        five_intron_pc_arr = np.array(five_intron_pc_list)
                        first_edge_pc_arr = np.array(first_edge_pc_list)
                        last_edge_pc_arr = np.array(last_edge_pc_list)
                        three_intron_pc_arr = np.array(three_intron_pc_list)
                        
                        combined_phastcons = np.concatenate([
                            five_intron_pc_arr, 
                            first_edge_pc_arr, 
                            last_edge_pc_arr, 
                            three_intron_pc_arr
                        ], axis=1)
                        
                        mean_pc = np.nanmean(combined_phastcons, axis=0)
                        median_pc = np.nanmedian(combined_phastcons, axis=0)
                        sd_pc = np.nanstd(combined_phastcons, axis=0, ddof=1)
                        
                        phase_summary_df['phastcons_mean'] = mean_pc
                        phase_summary_df['phastcons_median'] = median_pc
                        phase_summary_df['phastcons_sd'] = sd_pc
                    
                    # PhyloP statistics if available
                    has_phylop = 'fiveprime_intron_phylop' in phase_profiles[0]
                    if has_phylop:
                        # Collect arrays with padding
                        five_intron_pp_list = []
                        first_edge_pp_list = []
                        last_edge_pp_list = []
                        three_intron_pp_list = []
                        
                        for p in phase_profiles:
                            five_intron_pp_list.append(p['fiveprime_intron_phylop'])
                            
                            # Pad first edge - always check actual length
                            first_edge = np.array(p['first_exon_edge_phylop'])
                            actual_first_len = len(first_edge)
                            if actual_first_len < max_first_edge:
                                padding = np.full(max_first_edge - actual_first_len, np.nan)
                                first_edge = np.concatenate([first_edge, padding])
                            elif actual_first_len > max_first_edge:
                                # Truncate if somehow longer
                                first_edge = first_edge[:max_first_edge]
                            first_edge_pp_list.append(first_edge)
                            
                            # Pad last edge - always check actual length  
                            last_edge = np.array(p['last_exon_edge_phylop'])
                            actual_last_len = len(last_edge)
                            if actual_last_len < max_last_edge:
                                padding = np.full(max_last_edge - actual_last_len, np.nan)
                                last_edge = np.concatenate([padding, last_edge])
                            elif actual_last_len > max_last_edge:
                                # Truncate if somehow longer
                                last_edge = last_edge[-max_last_edge:]
                            last_edge_pp_list.append(last_edge)
                            
                            three_intron_pp_list.append(p['threeprime_intron_phylop'])
                        
                        # Verify all arrays are the correct length
                        first_edge_lens = [len(x) for x in first_edge_pp_list]
                        last_edge_lens = [len(x) for x in last_edge_pp_list]
                        if len(set(first_edge_lens)) > 1 or len(set(last_edge_lens)) > 1:
                            print(f"    ERROR: Padding failed!")
                            print(f"      First edge lengths after padding: {set(first_edge_lens)}")
                            print(f"      Last edge lengths after padding: {set(last_edge_lens)}")
                            print(f"      Expected: first={max_first_edge}, last={max_last_edge}")
                            continue
                        
                        five_intron_pp_arr = np.array(five_intron_pp_list)
                        first_edge_pp_arr = np.array(first_edge_pp_list)
                        last_edge_pp_arr = np.array(last_edge_pp_list)
                        three_intron_pp_arr = np.array(three_intron_pp_list)
                        
                        combined_phylop = np.concatenate([
                            five_intron_pp_arr,
                            first_edge_pp_arr,
                            last_edge_pp_arr,
                            three_intron_pp_arr
                        ], axis=1)
                        
                        mean_pp = np.nanmean(combined_phylop, axis=0)
                        median_pp = np.nanmedian(combined_phylop, axis=0)
                        sd_pp = np.nanstd(combined_phylop, axis=0, ddof=1)
                        
                        phase_summary_df['phylop_mean'] = mean_pp
                        phase_summary_df['phylop_median'] = median_pp
                        phase_summary_df['phylop_sd'] = sd_pp
                    
                    phase_summary_file = f"{output_prefix}_summary_phase{phase}.tsv"
                    phase_summary_df.to_csv(phase_summary_file, sep='\t', index=False)
                    print(f"  Saved phase {phase} summary to: {phase_summary_file}")
            
            # Also save exons with NA or ambiguous phase
            na_phase_df = valid_exons_df[valid_exons_df['primary_phase'].isna()].copy()
            if len(na_phase_df) > 0:
                na_phase_df = na_phase_df.drop('primary_phase', axis=1)
                na_file = f"{output_prefix}_per_exon_phase_NA.tsv"
                na_phase_df.to_csv(na_file, sep='\t', index=False)
                print(f"\nPhase NA/ambiguous:")
                print(f"  Exons: {len(na_phase_df)}")
                print(f"  Saved to: {na_file}")
            
            # Clean up helper column from original df
            valid_exons_df = valid_exons_df.drop('primary_phase', axis=1)
    
    print("\n=== ANALYSIS COMPLETE ===")
    print(f"Total exons analyzed: {len(valid_exons)}")
    metrics = []
    if phastcons_file:
        metrics.append("PhastCons")
    if phylop_file:
        metrics.append("PhyloP")
    print(f"Metrics: {' + '.join(metrics)}")
    
    return summary_df, valid_exons_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Calculate PhastCons and PhyloP conservation around AND within internal exons',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # PhastCons only
  python phastcons_phylop_exon_intron.py \\
      -b exons.bed -w hg38.phastCons100way.bw -o output/chr1
  
  # PhastCons + PhyloP with phase splitting
  python phastcons_phylop_exon_intron.py \\
      -b exons.bed \\
      -w hg38.phastCons100way.bw \\
      -p hg38.phyloP100way.bw \\
      -o output/chr1 -c chr1 \\
      --split-by-phase
        """
    )
    parser.add_argument('-b', '--bed', required=True,
                       help='Input BED12 file with exon coordinates')
    parser.add_argument('-w', '--phastcons', default=None,
                       help='PhastCons bigWig file (e.g., hg38.phastCons100way.bw)')
    parser.add_argument('-p', '--phylop', default=None,
                       help='PhyloP bigWig file (e.g., hg38.phyloP100way.bw)')
    parser.add_argument('-o', '--output', required=True,
                       help='Output prefix for result files')
    parser.add_argument('-c', '--chrom', default=None,
                       help='Chromosome to analyze (default: ALL)')
    parser.add_argument('-f', '--flank', type=int, default=50,
                       help='Flanking intron size in nt (default: 50)')
    parser.add_argument('-e', '--exon-edge', type=int, default=39,
                       help='Exon edge size in bp (default: 39, divisible by 3 for phase alignment)')
    parser.add_argument('-m', '--min-length', type=int, default=50,
                       help='Minimum exon length in nt (default: 50)')
    parser.add_argument('--split-by-phase', action='store_true',
                       help='Create separate output files for each phase (0, 1, 2)')
    
    args = parser.parse_args()
    
    # Validate that at least one conservation metric is provided
    if not args.phastcons and not args.phylop:
        parser.error("At least one of --phastcons (-w) or --phylop (-p) is required")
    
    print("="*60)
    print("Conservation Analysis (PhastCons + PhyloP)")
    print("="*60)
    print(f"BED file: {args.bed}")
    if args.phastcons:
        print(f"PhastCons bigWig: {args.phastcons}")
    if args.phylop:
        print(f"PhyloP bigWig: {args.phylop}")
    print(f"Output prefix: {args.output}")
    print(f"Chromosome: {args.chrom if args.chrom else 'ALL'}")
    print(f"Intron flank: {args.flank} nt")
    print(f"Exon edge: {args.exon_edge} bp")
    print(f"Min exon length: {args.min_length} nt")
    print(f"Split by phase: {args.split_by_phase}")
    print("="*60 + "\n")
    
    summary_df, per_exon_df = calculate_conservation_profile(
        args.bed, 
        args.phastcons, 
        args.output,
        chrom=args.chrom,
        flank_size=args.flank,
        exon_edge_bp=args.exon_edge,
        min_length=args.min_length,
        phylop_file=args.phylop,
        split_by_phase=args.split_by_phase
    )
    
    if summary_df is not None:
        print("\n=== SAMPLE OUTPUT ===")
        print(summary_df.head(10))