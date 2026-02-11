#!/usr/bin/env python3
"""
Calculate MaxEntScan scores for splice sites using maxentpy.

Reads splice site BED files with sequences and adds MaxEnt scores.
Uses fast matrix-based scoring for better performance.

Dec 14, 2025
Dylan Stermer
"""

import pandas as pd
import argparse
import sys
from maxentpy import maxent
from maxentpy.maxent import load_matrix5, load_matrix3
from maxentpy.maxent_fast import score5 as score5_fast, score3 as score3_fast

# Load matrices once at startup for fast scoring
print("Loading MaxEnt scoring matrices...")
matrix5 = load_matrix5()
matrix3 = load_matrix3()
print("Matrices loaded successfully!")

def parse_splice_site_bed(bed_file):
    """
    Parse splice site BED file with sequences.
    
    Expected columns:
    0: chrom
    1: start
    2: end
    3: name
    4: score
    5: strand
    6: sequence
    7: site_type
    8: exon_bp
    9: intron_bp
    10: gene_name
    11: gene_type
    12: tx_type
    """
    
    sites = []
    
    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            
            if len(fields) < 13:
                print(f"WARNING: Line has only {len(fields)} columns, expected 13", file=sys.stderr)
                continue
            
            sites.append({
                'chrom': fields[0],
                'start': int(fields[1]),
                'end': int(fields[2]),
                'name': fields[3],
                'score': int(fields[4]),
                'strand': fields[5],
                'sequence': fields[6],
                'site_type': fields[7],
                'exon_bp': int(fields[8]),
                'intron_bp': int(fields[9]),
                'gene_name': fields[10],
                'gene_type': fields[11],
                'tx_type': fields[12]
            })
    
    return pd.DataFrame(sites)

def calculate_maxent_scores(sites_df, site_type):
    """
    Calculate MaxEntScan scores for splice sites using fast matrix-based scoring.
    
    Args:
        sites_df: DataFrame with splice site information
        site_type: 'donor' or 'acceptor'
    
    Returns:
        DataFrame with added maxent_score column
    """
    
    print(f"Calculating MaxEnt scores for {len(sites_df)} {site_type} sites...")
    
    maxent_scores = []
    failed_count = 0
    
    for idx, row in sites_df.iterrows():
        if idx % 10000 == 0 and idx > 0:
            print(f"  Processed {idx}/{len(sites_df)} sites...")
        
        sequence = row['sequence'].upper()
        
        try:
            if site_type == 'donor':
                # Donor: 3 exon + 6 intron bases (9 total)
                if len(sequence) != 9:
                    print(f"WARNING: Donor sequence length is {len(sequence)}, expected 9: {sequence}", file=sys.stderr)
                    maxent_scores.append(None)
                    failed_count += 1
                    continue
                # Use fast matrix-based scoring
                score = score5_fast(sequence, matrix=matrix5)
            else:  # acceptor
                # Acceptor: 20 intron + 3 exon bases (23 total)
                if len(sequence) != 23:
                    print(f"WARNING: Acceptor sequence length is {len(sequence)}, expected 23: {sequence}", file=sys.stderr)
                    maxent_scores.append(None)
                    failed_count += 1
                    continue
                # Use fast matrix-based scoring
                score = score3_fast(sequence, matrix=matrix3)
            
            maxent_scores.append(score)
            
        except Exception as e:
            print(f"ERROR calculating score for {row['name']}: {e}", file=sys.stderr)
            print(f"  Sequence: {sequence}", file=sys.stderr)
            maxent_scores.append(None)
            failed_count += 1
    
    sites_df['maxent_score'] = maxent_scores
    
    print(f"Successfully scored {len(sites_df) - failed_count}/{len(sites_df)} sites")
    if failed_count > 0:
        print(f"Failed to score {failed_count} sites")
    
    return sites_df

def write_scored_sites(sites_df, output_file, site_type):
    """Write splice sites with MaxEnt scores to file."""
    
    with open(output_file, 'w') as f:
        # Write header
        f.write(f"# Splice site type: {site_type}\n")
        f.write("# Columns: chrom, start, end, name, score, strand, sequence, maxent_score, site_type, exon_bp, intron_bp, gene_name, gene_type, tx_type\n")
        
        if site_type == 'donor':
            f.write("# Donor (5' splice site): 3 exon + 6 intron bases\n")
            f.write("# MaxEntScan score calculated with maxent.score5()\n")
        else:
            f.write("# Acceptor (3' splice site): 20 intron + 3 exon bases\n")
            f.write("# MaxEntScan score calculated with maxent.score3()\n")
        
        # Write data
        for _, row in sites_df.iterrows():
            f.write(f"{row['chrom']}\t{row['start']}\t{row['end']}\t{row['name']}\t")
            f.write(f"{row['score']}\t{row['strand']}\t{row['sequence']}\t")
            
            # Handle None values
            if pd.isna(row['maxent_score']):
                f.write("NA\t")
            else:
                f.write(f"{row['maxent_score']:.6f}\t")
            
            f.write(f"{row['site_type']}\t{row['exon_bp']}\t{row['intron_bp']}\t")
            f.write(f"{row['gene_name']}\t{row['gene_type']}\t{row['tx_type']}\n")
    
    print(f"Wrote {len(sites_df)} scored {site_type} sites to {output_file}")

def print_score_summary(sites_df, site_type):
    """Print summary statistics for MaxEnt scores."""
    
    valid_scores = sites_df['maxent_score'].dropna()
    
    if len(valid_scores) == 0:
        print(f"No valid scores for {site_type} sites")
        return
    
    print(f"\n{site_type.upper()} SITE MaxEnt SCORE SUMMARY:")
    print(f"  Count: {len(valid_scores)}")
    print(f"  Mean: {valid_scores.mean():.4f}")
    print(f"  Median: {valid_scores.median():.4f}")
    print(f"  Std Dev: {valid_scores.std():.4f}")
    print(f"  Min: {valid_scores.min():.4f}")
    print(f"  Max: {valid_scores.max():.4f}")
    print(f"  Q1 (25%): {valid_scores.quantile(0.25):.4f}")
    print(f"  Q3 (75%): {valid_scores.quantile(0.75):.4f}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Calculate MaxEntScan scores for splice sites',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
MaxEntScan Scoring:
  
  Uses FAST matrix-based scoring for better performance:
    - score5_fast() for donor sites (up to 10x faster)
    - score3_fast() for acceptor sites (up to 13x faster)
  
  Donor sites (5' splice sites):
    - Requires 9bp sequence: 3 exon + 6 intron bases
    - Higher scores indicate stronger splice sites
  
  Acceptor sites (3' splice sites):
    - Requires 23bp sequence: 20 intron + 3 exon bases
    - Higher scores indicate stronger splice sites

Input Format:
  BED file with columns including sequence in column 7

Output Format:
  Same as input but with maxent_score added after sequence column

Example:
  python calculate_maxent_scores.py \\
      -d input/donor_sites.bed \\
      -a input/acceptor_sites.bed \\
      -o output/scored_splice_sites
        """
    )
    
    parser.add_argument('-d', '--donor', required=True,
                       help='Input donor sites BED file with sequences')
    parser.add_argument('-a', '--acceptor', required=True,
                       help='Input acceptor sites BED file with sequences')
    parser.add_argument('-o', '--output', required=True,
                       help='Output prefix for scored splice site files')
    
    args = parser.parse_args()
    
    print("="*70)
    print("MaxEntScan SPLICE SITE SCORING")
    print("="*70)
    print(f"Donor sites: {args.donor}")
    print(f"Acceptor sites: {args.acceptor}")
    print(f"Output prefix: {args.output}")
    print("="*70 + "\n")
    
    # Process donor sites
    print("Processing donor (5') splice sites...")
    donor_df = parse_splice_site_bed(args.donor)
    print(f"Loaded {len(donor_df)} donor sites")
    
    donor_df = calculate_maxent_scores(donor_df, 'donor')
    
    # Process acceptor sites
    print("\nProcessing acceptor (3') splice sites...")
    acceptor_df = parse_splice_site_bed(args.acceptor)
    print(f"Loaded {len(acceptor_df)} acceptor sites")
    
    acceptor_df = calculate_maxent_scores(acceptor_df, 'acceptor')
    
    # Create output directory if needed
    import os
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"\nCreated output directory: {output_dir}")
    
    # Write output files
    print("\n" + "="*70)
    print("WRITING OUTPUT FILES")
    print("="*70)
    
    donor_file = f"{args.output}_donor_sites_scored.bed"
    write_scored_sites(donor_df, donor_file, "donor")
    
    print()
    
    acceptor_file = f"{args.output}_acceptor_sites_scored.bed"
    write_scored_sites(acceptor_df, acceptor_file, "acceptor")
    
    # Print summaries
    print("\n" + "="*70)
    print("SCORE SUMMARIES")
    print("="*70)
    
    print_score_summary(donor_df, "donor")
    print_score_summary(acceptor_df, "acceptor")
    
    print("\n" + "="*70)
    print("SCORING COMPLETE")
    print("="*70)
    print(f"\nOutput files:")
    print(f"  1. {donor_file}")
    print(f"  2. {acceptor_file}")