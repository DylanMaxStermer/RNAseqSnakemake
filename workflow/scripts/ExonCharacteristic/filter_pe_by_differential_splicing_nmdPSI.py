#!/usr/bin/env python3
"""
Filter poison exons to keep only those with differentially spliced flanking introns.

Works with minimal format (coding_junction, upstream_intron, downstream_intron)
or full BED format with additional columns.

Filtering criteria:
1. Flanking intron must be in a cluster with significant differential splicing (p.adjust <= cutoff)
2. Flanking intron must have deltapsi > cutoff (increased inclusion upon NMD inhibition)
3. Optional: Flanking intron must have NMD_enriched >= threshold (minimum PSI in NMD-inhibited condition)

Author: Dylan Stermer
Date: December 2024
Updated: January 2026 - Support minimal format, added NMD_enriched filtering
"""

import pandas as pd
import argparse
import sys

def parse_intron_coord(intron_str):
    """
    Parse intron coordinate string to chr:start-end format.
    
    Input: 'chr1:32894879-32895305' or 'chr1:32894879:32895305'
    Output: 'chr1:32894879-32895305'
    """
    if not intron_str or intron_str == 'NA':
        return None
    
    # Handle both formats
    parts = str(intron_str).replace(':', '-').split('-')
    if len(parts) >= 3:
        chrom = parts[0]
        start = parts[1]
        end = parts[2]
        return f"{chrom}:{start}-{end}"
    return None

def load_differential_splicing_results(results_files, p_cutoff=0.05, deltapsi_cutoff=0, nmd_psi_cutoff=None):
    """
    Load all differential splicing results and extract significant introns.
    
    Args:
        results_files: List of paths to Results.tsv.gz files
        p_cutoff: Adjusted p-value cutoff (default: 0.05)
        deltapsi_cutoff: Delta PSI cutoff (default: 0)
        nmd_psi_cutoff: Minimum NMD_enriched value (default: None, no filtering)
    
    Returns:
        Set of significant intron coordinates (chr:start-end format)
    """
    
    significant_introns = set()
    
    for results_file in results_files:
        print(f"Loading {results_file}...")
        
        # Read gzipped file
        if str(results_file).endswith('.gz'):
            df = pd.read_csv(results_file, sep='\t', compression='gzip')
        else:
            df = pd.read_csv(results_file, sep='\t')
        
        print(f"  Total junctions: {len(df)}")
        
        # Check for required columns
        if 'NMD_enriched' not in df.columns and nmd_psi_cutoff is not None:
            print(f"  WARNING: NMD_enriched column not found, skipping NMD PSI filter for this file")
            nmd_filter_available = False
        else:
            nmd_filter_available = True
        
        # Build filter mask
        # Base filters: p.adjust <= cutoff AND deltapsi > cutoff
        sig_mask = (df['p.adjust'] <= p_cutoff) & (df['deltapsi'] > deltapsi_cutoff)
        
        # Optional NMD_enriched filter
        if nmd_psi_cutoff is not None and nmd_filter_available:
            sig_mask = sig_mask & (df['NMD_enriched'] >= nmd_psi_cutoff)
        
        sig_df = df[sig_mask]
        
        # Print filtering stats
        filter_desc = f"p.adjust <= {p_cutoff}, deltapsi > {deltapsi_cutoff}"
        if nmd_psi_cutoff is not None and nmd_filter_available:
            filter_desc += f", NMD_enriched >= {nmd_psi_cutoff}"
        
        print(f"  Significant ({filter_desc}): {len(sig_df)}")
        
        # Extract intron coordinates from Intron_coord column
        for intron_coord in sig_df['Intron_coord']:
            if pd.notna(intron_coord):
                significant_introns.add(intron_coord)
    
    print(f"\nTotal unique significant intron coordinates: {len(significant_introns)}")
    return significant_introns

def filter_poison_exons(pe_file, significant_introns, output_file):
    """
    Filter poison exon file to keep only exons with differentially spliced flanking introns.
    
    Works with minimal format (coding_junction, upstream_intron, downstream_intron)
    or full BED format.
    
    Args:
        pe_file: Path to poison exon file
        significant_introns: Set of significant intron coordinates
        output_file: Path to output filtered file
    """
    
    print(f"\nReading {pe_file}...")
    
    # Try reading with pandas (handles headers automatically)
    df = pd.read_csv(pe_file, sep='\t', comment='#')
    
    print(f"Columns: {list(df.columns)}")
    print(f"Total poison exons: {len(df)}")
    
    # Check for required columns
    required_cols = ['upstream_intron', 'downstream_intron']
    missing_cols = [col for col in required_cols if col not in df.columns]
    
    if missing_cols:
        print(f"ERROR: Missing required columns: {missing_cols}")
        sys.exit(1)
    
    # Parse and filter
    kept_mask = []
    filter_reasons = {'upstream': 0, 'downstream': 0, 'both': 0, 'neither': 0}
    
    for idx, row in df.iterrows():
        upstream_intron = parse_intron_coord(row['upstream_intron'])
        downstream_intron = parse_intron_coord(row['downstream_intron'])
        
        upstream_sig = upstream_intron in significant_introns
        downstream_sig = downstream_intron in significant_introns
        
        # Keep exon if either flanking intron is significantly differential
        if upstream_sig or downstream_sig:
            kept_mask.append(True)
            if upstream_sig and downstream_sig:
                filter_reasons['both'] += 1
            elif upstream_sig:
                filter_reasons['upstream'] += 1
            else:
                filter_reasons['downstream'] += 1
        else:
            kept_mask.append(False)
            filter_reasons['neither'] += 1
    
    # Create filtered dataframe
    filtered_df = df[kept_mask].copy()
    
    print(f"\nFiltered poison exons: {len(filtered_df)}")
    print(f"Removed: {len(df) - len(filtered_df)}")
    print(f"Retention rate: {len(filtered_df)/len(df)*100:.1f}%")
    
    print(f"\nFilter breakdown:")
    print(f"  Kept (upstream intron significant): {filter_reasons['upstream']}")
    print(f"  Kept (downstream intron significant): {filter_reasons['downstream']}")
    print(f"  Kept (both introns significant): {filter_reasons['both']}")
    print(f"  Removed (neither intron significant): {filter_reasons['neither']}")
    
    # Write output in same format as input
    filtered_df.to_csv(output_file, sep='\t', index=False)
    
    print(f"\nWrote filtered exons to {output_file}")
    
    # Print summary statistics (only if columns are available)
    print(f"\n{'='*70}")
    print("SUMMARY STATISTICS")
    print(f"{'='*70}")
    
    # Basic stats always available
    print(f"\nFiltered exons: {len(filtered_df)}")
    print(f"Retention rate: {len(filtered_df)/len(df)*100:.1f}%")
    
    # Optional stats if columns exist
    if 'coding_junction' in filtered_df.columns:
        print(f"\nCoding junctions:")
        print(f"  Total unique: {filtered_df['coding_junction'].nunique()}")
    
    if 'name' in filtered_df.columns:
        print(f"\nExons per gene:")
        gene_counts = filtered_df['name'].value_counts()
        print(f"  Total genes: {len(gene_counts)}")
        print(f"  Mean exons per gene: {gene_counts.mean():.1f}")
        print(f"  Median exons per gene: {gene_counts.median():.0f}")
        print(f"  Max exons in one gene: {gene_counts.max()}")
        
        print(f"\nTop 10 genes by exon count:")
        for gene, count in gene_counts.head(10).items():
            print(f"  {gene}: {count}")
    
    if 'length' in filtered_df.columns:
        print(f"\nLength distribution:")
        print(f"  Mean: {filtered_df['length'].mean():.1f} nt")
        print(f"  Median: {filtered_df['length'].median():.0f} nt")
        print(f"  Range: {filtered_df['length'].min()}-{filtered_df['length'].max()} nt")
    
    if 'strand' in filtered_df.columns:
        print(f"\nStrand distribution:")
        for strand in ['+', '-']:
            count = (filtered_df['strand'] == strand).sum()
            if count > 0:
                pct = count / len(filtered_df) * 100
                print(f"  {strand}: {count} ({pct:.1f}%)")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Filter poison exons by differential splicing of flanking introns',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:

Basic filtering (p-value and deltapsi only):
  python filter_pe_by_differential_splicing.py \\
      -i input/poison_exons.tsv \\
      -r results/differential_splicing_tidy/*/Results.tsv.gz \\
      -o output/poison_exons_filtered.tsv

With NMD PSI threshold (require minimum inclusion in NMD-inhibited condition):
  python filter_pe_by_differential_splicing.py \\
      -i input/poison_exons.tsv \\
      -r results/differential_splicing_tidy/*/Results.tsv.gz \\
      -o output/poison_exons_filtered.tsv \\
      --nmd-psi 0.05

Stringent filtering:
  python filter_pe_by_differential_splicing.py \\
      -i input/poison_exons.tsv \\
      -r results/differential_splicing_tidy/*/Results.tsv.gz \\
      -o output/poison_exons_filtered_stringent.tsv \\
      -p 0.01 \\
      -d 0.05 \\
      --nmd-psi 0.10

This will keep poison exons where upstream_intron or downstream_intron
appears in a cluster with:
  - Adjusted p-value <= 0.05 (or specified -p value)
  - Delta PSI > 0 (or specified -d value; increased inclusion upon NMD inhibition)
  - NMD_enriched >= threshold (if --nmd-psi specified; minimum PSI in NMD condition)
        """
    )
    
    parser.add_argument('-i', '--input', required=True,
                       help='Input poison exon file (TSV with upstream_intron and downstream_intron columns)')
    parser.add_argument('-r', '--results', nargs='+', required=True,
                       help='Differential splicing Results.tsv.gz file(s)')
    parser.add_argument('-o', '--output', required=True,
                       help='Output filtered file (same format as input)')
    parser.add_argument('-p', '--p-cutoff', type=float, default=0.05,
                       help='Adjusted p-value cutoff (default: 0.05)')
    parser.add_argument('-d', '--deltapsi-cutoff', type=float, default=0,
                       help='Delta PSI cutoff (default: 0)')
    parser.add_argument('--nmd-psi', type=float, default=None,
                       help='Minimum NMD_enriched value (PSI in NMD-inhibited condition). '
                            'If not specified, no NMD PSI filtering is applied. '
                            'Recommended values: 0.02-0.10 depending on stringency.')
    
    args = parser.parse_args()
    
    print("="*70)
    print("FILTER POISON EXONS BY DIFFERENTIAL SPLICING")
    print("="*70)
    print(f"Input file: {args.input}")
    print(f"Number of Results files: {len(args.results)}")
    for r in args.results:
        print(f"  - {r}")
    print(f"Output: {args.output}")
    print(f"\nFilter cutoffs:")
    print(f"  Adjusted p-value: <= {args.p_cutoff}")
    print(f"  Delta PSI: > {args.deltapsi_cutoff}")
    if args.nmd_psi is not None:
        print(f"  NMD_enriched (PSI in NMD condition): >= {args.nmd_psi}")
    else:
        print(f"  NMD_enriched: No filtering (not specified)")
    print("="*70 + "\n")
    
    # Load differential splicing results
    significant_introns = load_differential_splicing_results(
        args.results, 
        args.p_cutoff, 
        args.deltapsi_cutoff,
        args.nmd_psi
    )
    
    # Filter poison exons
    filter_poison_exons(args.input, significant_introns, args.output)
    
    print("\n" + "="*70)
    print("FILTERING COMPLETE")
    print("="*70)