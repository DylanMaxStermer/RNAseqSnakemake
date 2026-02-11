#!/usr/bin/env python3
"""
Extract splice site coordinates and sequences from exon BED files.

5' splice site: 3 bases in exon + 6 bases in intron (donor site)
3' splice site: 20 bases in intron + 3 bases in exon (acceptor site)

Strand-aware extraction (genomic coordinates are NOT pre-flipped).

Used code from:
https://github.com/kepbod/maxentpy?tab=readme-ov-file

Dec 14, 2025
Dylan Stermer
"""

import pandas as pd
import argparse
import sys
import pyfaidx

def parse_exon_bed(bed_file):
    """
    Parse exon BED file - handles both old format and standardized TSV format.
    
    Standardized format columns (first 10):
    chrom, start, end, name, strand, length, gene_name, gene_type, tx_type, exon_position
    """
    
    exons = []
    header_skipped = False
    column_indices = None
    
    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            
            # Skip header row (detect by checking if 'chrom' or 'start' is in first two fields)
            if not header_skipped and (fields[0] == 'chrom' or fields[1] == 'start'):
                # Use header to build column index map
                column_indices = {name: idx for idx, name in enumerate(fields)}
                header_skipped = True
                continue
            
            # If no header was found, assume old 18-column format
            if column_indices is None:
                if len(fields) < 18:
                    print(f"WARNING: Line has only {len(fields)} columns, expected 18", file=sys.stderr)
                    continue
                
                exons.append({
                    'chrom': fields[0],
                    'start': int(fields[1]),
                    'end': int(fields[2]),
                    'name': fields[3],
                    'score': int(fields[4]),
                    'strand': fields[5],
                    'length': int(fields[6]),
                    'exon_position': fields[10],
                    'gene_name': fields[15],
                    'gene_type': fields[16],
                    'tx_type': fields[17]
                })
            else:
                # Use standardized format with header-based indexing
                try:
                    exons.append({
                        'chrom': fields[column_indices['chrom']],
                        'start': int(fields[column_indices['start']]),
                        'end': int(fields[column_indices['end']]),
                        'name': fields[column_indices['name']],
                        'score': int(fields[column_indices.get('original_score', column_indices.get('score', -1))]) if column_indices.get('original_score') or column_indices.get('score') else 0,
                        'strand': fields[column_indices['strand']],
                        'length': int(fields[column_indices['length']]),
                        'exon_position': fields[column_indices['exon_position']],
                        'gene_name': fields[column_indices['gene_name']],
                        'gene_type': fields[column_indices['gene_type']],
                        'tx_type': fields[column_indices['tx_type']]
                    })
                except (KeyError, ValueError) as e:
                    print(f"WARNING: Failed to parse line: {e}", file=sys.stderr)
                    continue
    
    return pd.DataFrame(exons)

def extract_splice_sites(exons_df, fasta_file, donor_exon_bp=3, donor_intron_bp=6, 
                         acceptor_intron_bp=20, acceptor_exon_bp=3):
    """
    Extract splice site coordinates and sequences from exons.
    
    5' splice site (donor):
        - Last 3bp of exon + first 6bp of intron (genomic downstream)
        - Only for exons that are NOT last/single
        
    3' splice site (acceptor):
        - Last 20bp of intron + first 3bp of exon (genomic upstream)
        - Only for exons that are NOT first/single
    
    Args:
        exons_df: DataFrame with exon coordinates
        fasta_file: Path to reference genome FASTA file
        donor_exon_bp: Bases from exon for donor (default: 3)
        donor_intron_bp: Bases from intron for donor (default: 6)
        acceptor_intron_bp: Bases from intron for acceptor (default: 20)
        acceptor_exon_bp: Bases from exon for acceptor (default: 3)
    
    Returns:
        (donor_sites_df, acceptor_sites_df): DataFrames with splice site coordinates and sequences
    """
    
    # Load reference genome
    print(f"Loading reference genome: {fasta_file}")
    genome = pyfaidx.Fasta(fasta_file)
    
    donor_sites = []
    acceptor_sites = []
    
    for idx, exon in exons_df.iterrows():
        if idx % 10000 == 0 and idx > 0:
            print(f"  Processed {idx}/{len(exons_df)} exons...")
        
        chrom = exon['chrom']
        start = exon['start']
        end = exon['end']
        strand = exon['strand']
        position = exon['exon_position']
        
        # Check if chromosome exists in genome
        if chrom not in genome:
            continue
        
        # 5' splice site (donor) - skip last and single exons
        if position not in ['last', 'single']:
            if strand == '+':
                # Plus strand: donor is at exon END (3' end of exon)
                # Get last 3bp of exon + first 6bp of downstream intron
                donor_start = end - donor_exon_bp
                donor_end = end + donor_intron_bp
                
                # Extract sequence (pyfaidx uses 0-based coordinates)
                seq = str(genome[chrom][donor_start:donor_end])
                
            else:
                # Minus strand: donor is at exon START (genomically)
                # which is the 3' end of the exon in transcript orientation
                # Get 6bp upstream intron + first 3bp of exon
                donor_start = start - donor_intron_bp
                donor_end = start + donor_exon_bp
                
                # Extract sequence and reverse complement
                seq = str(genome[chrom][donor_start:donor_end].reverse.complement)
            
            donor_sites.append({
                'chrom': chrom,
                'start': donor_start,
                'end': donor_end,
                'name': f"{exon['name']}_donor",
                'score': exon['score'],
                'strand': strand,
                'sequence': seq,
                'site_type': 'donor',
                'exon_bp': donor_exon_bp,
                'intron_bp': donor_intron_bp,
                'gene_name': exon['gene_name'],
                'gene_type': exon['gene_type'],
                'tx_type': exon['tx_type']
            })
        
        # 3' splice site (acceptor) - skip first and single exons
        if position not in ['first', 'single']:
            if strand == '+':
                # Plus strand: acceptor is at exon START (5' end of exon)
                # Get last 20bp of upstream intron + first 3bp of exon
                acceptor_start = start - acceptor_intron_bp
                acceptor_end = start + acceptor_exon_bp
                
                # Extract sequence
                seq = str(genome[chrom][acceptor_start:acceptor_end])
                
            else:
                # Minus strand: acceptor is at exon END (genomically)
                # which is the 5' end of the exon in transcript orientation
                # Get last 3bp of exon + 20bp downstream intron
                acceptor_start = end - acceptor_exon_bp
                acceptor_end = end + acceptor_intron_bp
                
                # Extract sequence and reverse complement
                seq = str(genome[chrom][acceptor_start:acceptor_end].reverse.complement)
            
            acceptor_sites.append({
                'chrom': chrom,
                'start': acceptor_start,
                'end': acceptor_end,
                'name': f"{exon['name']}_acceptor",
                'score': exon['score'],
                'strand': strand,
                'sequence': seq,
                'site_type': 'acceptor',
                'exon_bp': acceptor_exon_bp,
                'intron_bp': acceptor_intron_bp,
                'gene_name': exon['gene_name'],
                'gene_type': exon['gene_type'],
                'tx_type': exon['tx_type']
            })
    
    donor_df = pd.DataFrame(donor_sites)
    acceptor_df = pd.DataFrame(acceptor_sites)
    
    return donor_df, acceptor_df

def write_splice_sites(sites_df, output_file, site_type):
    """Write splice sites to BED format with sequences."""
    
    if len(sites_df) == 0:
        print(f"WARNING: No {site_type} sites to write")
        return
    
    with open(output_file, 'w') as f:
        # Write header
        f.write(f"# Splice site type: {site_type}\n")
        f.write("# Columns: chrom, start, end, name, score, strand, sequence, site_type, exon_bp, intron_bp, gene_name, gene_type, tx_type\n")
        
        if site_type == 'donor':
            f.write("# Donor (5' splice site): Last 3bp of exon + first 6bp of intron\n")
            f.write("#   Plus strand: exon_end - 3 to exon_end + 6\n")
            f.write("#   Minus strand: exon_start - 6 to exon_start + 3\n")
            f.write("#   Sequence is in transcript orientation (reverse complemented for minus strand)\n")
        else:
            f.write("# Acceptor (3' splice site): Last 20bp of intron + first 3bp of exon\n")
            f.write("#   Plus strand: exon_start - 20 to exon_start + 3\n")
            f.write("#   Minus strand: exon_end - 3 to exon_end + 20\n")
            f.write("#   Sequence is in transcript orientation (reverse complemented for minus strand)\n")
        
        # Write data
        for _, row in sites_df.iterrows():
            f.write(f"{row['chrom']}\t{row['start']}\t{row['end']}\t{row['name']}\t")
            f.write(f"{row['score']}\t{row['strand']}\t{row['sequence']}\t{row['site_type']}\t")
            f.write(f"{row['exon_bp']}\t{row['intron_bp']}\t")
            f.write(f"{row['gene_name']}\t{row['gene_type']}\t{row['tx_type']}\n")
    
    print(f"Wrote {len(sites_df)} {site_type} sites to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Extract splice site coordinates from exon BED files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Splice Site Definitions:

5' splice site (DONOR):
  - Last 3bp of exon + first 6bp of intron
  - Plus strand: exon_end - 3 to exon_end + 6
  - Minus strand: exon_start - 6 to exon_start + 3
  - Canonical sequence: AG|GTRAGT (exon|intron)
  - Extracted for all exons EXCEPT last and single

3' splice site (ACCEPTOR):
  - Last 20bp of intron + first 3bp of exon
  - Plus strand: exon_start - 20 to exon_start + 3
  - Minus strand: exon_end - 3 to exon_end + 20
  - Canonical sequence: (Y)nNYAG|G (intron|exon)
  - Extracted for all exons EXCEPT first and single

Output Format:
  BED6 + additional columns:
    - site_type: "donor" or "acceptor"
    - exon_bp: number of bases from exon
    - intron_bp: number of bases from intron
    - gene_name: gene name
    - gene_type: gene biotype
    - tx_type: transcript biotype

Example:
  python extract_splice_sites.py \\
      -i cummulative_constitutive_exons.bed \\
      -o output/constitutive_splice_sites
        """
    )
    
    parser.add_argument('-i', '--input', required=True,
                       help='Input exon BED file')
    parser.add_argument('-f', '--fasta', required=True,
                       help='Reference genome FASTA file (e.g., hg38.fa)')
    parser.add_argument('-o', '--output', required=True,
                       help='Output prefix for splice site files')
    
    args = parser.parse_args()
    
    print("="*70)
    print("SPLICE SITE EXTRACTION")
    print("="*70)
    print(f"Input: {args.input}")
    print(f"Reference genome: {args.fasta}")
    print(f"Output prefix: {args.output}")
    print(f"Donor site: 3bp exon + 6bp intron")
    print(f"Acceptor site: 20bp intron + 3bp exon")
    print("="*70 + "\n")
    
    # Parse exons
    print("Parsing exon BED file...")
    exons_df = parse_exon_bed(args.input)
    print(f"Found {len(exons_df)} exons")
    
    # Count exon positions
    print("\nExon position distribution:")
    for pos in ['first', 'internal', 'last', 'single']:
        count = (exons_df['exon_position'] == pos).sum()
        if count > 0:
            print(f"  {pos.capitalize()}: {count}")
    
    # Extract splice sites
    print("\nExtracting splice sites...")
    donor_df, acceptor_df = extract_splice_sites(
        exons_df,
        args.fasta
    )
    
    print(f"Extracted {len(donor_df)} donor sites (5' splice sites)")
    print(f"Extracted {len(acceptor_df)} acceptor sites (3' splice sites)")
    
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
    
    donor_file = f"{args.output}_donor_sites.bed"
    write_splice_sites(donor_df, donor_file, "donor")
    
    print()
    
    acceptor_file = f"{args.output}_acceptor_sites.bed"
    write_splice_sites(acceptor_df, acceptor_file, "acceptor")
    
    print("\n" + "="*70)
    print("EXTRACTION COMPLETE")
    print("="*70)
    print(f"\nOutput files:")
    print(f"  1. {donor_file}")
    print(f"  2. {acceptor_file}")
    
    # Summary statistics
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"Total exons processed: {len(exons_df)}")
    print(f"  First exons (no acceptor): {(exons_df['exon_position'] == 'first').sum()}")
    print(f"  Internal exons (both sites): {(exons_df['exon_position'] == 'internal').sum()}")
    print(f"  Last exons (no donor): {(exons_df['exon_position'] == 'last').sum()}")
    print(f"  Single exons (no sites): {(exons_df['exon_position'] == 'single').sum()}")
    print(f"\nSplice sites extracted:")
    print(f"  Donor sites: {len(donor_df)}")
    print(f"  Acceptor sites: {len(acceptor_df)}")