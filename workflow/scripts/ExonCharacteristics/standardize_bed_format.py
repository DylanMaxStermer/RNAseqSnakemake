

#!/usr/bin/env python3
"""
Standardize diverse exon annotation files into uniform format for PhastCons analysis.

Reorders columns to put standard format first, keeps all other columns after.

Standard format columns (always first):
  chrom, start, end, name, strand, length, gene_name, gene_type, tx_type, exon_position

Author: Dylan Stermer
Date: January 2026
"""

import pandas as pd
import argparse
import sys

def remove_unwanted_columns(df):
    """
    Remove columns not needed for downstream PhastCons analysis.
    Only called AFTER standardization uses them to populate name column.
    
    Args:
        df: DataFrame with standardized columns
    
    Returns:
        DataFrame with unwanted columns removed
    """
    columns_to_remove = [
        'phase_consistent',
        'cumulative_phase_consistent', 
        'exon_num',
        'total_exons',
        'constitutive',
        'pe_name',
        'cassette_id',
        'input_gene_name',
        'n_transcripts'
    ]
    
    # Only remove columns that actually exist
    existing_to_remove = [col for col in columns_to_remove if col in df.columns]
    
    if existing_to_remove:
        print(f"  Removing {len(existing_to_remove)} unwanted columns: {existing_to_remove}")
        df = df.drop(columns=existing_to_remove)
    else:
        print(f"  No unwanted columns found to remove")
    
    return df

def detect_file_format(df):
    """
    Detect which format the input file is in based on column names.
    
    Returns:
        format_type: str identifying the format
    """
    cols = set(df.columns)
    
    # Check for specific column patterns
    if 'mazin_tag' in cols:
        return 'ortholog_match'
    elif 'pe_name' in cols and 'position_relative_pe' in cols:
        return 'flanking_match'
    elif 'pe_name' not in cols and 'position_relative_pe' not in cols and 'mazin_tag' not in cols and 'cassette_id' not in cols:
        # Check if it's poison exon annotated (has phase but simpler structure)
        if 'feature_types' in cols and df['feature_types'].iloc[0] == 'poisonExon':
            return 'poison_annotated'
    elif 'cassette_id' in cols and 'input_gene_name' in cols:
        return 'cassette_match'
    
    return 'unknown'

def get_exon_position_from_columns(row, df):
    """
    Extract exon_position from wherever it exists in the row.
    Searches common column names.
    """
    # Try direct column
    if 'exon_position' in df.columns and pd.notna(row.get('exon_position')):
        return row['exon_position']
    
    # Try position_relative_pe (for flanking exons)
    if 'position_relative_pe' in df.columns and pd.notna(row.get('position_relative_pe')):
        return row['position_relative_pe']
    
    # Default
    return 'NA'

def get_gene_name_from_columns(row, df):
    """
    Extract gene_name from wherever it exists.
    Priority: gtf_gene_name > gene_name > input_gene_name
    """
    # Try GTF gene name first (most authoritative)
    if 'gtf_gene_name' in df.columns and pd.notna(row.get('gtf_gene_name')) and row.get('gtf_gene_name') != 'NA':
        return row['gtf_gene_name']
    
    # Try standard gene_name
    if 'gene_name' in df.columns and pd.notna(row.get('gene_name')) and row.get('gene_name') != 'NA':
        return row['gene_name']
    
    # Try input_gene_name (for cassettes)
    if 'input_gene_name' in df.columns and pd.notna(row.get('input_gene_name')) and row.get('input_gene_name') != 'NA':
        return row['input_gene_name']
    
    return 'NA'

def get_name_from_columns(row, df):
    """
    Extract a name/identifier from wherever it exists.
    Priority order depends on file type.
    """
    # Try specific identifiers first
    if 'mazin_tag' in df.columns and pd.notna(row.get('mazin_tag')):
        return row['mazin_tag']
    
    if 'pe_name' in df.columns and pd.notna(row.get('pe_name')):
        return row['pe_name']
    
    if 'cassette_id' in df.columns and pd.notna(row.get('cassette_id')):
        return row['cassette_id']
    
    if 'name' in df.columns and pd.notna(row.get('name')):
        return row['name']
    
    # Fallback: create from coordinates
    return f"{row['chrom']}:{row['start']}-{row['end']}_{row.get('strand', '.')}"

def parse_column_names_from_comments(input_file):
    """
    Parse column names from comment header lines.
    Looks for lines like: # Columns: chrom, start, end, ...
    
    Returns:
        list of column names or None if not found
    """
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('# Columns:'):
                # Extract column names after "# Columns:"
                cols_str = line.replace('# Columns:', '').strip()
                # Split by comma and strip whitespace
                col_names = [col.strip() for col in cols_str.split(',')]
                return col_names
            # Stop at first non-comment line
            if not line.startswith('#'):
                break
    return None

def deduplicate_dataframe(df, format_type):
    """
    Remove duplicate rows from dataframe.
    
    For cassette_match format: if rows are identical except for skip_detected,
    keep the row where skip_detected=True.
    
    For all formats: remove exact duplicate rows.
    
    Args:
        df: DataFrame to deduplicate
        format_type: detected file format
    
    Returns:
        Deduplicated DataFrame
    """
    n_before = len(df)
    
    if format_type == 'cassette_match' and 'skip_detected' in df.columns:
        # For cassette files: prioritize skip_detected=True
        
        # Normalize skip_detected to boolean for reliable sorting
        df['_skip_bool'] = df['skip_detected'].apply(
            lambda x: True if x == True or x == 'True' or x == 'true' or x == 1 else False
        )
        
        # Define key columns (everything except skip_detected and our temp column)
        key_columns = [c for c in df.columns if c not in ['skip_detected', '_skip_bool']]
        
        # Sort so True comes first (ascending=False puts True before False)
        df = df.sort_values('_skip_bool', ascending=False)
        
        # Drop duplicates based on key columns, keeping first (which is True if exists)
        df = df.drop_duplicates(subset=key_columns, keep='first')
        
        # Remove temp column
        df = df.drop(columns=['_skip_bool'])
        
        # Also remove any exact full duplicates
        df = df.drop_duplicates()
        
    else:
        # For other formats: simple exact duplicate removal
        df = df.drop_duplicates()
    
    n_after = len(df)
    
    if n_before != n_after:
        print(f"  Deduplication: removed {n_before - n_after} duplicate rows ({n_before} -> {n_after})")
    else:
        print(f"  Deduplication: no duplicates found ({n_after} rows)")
    
    return df

def standardize_file(input_file, output_file):
    """
    Standardize input file to uniform format.
    
    Reorders columns to put standard format first:
      chrom, start, end, name, strand, length, gene_name, gene_type, tx_type, exon_position
    
    Keeps all other columns after these.
    Removes duplicate rows.
    """
    print(f"Reading input file: {input_file}")
    
    # Try to parse column names from comments
    col_names = parse_column_names_from_comments(input_file)
    
    if col_names:
        print(f"  Parsed column names from comments: {len(col_names)} columns")
        # Read file with parsed column names
        df = pd.read_csv(input_file, sep='\t', comment='#', header=None, names=col_names)
    else:
        print(f"  No column names in comments, reading with header row")
        # Read file, handling comment lines
        df = pd.read_csv(input_file, sep='\t', comment='#')
    
    print(f"  Rows: {len(df)}")
    print(f"  Columns: {len(df.columns)}")
    
    # Detect format
    format_type = detect_file_format(df)
    print(f"  Detected format: {format_type}")
    
    # Verify required columns exist
    required = ['chrom', 'start', 'end', 'strand']
    missing = [col for col in required if col not in df.columns]
    if missing:
        print(f"ERROR: Missing required columns: {missing}")
        sys.exit(1)
    
    # Create standardized columns
    print("\nStandardizing columns...")
    
    standardized_df = pd.DataFrame()
    
    # Standard columns (always first, in this order)
    standardized_df['chrom'] = df['chrom']
    standardized_df['start'] = df['start']
    standardized_df['end'] = df['end']
    
    # Name - extract from various possible columns
    standardized_df['name'] = df.apply(lambda row: get_name_from_columns(row, df), axis=1)
    
    standardized_df['strand'] = df['strand']
    
    # Length - calculate if not present
    if 'length' in df.columns:
        standardized_df['length'] = df['length']
    elif 'cassette_length' in df.columns:
        standardized_df['length'] = df['cassette_length']
    else:
        standardized_df['length'] = df['end'] - df['start']
    
    # Gene name - extract from various columns
    standardized_df['gene_name'] = df.apply(lambda row: get_gene_name_from_columns(row, df), axis=1)
    
    # Gene type
    if 'gene_type' in df.columns:
        standardized_df['gene_type'] = df['gene_type']
    else:
        standardized_df['gene_type'] = 'NA'
    
    # Transcript type
    if 'tx_type' in df.columns:
        standardized_df['tx_type'] = df['tx_type']
    else:
        standardized_df['tx_type'] = 'NA'
    
    # Exon position
    standardized_df['exon_position'] = df.apply(lambda row: get_exon_position_from_columns(row, df), axis=1)
    
    # Add all remaining columns (that aren't already in standardized_df)
    standard_cols = set(standardized_df.columns)
    remaining_cols = [col for col in df.columns if col not in standard_cols]
    
    for col in remaining_cols:
        standardized_df[col] = df[col]
    
    print(f"\nStandardized format:")
    print(f"  First 10 columns (standard): {list(standardized_df.columns[:10])}")
    print(f"  Total columns: {len(standardized_df.columns)}")
    print(f"  Additional columns kept: {len(remaining_cols)}")
    
    # Remove unwanted columns AFTER name has been populated
    print("\nRemoving unwanted columns...")
    standardized_df = remove_unwanted_columns(standardized_df)
    print(f"  Columns after removal: {len(standardized_df.columns)}")
    
    # Deduplicate
    print("\nRemoving duplicates...")
    standardized_df = deduplicate_dataframe(standardized_df, format_type)
    
    # Write output
    print(f"\nWriting to: {output_file}")
    
    with open(output_file, 'w') as f:
        # Write header comments
        f.write("# Standardized exon annotation file\n")
        f.write("# Standard columns (first 10): chrom, start, end, name, strand, length, ")
        f.write("gene_name, gene_type, tx_type, exon_position\n")
        f.write("# Coordinates: 0-based half-open (BED format)\n")
        f.write(f"# Original format: {format_type}\n")
        f.write(f"# Total columns: {len(standardized_df.columns)}\n")
        
        # Write data
        standardized_df.to_csv(f, sep='\t', index=False)
    
    print(f"✓ Standardization complete!")
    
    # Show sample
    print(f"\nSample output (first 3 rows, standard columns only):")
    print(standardized_df[['chrom', 'start', 'end', 'name', 'strand', 'length', 
                          'gene_name', 'gene_type', 'tx_type', 'exon_position']].head(3))
    
    return standardized_df


def batch_standardize(input_files, output_dir):
    """
    Standardize multiple files at once.
    
    Args:
        input_files: List of input file paths
        output_dir: Output directory for standardized files
    """
    import os
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")
    
    print("="*70)
    print("BATCH STANDARDIZATION")
    print("="*70)
    print(f"Input files: {len(input_files)}")
    print(f"Output directory: {output_dir}")
    print("="*70 + "\n")
    
    for i, input_file in enumerate(input_files, 1):
        print(f"\n{'='*70}")
        print(f"Processing file {i}/{len(input_files)}")
        print(f"{'='*70}")
        
        # Generate output filename
        basename = os.path.basename(input_file)
        # Remove .bed extension if present
        if basename.endswith('.bed'):
            basename = basename[:-4]
        output_file = os.path.join(output_dir, f"{basename}_standardized.tsv")
        
        try:
            standardize_file(input_file, output_file)
        except Exception as e:
            print(f"ERROR processing {input_file}: {e}")
            continue
    
    print("\n" + "="*70)
    print("BATCH PROCESSING COMPLETE")
    print("="*70)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Standardize exon annotation files for PhastCons analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Standardizes diverse exon annotation file formats into uniform format.

Standard format (first 10 columns):
  chrom, start, end, name, strand, length, gene_name, gene_type, tx_type, exon_position

All other columns are preserved and placed after the standard columns.
Duplicate rows are removed. For cassette files, rows with skip_detected=True
are prioritized when duplicates differ only in that column.

Supported input formats:
  - Ortholog GTF matches (mazin_tag)
  - Flanking exon matches (pe_name, position_relative_pe)
  - Poison exon annotated (feature_types=poisonExon)
  - Cassette exon matches (cassette_id)

Examples:
  # Single file
  python standardize_exon_files.py \\
      -i input.bed \\
      -o standardized_output.tsv
  
  # Batch processing
  python standardize_exon_files.py \\
      -i file1.bed file2.bed file3.bed \\
      -d output_directory/
        """
    )
    
    parser.add_argument('-i', '--input', nargs='+', required=True,
                       help='Input file(s) to standardize')
    parser.add_argument('-o', '--output', 
                       help='Output file (for single input) or use -d for batch')
    parser.add_argument('-d', '--output-dir',
                       help='Output directory for batch processing')
    
    args = parser.parse_args()
    
    # Validate arguments
    if len(args.input) == 1 and args.output:
        # Single file mode
        standardize_file(args.input[0], args.output)
    elif args.output_dir:
        # Batch mode
        batch_standardize(args.input, args.output_dir)
    else:
        print("ERROR: Either provide single input with -o, or multiple inputs with -d")
        sys.exit(1)