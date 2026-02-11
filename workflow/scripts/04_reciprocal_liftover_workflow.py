#!/usr/bin/env python3
"""
Reciprocal liftOver workflow for identifying high-confidence orthologous exons.

Performs forward liftover, reverse liftover, calculates coordinate tolerance,
and classifies exons into 1:1 and 1:many ortholog relationships.

Usage:
    python 04_reciprocal_liftover_workflow.py \
        --input exons.tsv \
        --forward-chain hg38ToMm39.chain.gz \
        --reverse-chain mm39ToHg38.chain.gz \
        --output-dir ./output/ \
        --source-assembly hg38 \
        --target-assembly mm39
"""

import argparse
import csv
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from collections import defaultdict

# ============================================================================
# PHASE 1: SETUP AND BASE FUNCTIONS - Status: [X]
# ============================================================================

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Reciprocal liftOver workflow for orthologous exon mapping.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
    python 04_reciprocal_liftover_workflow.py \\
        --input /path/to/exons.tsv \\
        --forward-chain /path/to/hg38ToMm39.chain.gz \\
        --reverse-chain /path/to/mm39ToHg38.chain.gz \\
        --output-dir /path/to/output/ \\
        --source-assembly hg38 \\
        --target-assembly mm39
        """
    )
    parser.add_argument('--input', '-i', required=True,
                        help='Input TSV file with exon coordinates')
    parser.add_argument('--forward-chain', '-f', required=True,
                        help='Path to forward chain file (.chain or .chain.gz)')
    parser.add_argument('--reverse-chain', '-r', required=True,
                        help='Path to reverse chain file (.chain or .chain.gz)')
    parser.add_argument('--output-dir', '-o', required=True,
                        help='Output directory for results')
    parser.add_argument('--source-assembly', '-s', required=True,
                        help='Source assembly name (e.g., hg38)')
    parser.add_argument('--target-assembly', '-t', required=True,
                        help='Target assembly name (e.g., mm39)')
    parser.add_argument('--conda-env', default=None,
                        help='Path to conda environment with liftOver (optional)')
    parser.add_argument('--min-match', type=float, default=0.7,
                        help='Minimum ratio of bases that must remap (default: 0.7)')
    return parser.parse_args()


def parse_input_tsv(input_path):
    """
    Parse input TSV file, extracting exon coordinates.

    Args:
        input_path: Path to input TSV file

    Returns:
        List of dicts with keys: chrom, start, end, strand, gene_name, original_name
    """
    exons = []

    with open(input_path, 'r') as f:
        # Find header line (skip comments)
        header = None
        for line in f:
            if line.startswith('#'):
                continue
            header = line.strip().split('\t')
            break

        if header is None:
            raise ValueError("No header line found in input file")

        # Find column indices
        try:
            chrom_idx = header.index('chrom')
            start_idx = header.index('start')
            end_idx = header.index('end')
            strand_idx = header.index('strand')
            gene_name_idx = header.index('gene_name')
        except ValueError as e:
            raise ValueError(f"Missing required column: {e}")

        # Optional name column
        name_idx = header.index('name') if 'name' in header else None

        # Parse data rows
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) <= max(chrom_idx, start_idx, end_idx, strand_idx, gene_name_idx):
                continue

            exon = {
                'chrom': row[chrom_idx],
                'start': int(row[start_idx]),
                'end': int(row[end_idx]),
                'strand': row[strand_idx],
                'gene_name': row[gene_name_idx],
                'original_name': row[name_idx] if name_idx is not None else None
            }
            exons.append(exon)

    return exons


def deduplicate_exons(exons):
    """
    Deduplicate exons by (chrom, start, end, strand).

    Args:
        exons: List of exon dicts

    Returns:
        List of unique exon dicts with added 'exon_id' field
    """
    seen = {}
    unique_exons = []

    for exon in exons:
        key = (exon['chrom'], exon['start'], exon['end'], exon['strand'])
        if key not in seen:
            seen[key] = True
            unique_exons.append(exon)

    # Add sequential exon IDs
    for i, exon in enumerate(unique_exons, start=1):
        exon['exon_id'] = f"exon_{i:06d}"

    return unique_exons


def write_bed_file(exons, bed_path):
    """
    Write exons to BED6 format for liftOver.

    Args:
        exons: List of exon dicts with exon_id or coordinate info
        bed_path: Path to output BED file
    """
    with open(bed_path, 'w') as f:
        for exon in exons:
            # Determine name field (exon_id or create from coordinates)
            if 'exon_id' in exon:
                name = exon['exon_id']
            else:
                # For reverse lift, create unique identifier from target coords
                name = f"{exon['target_chr']}:{exon['target_start']}-{exon['target_end']}"

            # BED6: chrom, start, end, name, score, strand
            chrom = exon.get('chrom') or exon.get('target_chr')
            start = exon.get('start') or exon.get('target_start')
            end = exon.get('end') or exon.get('target_end')
            strand = exon.get('strand') or exon.get('target_strand')

            f.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n")


def run_liftover(input_bed, chain_file, output_bed, unmapped_bed, min_match, conda_env=None, use_multiple=False):
    """
    Run UCSC liftOver tool.

    Args:
        input_bed: Path to input BED file
        chain_file: Path to chain file
        output_bed: Path to output lifted BED file
        unmapped_bed: Path to output unmapped BED file
        min_match: Minimum match ratio
        conda_env: Path to conda environment (optional)
        use_multiple: Whether to use -multiple flag (optional)

    Returns:
        True if successful, False otherwise
    """
    # Build command
    if conda_env:
        # Use conda run to execute in environment
        cmd = [
            'conda', 'run', '-p', conda_env,
            'liftOver',
            str(input_bed),
            str(chain_file),
            str(output_bed),
            str(unmapped_bed),
            f'-minMatch={min_match}'
        ]
    else:
        cmd = [
            'liftOver',
            str(input_bed),
            str(chain_file),
            str(output_bed),
            str(unmapped_bed),
            f'-minMatch={min_match}'
        ]

    # Add -multiple flag if requested
    if use_multiple:
        cmd.append('-multiple')

    print(f"Running: {' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"liftOver stderr: {result.stderr}", file=sys.stderr)
            return False
        return True
    except FileNotFoundError:
        print("Error: liftOver not found. Please specify --conda-env or ensure liftOver is in PATH.", file=sys.stderr)
        return False


def parse_lifted_bed(bed_path, allow_multiple=False):
    """
    Parse lifted BED file.

    Args:
        bed_path: Path to BED file
        allow_multiple: If True, handle multiple mappings per ID (returns list per key)

    Returns:
        If allow_multiple=False: Dict mapping name to coordinates
        If allow_multiple=True: Dict mapping name to list of coordinates
    """
    if not os.path.exists(bed_path):
        return {} if not allow_multiple else defaultdict(list)

    if allow_multiple:
        lifted = defaultdict(list)
        with open(bed_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 6:
                    name = parts[3]
                    lifted[name].append({
                        'chr': parts[0],
                        'start': int(parts[1]),
                        'end': int(parts[2]),
                        'strand': parts[5]
                    })
        return dict(lifted)
    else:
        lifted = {}
        with open(bed_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 6:
                    name = parts[3]
                    # Keep only first mapping if multiple exist
                    if name not in lifted:
                        lifted[name] = {
                            'chr': parts[0],
                            'start': int(parts[1]),
                            'end': int(parts[2]),
                            'strand': parts[5]
                        }
        return lifted


def parse_unmapped_bed(bed_path, min_match=0.7):
    """
    Parse unmapped BED file to extract failure reasons.

    Args:
        bed_path: Path to unmapped BED file
        min_match: Minimum match threshold used (for labeling partial matches)

    Returns:
        Dict mapping name to unmapped_reason
    """
    unmapped = {}
    current_reason = "Unknown"

    if not os.path.exists(bed_path):
        return unmapped

    with open(bed_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                # Reason line (e.g., "#Deleted in new")
                current_reason = line[1:].strip()
                # Replace "Partially deleted in new" with threshold-specific message
                if current_reason == "Partially deleted in new":
                    current_reason = f"match < {int(min_match*100)}%"
            elif line:
                parts = line.split('\t')
                if len(parts) >= 4:
                    name = parts[3]
                    unmapped[name] = current_reason

    return unmapped


# ============================================================================
# PHASE 2: CORE RECIPROCAL LIFTOVER FUNCTIONS - Status: [X]
# ============================================================================

def calculate_coordinate_tolerance(original_start, original_end, reciprocal_start, reciprocal_end):
    """
    Calculate the maximum coordinate shift after reciprocal lift.

    Args:
        original_start: Original start coordinate
        original_end: Original end coordinate
        reciprocal_start: Reciprocal-lifted start coordinate
        reciprocal_end: Reciprocal-lifted end coordinate

    Returns:
        Tuple of (tolerance_str, start_diff, end_diff)
        - tolerance_str: "0", "1", ..., "10", or ">10"
        - start_diff: Absolute difference in start coordinates
        - end_diff: Absolute difference in end coordinates
    """
    start_diff = abs(original_start - reciprocal_start)
    end_diff = abs(original_end - reciprocal_end)
    tolerance = max(start_diff, end_diff)

    if tolerance > 10:
        tolerance_str = ">10"
    else:
        tolerance_str = str(tolerance)

    return tolerance_str, start_diff, end_diff


def run_reciprocal_liftover(forward_lifted, reverse_chain, temp_dir, min_match, conda_env):
    """
    Takes forward-lifted exons and lifts them back to original assembly.

    Args:
        forward_lifted: Dict of exon_id -> list of target coordinates from forward lift
        reverse_chain: Path to reverse chain file
        temp_dir: Temporary directory for BED files
        min_match: Minimum match threshold
        conda_env: Conda environment path

    Returns:
        Dict mapping exon_id -> list of reverse-lifted coordinates
    """
    # Create BED file from forward-lifted coordinates
    reverse_input_bed = Path(temp_dir) / "forward_lifted.bed"
    reverse_output_bed = Path(temp_dir) / "reverse_lifted.bed"
    reverse_unmapped_bed = Path(temp_dir) / "reverse_unmapped.bed"

    # Convert forward_lifted to list format for BED writing
    # IMPORTANT: Don't include exon_id so write_bed_file() creates coordinate-based names
    # that match the lookup keys used in classify_ortholog_mapping()
    bed_entries = []
    for exon_id, target_coords_list in forward_lifted.items():
        for target_coords in target_coords_list:
            bed_entries.append({
                'target_chr': target_coords['chr'],
                'target_start': target_coords['start'],
                'target_end': target_coords['end'],
                'target_strand': target_coords['strand']
            })

    write_bed_file(bed_entries, reverse_input_bed)

    # Run liftOver with -multiple flag
    success = run_liftover(
        reverse_input_bed, reverse_chain, reverse_output_bed, reverse_unmapped_bed,
        min_match, conda_env, use_multiple=True
    )

    if not success:
        print("Warning: Reverse liftOver failed", file=sys.stderr)
        return {}

    # Parse reverse lift results
    reverse_lifted = parse_lifted_bed(reverse_output_bed, allow_multiple=True)

    return reverse_lifted


def classify_ortholog_mapping(unique_exons, forward_lifted, reverse_lifted):
    """
    Classify exons into reciprocal vs non-reciprocal, and 1:1 vs 1:many.

    Args:
        unique_exons: List of original exon dicts
        forward_lifted: Dict of exon_id -> list of target coordinates
        reverse_lifted: Dict of target_coord_key -> list of source coordinates

    Returns:
        Tuple of (reciprocal_hits, non_reciprocal, one_to_one, one_to_many, quality_metrics)
    """
    reciprocal_hits = []
    non_reciprocal = []
    one_to_one = []
    one_to_many = []
    quality_metrics = []

    # Create lookup dict for original exons
    exon_lookup = {exon['exon_id']: exon for exon in unique_exons}

    # Track mapping groups for 1:many
    mapping_group_counter = 1

    for exon in unique_exons:
        exon_id = exon['exon_id']
        original_chr = exon['chrom']
        original_start = exon['start']
        original_end = exon['end']
        original_strand = exon['strand']
        gene_name = exon['gene_name']

        # Check if this exon was forward-lifted
        if exon_id not in forward_lifted:
            # Failed forward lift
            quality_metrics.append({
                'exon_id': exon_id,
                'gene_name': gene_name,
                'forward_lift_success': False,
                'reverse_lift_success': False,
                'reciprocal_match': False,
                'tolerance': 'N/A',
                'mapping_type': 'failed',
                'source_coords': f"{original_chr}:{original_start}-{original_end}",
                'target_coords': 'N/A',
                'reciprocal_coords': 'N/A'
            })
            non_reciprocal.append(exon_id)
            continue

        # Get forward-lifted coordinates (can be multiple)
        target_coords_list = forward_lifted[exon_id]
        mapping_type = "1:1" if len(target_coords_list) == 1 else "1:many"

        # Check each forward-lifted target for reciprocal mapping
        has_reciprocal = False
        best_quality_metric = None  # Track best/first quality metric for this exon

        for target_coords in target_coords_list:
            target_chr = target_coords['chr']
            target_start = target_coords['start']
            target_end = target_coords['end']
            target_strand = target_coords['strand']

            # Create key for reverse lookup
            target_key = f"{target_chr}:{target_start}-{target_end}"

            # Check if this target was reverse-lifted
            if target_key in reverse_lifted:
                reverse_coords_list = reverse_lifted[target_key]

                # Check if any reverse-lifted coords match original
                chr_matched = False
                for reverse_coords in reverse_coords_list:
                    reciprocal_chr = reverse_coords['chr']
                    reciprocal_start = reverse_coords['start']
                    reciprocal_end = reverse_coords['end']

                    # Calculate tolerance
                    tolerance_str, start_diff, end_diff = calculate_coordinate_tolerance(
                        original_start, original_end, reciprocal_start, reciprocal_end
                    )

                    # Check if chromosomes match
                    if reciprocal_chr == original_chr:
                        chr_matched = True
                        has_reciprocal = True

                        # Create reciprocal hit record
                        hit_record = {
                            'exon_id': exon_id,
                            'gene_name': gene_name,
                            'source_chr': original_chr,
                            'source_start': original_start,
                            'source_end': original_end,
                            'strand': original_strand,
                            'target_chr': target_chr,
                            'target_start': target_start,
                            'target_end': target_end,
                            'reciprocal_chr': reciprocal_chr,
                            'reciprocal_start': reciprocal_start,
                            'reciprocal_end': reciprocal_end,
                            'tolerance': tolerance_str,
                            'start_diff': start_diff,
                            'end_diff': end_diff,
                            'mapping_type': mapping_type
                        }

                        reciprocal_hits.append(hit_record)

                        # Add to 1:1 or 1:many
                        ortholog_record = hit_record.copy()
                        if mapping_type == "1:1":
                            one_to_one.append(ortholog_record)
                        else:
                            ortholog_record['mapping_group'] = mapping_group_counter
                            one_to_many.append(ortholog_record)

                        # Store best quality metric (first reciprocal hit for this exon)
                        if best_quality_metric is None:
                            best_quality_metric = {
                                'exon_id': exon_id,
                                'gene_name': gene_name,
                                'forward_lift_success': True,
                                'reverse_lift_success': True,
                                'reciprocal_match': True,
                                'tolerance': tolerance_str,
                                'mapping_type': mapping_type,
                                'source_coords': f"{original_chr}:{original_start}-{original_end}",
                                'target_coords': f"{target_chr}:{target_start}-{target_end}",
                                'reciprocal_coords': f"{reciprocal_chr}:{reciprocal_start}-{reciprocal_end}"
                            }

                # If reverse lifted but chromosome didn't match, record that
                if not chr_matched and best_quality_metric is None:
                    best_quality_metric = {
                        'exon_id': exon_id,
                        'gene_name': gene_name,
                        'forward_lift_success': True,
                        'reverse_lift_success': True,
                        'reciprocal_match': False,
                        'tolerance': 'N/A',
                        'mapping_type': 'failed',
                        'source_coords': f"{original_chr}:{original_start}-{original_end}",
                        'target_coords': f"{target_chr}:{target_start}-{target_end}",
                        'reciprocal_coords': f"{reverse_coords_list[0]['chr']}:{reverse_coords_list[0]['start']}-{reverse_coords_list[0]['end']}"
                    }
            else:
                # Forward lifted but reverse lift failed
                if best_quality_metric is None:
                    best_quality_metric = {
                        'exon_id': exon_id,
                        'gene_name': gene_name,
                        'forward_lift_success': True,
                        'reverse_lift_success': False,
                        'reciprocal_match': False,
                        'tolerance': 'N/A',
                        'mapping_type': 'failed',
                        'source_coords': f"{original_chr}:{original_start}-{original_end}",
                        'target_coords': f"{target_chr}:{target_start}-{target_end}",
                        'reciprocal_coords': 'N/A'
                    }

        # Add quality metric for this exon (only once per exon)
        if best_quality_metric is not None:
            quality_metrics.append(best_quality_metric)

        if not has_reciprocal:
            non_reciprocal.append(exon_id)

        # Increment mapping group for 1:many
        if mapping_type == "1:many":
            mapping_group_counter += 1

    return reciprocal_hits, non_reciprocal, one_to_one, one_to_many, quality_metrics


# ============================================================================
# PHASE 3: OUTPUT FILE GENERATION - Status: [X]
# ============================================================================

def write_one_way_outputs(unique_exons, forward_lifted, forward_unmapped, output_dir, source_asm, target_asm):
    """
    Write standard one-way liftover output files.

    Args:
        unique_exons: List of original exon dicts
        forward_lifted: Dict of exon_id -> list of target coords
        forward_unmapped: Dict of exon_id -> unmapped reason
        output_dir: Output directory path
        source_asm: Source assembly name
        target_asm: Target assembly name

    Returns:
        Tuple of (lifted_count, unmapped_count)
    """
    one_way_dir = Path(output_dir) / "one_way"
    one_way_dir.mkdir(parents=True, exist_ok=True)

    lifted_path = one_way_dir / f"{source_asm}_to_{target_asm}_lifted.tsv"
    unmapped_path = one_way_dir / f"{source_asm}_to_{target_asm}_unmapped.tsv"

    lifted_count = 0
    unmapped_count = 0

    # Write lifted file (take first mapping if multiple)
    with open(lifted_path, 'w') as f:
        f.write("exon_id\tgene_name\tsource_chr\tsource_start\tsource_end\tstrand\ttarget_chr\ttarget_start\ttarget_end\n")
        for exon in unique_exons:
            if exon['exon_id'] in forward_lifted:
                target_coords_list = forward_lifted[exon['exon_id']]
                # Write all mappings (for 1:many cases)
                for target_coords in target_coords_list:
                    f.write(f"{exon['exon_id']}\t{exon['gene_name']}\t{exon['chrom']}\t{exon['start']}\t{exon['end']}\t{exon['strand']}\t{target_coords['chr']}\t{target_coords['start']}\t{target_coords['end']}\n")
                    lifted_count += 1

    # Write unmapped file
    with open(unmapped_path, 'w') as f:
        f.write("exon_id\tgene_name\tsource_chr\tsource_start\tsource_end\tstrand\tunmapped_reason\n")
        for exon in unique_exons:
            if exon['exon_id'] in forward_unmapped:
                reason = forward_unmapped[exon['exon_id']]
                f.write(f"{exon['exon_id']}\t{exon['gene_name']}\t{exon['chrom']}\t{exon['start']}\t{exon['end']}\t{exon['strand']}\t{reason}\n")
                unmapped_count += 1

    print(f"\nOne-way output files written:")
    print(f"  Lifted:   {lifted_path} ({lifted_count} mappings)")
    print(f"  Unmapped: {unmapped_path} ({unmapped_count} exons)")

    return lifted_count, unmapped_count


def write_reciprocal_outputs(reciprocal_hits, non_reciprocal, unique_exons, output_dir, source_asm, target_asm):
    """
    Write reciprocal liftover output files.

    Args:
        reciprocal_hits: List of reciprocal hit records
        non_reciprocal: List of exon_ids that didn't reciprocally map
        unique_exons: List of original exon dicts
        output_dir: Output directory path
        source_asm: Source assembly name
        target_asm: Target assembly name

    Returns:
        Tuple of (reciprocal_count, non_reciprocal_count)
    """
    reciprocal_dir = Path(output_dir) / "reciprocal"
    reciprocal_dir.mkdir(parents=True, exist_ok=True)

    lifted_path = reciprocal_dir / f"{source_asm}_to_{target_asm}_reciprocal_lifted.tsv"
    unmapped_path = reciprocal_dir / f"{source_asm}_to_{target_asm}_reciprocal_unmapped.tsv"

    # Write reciprocal lifted file
    with open(lifted_path, 'w') as f:
        f.write("exon_id\tgene_name\tsource_chr\tsource_start\tsource_end\tstrand\t")
        f.write("target_chr\ttarget_start\ttarget_end\t")
        f.write("reciprocal_chr\treciprocal_start\treciprocal_end\t")
        f.write("tolerance\tstart_diff\tend_diff\tmapping_type\n")

        for hit in reciprocal_hits:
            f.write(f"{hit['exon_id']}\t{hit['gene_name']}\t")
            f.write(f"{hit['source_chr']}\t{hit['source_start']}\t{hit['source_end']}\t{hit['strand']}\t")
            f.write(f"{hit['target_chr']}\t{hit['target_start']}\t{hit['target_end']}\t")
            f.write(f"{hit['reciprocal_chr']}\t{hit['reciprocal_start']}\t{hit['reciprocal_end']}\t")
            f.write(f"{hit['tolerance']}\t{hit['start_diff']}\t{hit['end_diff']}\t{hit['mapping_type']}\n")

    # Write reciprocal unmapped file
    exon_lookup = {exon['exon_id']: exon for exon in unique_exons}
    with open(unmapped_path, 'w') as f:
        f.write("exon_id\tgene_name\tsource_chr\tsource_start\tsource_end\tstrand\tunmapped_reason\n")
        for exon_id in non_reciprocal:
            exon = exon_lookup[exon_id]
            f.write(f"{exon_id}\t{exon['gene_name']}\t{exon['chrom']}\t{exon['start']}\t{exon['end']}\t{exon['strand']}\tFailed reciprocal lift\n")

    print(f"\nReciprocal output files written:")
    print(f"  Lifted:   {lifted_path} ({len(reciprocal_hits)} reciprocal mappings)")
    print(f"  Unmapped: {unmapped_path} ({len(non_reciprocal)} exons)")

    return len(reciprocal_hits), len(non_reciprocal)


def write_ortholog_outputs(one_to_one, one_to_many, quality_metrics, output_dir, source_asm, target_asm):
    """
    Write ortholog classification and quality metrics files.

    Args:
        one_to_one: List of 1:1 ortholog records
        one_to_many: List of 1:many ortholog records
        quality_metrics: List of quality metrics records
        output_dir: Output directory path
        source_asm: Source assembly name
        target_asm: Target assembly name

    Returns:
        Tuple of (one_to_one_count, one_to_many_count)
    """
    ortholog_dir = Path(output_dir) / "ortholog_mapping"
    ortholog_dir.mkdir(parents=True, exist_ok=True)

    one_to_one_path = ortholog_dir / f"{source_asm}_to_{target_asm}_1to1_orthologs.tsv"
    one_to_many_path = ortholog_dir / f"{source_asm}_to_{target_asm}_1tomany_orthologs.tsv"
    metrics_path = ortholog_dir / f"{source_asm}_to_{target_asm}_quality_metrics.tsv"

    # Write 1:1 orthologs
    with open(one_to_one_path, 'w') as f:
        f.write("exon_id\tgene_name\tsource_chr\tsource_start\tsource_end\tstrand\t")
        f.write("target_chr\ttarget_start\ttarget_end\ttolerance\tmapping_type\n")

        for record in one_to_one:
            f.write(f"{record['exon_id']}\t{record['gene_name']}\t")
            f.write(f"{record['source_chr']}\t{record['source_start']}\t{record['source_end']}\t{record['strand']}\t")
            f.write(f"{record['target_chr']}\t{record['target_start']}\t{record['target_end']}\t")
            f.write(f"{record['tolerance']}\t{record['mapping_type']}\n")

    # Write 1:many orthologs
    with open(one_to_many_path, 'w') as f:
        f.write("exon_id\tgene_name\tsource_chr\tsource_start\tsource_end\tstrand\t")
        f.write("target_chr\ttarget_start\ttarget_end\ttolerance\tmapping_group\tmapping_type\n")

        for record in one_to_many:
            f.write(f"{record['exon_id']}\t{record['gene_name']}\t")
            f.write(f"{record['source_chr']}\t{record['source_start']}\t{record['source_end']}\t{record['strand']}\t")
            f.write(f"{record['target_chr']}\t{record['target_start']}\t{record['target_end']}\t")
            f.write(f"{record['tolerance']}\t{record['mapping_group']}\t{record['mapping_type']}\n")

    # Write quality metrics
    with open(metrics_path, 'w') as f:
        f.write("exon_id\tgene_name\tforward_lift_success\treverse_lift_success\treciprocal_match\t")
        f.write("tolerance\tmapping_type\tsource_coords\ttarget_coords\treciprocal_coords\n")

        for record in quality_metrics:
            f.write(f"{record['exon_id']}\t{record['gene_name']}\t")
            f.write(f"{record['forward_lift_success']}\t{record['reverse_lift_success']}\t{record['reciprocal_match']}\t")
            f.write(f"{record['tolerance']}\t{record['mapping_type']}\t")
            f.write(f"{record['source_coords']}\t{record['target_coords']}\t{record['reciprocal_coords']}\n")

    print(f"\nOrtholog mapping files written:")
    print(f"  1:1 orthologs:    {one_to_one_path} ({len(one_to_one)} mappings)")
    print(f"  1:many orthologs: {one_to_many_path} ({len(one_to_many)} mappings)")
    print(f"  Quality metrics:  {metrics_path} ({len(quality_metrics)} exons)")

    return len(one_to_one), len(one_to_many)

# ============================================================================
# PHASE 4: MAIN WORKFLOW INTEGRATION - Status: [X]
# ============================================================================

def main():
    """Main execution function."""
    args = parse_args()

    print("=" * 80)
    print("RECIPROCAL LIFTOVER WORKFLOW")
    print("=" * 80)
    print(f"Input:            {args.input}")
    print(f"Forward chain:    {args.forward_chain}")
    print(f"Reverse chain:    {args.reverse_chain}")
    print(f"Output dir:       {args.output_dir}")
    print(f"Source assembly:  {args.source_assembly}")
    print(f"Target assembly:  {args.target_assembly}")
    print(f"Min match:        {args.min_match}")
    if args.conda_env:
        print(f"Conda env:        {args.conda_env}")
    print()

    # Step 1: Parse input
    print("Step 1: Parsing input TSV...")
    exons = parse_input_tsv(args.input)
    print(f"  Total rows read: {len(exons):,}")

    # Step 2: Deduplicate
    print("\nStep 2: Deduplicating exons...")
    unique_exons = deduplicate_exons(exons)
    print(f"  Unique exons: {len(unique_exons):,}")
    print(f"  Duplicates removed: {len(exons) - len(unique_exons):,}")

    # Step 3: Create output directory structure
    print("\nStep 3: Creating output directories...")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"  Output directory: {output_dir}")

    # Use temporary directory for intermediate files
    with tempfile.TemporaryDirectory() as tmpdir:
        print(f"\nStep 4: Running forward liftOver (with -multiple flag)...")

        # Create BED file for forward lift
        forward_input_bed = Path(tmpdir) / "input.bed"
        forward_output_bed = Path(tmpdir) / "forward_output.bed"
        forward_unmapped_bed = Path(tmpdir) / "forward_unmapped.bed"

        write_bed_file(unique_exons, forward_input_bed)

        # Run forward liftOver with -multiple flag
        success = run_liftover(
            forward_input_bed, args.forward_chain, forward_output_bed, forward_unmapped_bed,
            args.min_match, args.conda_env, use_multiple=True
        )

        if not success:
            print("Error: Forward liftOver failed", file=sys.stderr)
            sys.exit(1)

        # Parse forward lift results
        print("\nStep 5: Parsing forward liftOver results...")
        forward_lifted = parse_lifted_bed(forward_output_bed, allow_multiple=True)
        forward_unmapped = parse_unmapped_bed(forward_unmapped_bed, args.min_match)

        unique_forward_exons = len(forward_lifted)
        total_forward_mappings = sum(len(coords_list) for coords_list in forward_lifted.values())
        print(f"  Unique exons lifted: {unique_forward_exons:,}")
        print(f"  Total mappings (including 1:many): {total_forward_mappings:,}")
        print(f"  Unmapped: {len(forward_unmapped):,}")

        # Step 6: Write one-way outputs
        print("\nStep 6: Writing one-way liftOver outputs...")
        write_one_way_outputs(
            unique_exons, forward_lifted, forward_unmapped,
            args.output_dir, args.source_assembly, args.target_assembly
        )

        # Step 7: Run reverse liftOver
        print("\nStep 7: Running reverse liftOver (reciprocal mapping)...")
        reverse_lifted = run_reciprocal_liftover(
            forward_lifted, args.reverse_chain, tmpdir,
            args.min_match, args.conda_env
        )
        print(f"  Reverse mappings obtained: {len(reverse_lifted):,}")

        # Step 8: Classify orthologs
        print("\nStep 8: Classifying ortholog relationships...")
        reciprocal_hits, non_reciprocal, one_to_one, one_to_many, quality_metrics = \
            classify_ortholog_mapping(unique_exons, forward_lifted, reverse_lifted)

        print(f"  Reciprocal hits: {len(reciprocal_hits):,}")
        print(f"  Non-reciprocal: {len(non_reciprocal):,}")
        print(f"  1:1 orthologs: {len(one_to_one):,}")
        print(f"  1:many orthologs: {len(one_to_many):,}")

        # Step 9: Write reciprocal outputs
        print("\nStep 9: Writing reciprocal liftOver outputs...")
        write_reciprocal_outputs(
            reciprocal_hits, non_reciprocal, unique_exons,
            args.output_dir, args.source_assembly, args.target_assembly
        )

        # Step 10: Write ortholog outputs
        print("\nStep 10: Writing ortholog mapping outputs...")
        write_ortholog_outputs(
            one_to_one, one_to_many, quality_metrics,
            args.output_dir, args.source_assembly, args.target_assembly
        )

    # Step 11: Calculate and print summary statistics
    print("\n" + "=" * 80)
    print("SUMMARY STATISTICS")
    print("=" * 80)

    total = len(unique_exons)
    forward_success = len(forward_lifted)
    reciprocal_success = len(reciprocal_hits)

    print(f"\nTotal unique exons:          {total:,}")
    print(f"Forward lift success:        {forward_success:,} ({100*forward_success/total:.2f}%)")
    print(f"Reverse lift success:        {len(reverse_lifted):,}")
    print(f"Reciprocal matches:          {reciprocal_success:,} ({100*reciprocal_success/total:.2f}%)")

    # Tolerance distribution
    print(f"\nTolerance distribution (reciprocal hits):")
    tolerance_counts = {'0': 0, '1-2': 0, '3-5': 0, '6-10': 0, '>10': 0}

    for hit in reciprocal_hits:
        tol = hit['tolerance']
        if tol == '0':
            tolerance_counts['0'] += 1
        elif tol in ['1', '2']:
            tolerance_counts['1-2'] += 1
        elif tol in ['3', '4', '5']:
            tolerance_counts['3-5'] += 1
        elif tol in ['6', '7', '8', '9', '10']:
            tolerance_counts['6-10'] += 1
        elif tol == '>10':
            tolerance_counts['>10'] += 1

    if reciprocal_success > 0:
        print(f"  Exact match (0bp):           {tolerance_counts['0']:,} ({100*tolerance_counts['0']/reciprocal_success:.2f}%)")
        print(f"  ±1-2bp:                      {tolerance_counts['1-2']:,} ({100*tolerance_counts['1-2']/reciprocal_success:.2f}%)")
        print(f"  ±3-5bp:                      {tolerance_counts['3-5']:,} ({100*tolerance_counts['3-5']/reciprocal_success:.2f}%)")
        print(f"  ±6-10bp:                     {tolerance_counts['6-10']:,} ({100*tolerance_counts['6-10']/reciprocal_success:.2f}%)")
        print(f"  >10bp:                       {tolerance_counts['>10']:,} ({100*tolerance_counts['>10']/reciprocal_success:.2f}%)")

    # Ortholog classification
    print(f"\nOrtholog classification:")
    print(f"  1:1 orthologs:               {len(one_to_one):,} ({100*len(one_to_one)/total:.2f}% of total)")
    print(f"  1:many orthologs:            {len(one_to_many):,} ({100*len(one_to_many)/total:.2f}% of total)")
    print(f"  Failed mapping:              {len(non_reciprocal):,} ({100*len(non_reciprocal)/total:.2f}% of total)")

    # Sanity checks
    print(f"\nSanity checks:")
    if len(quality_metrics) == total:
        print(f"  ✓ Quality metrics covers all {total:,} input exons")
    else:
        print(f"  ✗ Warning: Quality metrics has {len(quality_metrics):,} entries but {total:,} input exons")

    if forward_success + len(forward_unmapped) == total:
        print(f"  ✓ Forward lift accounting correct")
    else:
        print(f"  ✗ Warning: Forward lift accounting mismatch")

    print("\n" + "=" * 80)
    print("RECIPROCAL LIFTOVER COMPLETE")
    print("=" * 80)
    print(f"\nOutput directory: {args.output_dir}")
    print(f"  - one_way/                 (standard liftover results)")
    print(f"  - reciprocal/              (reciprocal liftover results)")
    print(f"  - ortholog_mapping/        (1:1 and 1:many orthologs + quality metrics)")
    print()


if __name__ == '__main__':
    main()
