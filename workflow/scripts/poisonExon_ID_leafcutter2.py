#!/usr/bin/env python3
"""
merged_poison_exon_pipeline.py
Integrated Poison Exon Detection Pipeline - Leafcutter2 Version

Combines three sequential analysis steps:
1. oct2025_ID_plusCryptic.py - Identify poison exons and cryptic splice sites
2. complete_poison_junction.py - Analyze poison exon clusters
3. BED_color_convert.py - Convert to colored BED format

This version is optimized for Leafcutter2 cluster_ratios.gz format
and supports both GENCODE (gene_name) and Ensembl (gene_id) GTF annotations.

Author: Dylan Stermer
Date: October 31, 2025
Updated: February 2026 (Leafcutter2 compatibility)
"""

import argparse
import sys
import pandas as pd
import gzip
from collections import defaultdict
import os
from pathlib import Path


def parse_gtf_attribute(info_string, attribute_name, fallback_attribute=None):
    """
    Safely extract an attribute from GTF info field.
    Returns the attribute value if found, falls back to fallback_attribute, or None.

    This function is designed to handle both GENCODE and Ensembl GTF formats.
    GENCODE GTFs typically have gene_name for all entries, while Ensembl GTFs
    may only have gene_id for some entries.
    """
    pattern = f'{attribute_name} "'
    if pattern in info_string:
        return info_string.split(pattern)[1].split('"')[0]
    elif fallback_attribute:
        pattern = f'{fallback_attribute} "'
        if pattern in info_string:
            return info_string.split(pattern)[1].split('"')[0]
    return None


def parse_args(args=None):
    parser = argparse.ArgumentParser(
        description="Integrated poison exon detection and analysis pipeline (Leafcutter2 version)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Inputs
    parser.add_argument('-g', '--gtf', required=True,
                        help='Input GTF file path (e.g., GRCh38 Gencode or Ensembl)')
    parser.add_argument('-l', '--leafcutter2-annotations', required=True,
                        help='Leafcutter2 junction classification file path')
    parser.add_argument('-c', '--leafcutter-counts', required=True,
                        help='Path to Leafcutter2 cluster_ratios file (gzipped)')
    parser.add_argument('-o', '--output-prefix', required=True,
                        help='Output prefix for files (base directory and sample name)')

    # Parameters
    parser.add_argument('-d', '--max-distance', type=int, default=500,
                        help='Maximum distance between non-coding introns to consider as poison exon')

    return parser.parse_args(args)


def setup_directories(output_prefix):
    """Create organized output directory structure."""
    base_dir = os.path.dirname(output_prefix)
    if not base_dir:
        base_dir = '.'

    pe_dir = os.path.join(base_dir, 'PE')
    colored_bed_dir = os.path.join(base_dir, 'ColoredBed')

    os.makedirs(pe_dir, exist_ok=True)
    os.makedirs(colored_bed_dir, exist_ok=True)

    return pe_dir, colored_bed_dir


# ============================================================================
# STEP 1: oct2025_ID_plusCryptic.py functionality
# ============================================================================

def process_gtf_to_dict(gtf_path):
    """Process GTF file and return dictionary for poison exon detection.

    This function now supports both GENCODE and Ensembl GTF formats by attempting
    to extract gene_name first, and falling back to gene_id if gene_name is not present.
    """
    dict_cds = {}
    skipped_count = 0
    processed_count = 0

    print(f"[Step 1] Processing GTF file for poison exon detection: {gtf_path}", file=sys.stderr)
    for line in open(gtf_path):
        if line[0] == "#":
            continue
        chr, source, ann_type, start, stop, dot, strand, frame, info = line.split("\t")
        if ann_type != "CDS":
            continue

        # Try gene_name first, fallback to gene_id for Ensembl GTFs
        gene_identifier = parse_gtf_attribute(info, 'gene_name', 'gene_id')

        if gene_identifier is None:
            skipped_count += 1
            continue

        if gene_identifier not in dict_cds:
            dict_cds[gene_identifier] = set()
        dict_cds[gene_identifier].add((chr, int(start), int(stop), strand))
        processed_count += 1

    print(f"[Step 1] Processed {processed_count} CDS entries for {len(dict_cds)} genes", file=sys.stderr)
    if skipped_count > 0:
        print(f"[Step 1] Skipped {skipped_count} CDS entries without gene identifiers", file=sys.stderr)

    return dict_cds


def process_gtf_to_df(gtf_path):
    """Process GTF file and return DataFrame for filtering.

    Extracts all CDS regions for filtering poison exons that overlap with coding sequences.
    """
    cds_records = []

    print(f"[Step 1] Processing GTF file for CDS filtering: {gtf_path}", file=sys.stderr)
    for line in open(gtf_path):
        if line[0] == "#":
            continue
        chr, source, ann_type, start, stop, dot, strand, frame, info = line.split("\t")
        if ann_type != "CDS":
            continue
        cds_records.append({
            '#chr': chr,
            'start': int(start),
            'stop': int(stop)
        })

    print(f"[Step 1] Found {len(cds_records)} CDS regions for overlap filtering", file=sys.stderr)
    return pd.DataFrame(cds_records)


def process_leafcutter(leafcutter_path):
    """Process Leafcutter2 annotations."""
    dict_leafcutter2 = {}

    print(f"[Step 1] Processing Leafcutter2 annotations: {leafcutter_path}", file=sys.stderr)
    for line in open(leafcutter_path):
        if line.strip().startswith("Gene_name"):
            continue
        gene_name, intron_cord, strand, annotation, coding, UTR, GencodePC = line.strip().split("\t")
        chrom, positions = intron_cord.split(':')
        start, stop = positions.split('-')
        if coding == "False" and UTR == "False":
            if gene_name not in dict_leafcutter2:
                dict_leafcutter2[gene_name] = set()

            if strand == "+":
                stop = int(stop) + 1
                dict_leafcutter2[gene_name].add((chrom, int(start), stop, annotation))
            elif strand == "-":
                stop = int(stop) + 1
                dict_leafcutter2[gene_name].add((chrom, int(start), stop, annotation))

    return dict_leafcutter2


def detect_poison_exons(dict_cds, dict_leafcutter2, max_distance):
    """Detect potential poison exons."""
    dict_poison_5_junction = {}
    dict_poison_3_junction = {}

    print("[Step 1] Identifying potential poison exons...", file=sys.stderr)
    for gene_name in dict_leafcutter2:
        if gene_name not in dict_cds:
            continue

        coding_exons = dict_cds[gene_name]
        non_code_introns = dict_leafcutter2[gene_name]

        for non_code_chr, non_code_start, non_code_stop, non_code_annote in non_code_introns:
            for coding_chr, coding_start, coding_stop, coding_strand in coding_exons:
                if non_code_chr != coding_chr:
                    continue

                offset = 0

                if non_code_start == coding_stop:
                    if gene_name not in dict_poison_5_junction:
                        dict_poison_5_junction[gene_name] = set()
                    coding_5_exon = f"{coding_chr}:{coding_start}-{coding_stop}"
                    dict_poison_5_junction[gene_name].add(
                        (non_code_stop, non_code_start, non_code_chr, coding_5_exon, coding_strand)
                    )
                if non_code_stop == coding_start - offset:
                    if gene_name not in dict_poison_3_junction:
                        dict_poison_3_junction[gene_name] = set()
                    coding_3_exon = f"{coding_chr}:{coding_start}-{coding_stop}"
                    dict_poison_3_junction[gene_name].add(
                        (non_code_start, non_code_stop, non_code_chr, coding_3_exon, coding_strand)
                    )

    return dict_poison_5_junction, dict_poison_3_junction


def remove_overlapping_exons(poison_df, cds_df):
    """Remove poison exons that overlap any CDS region.

    This function filters out poison exons that overlap with coding sequences,
    as these are less likely to be true poison exons (which should be in non-coding regions).
    """
    # Handle empty inputs
    if poison_df.empty:
        return poison_df
    if cds_df.empty:
        print("[Step 1] Warning: No CDS regions found in GTF", file=sys.stderr)
        return poison_df

    filtered_exons = []

    for chrom in poison_df['#chr'].unique():
        poison_group = poison_df[poison_df['#chr'] == chrom]
        cds_group = cds_df[cds_df['#chr'] == chrom]

        if cds_group.empty:
            filtered_exons.append(poison_group)
            continue

        for _, exon in poison_group.iterrows():
            overlap = cds_group[(cds_group['start'] < exon['end']) & (cds_group['stop'] > exon['start'])]
            if overlap.empty:
                filtered_exons.append(exon)

    if not filtered_exons:
        return pd.DataFrame(columns=poison_df.columns)

    return pd.DataFrame(filtered_exons)


def write_poison_exons(output_path, dict_poison_5_junction, dict_poison_3_junction, max_distance):
    """Write poison exons to BED file and return as DataFrame."""
    seen_records = set()
    poison_records = []

    with open(output_path, 'w') as bed:
        bed.write('#chr\tstart\tend\tname\tlength\tstrand\tupstream_intron\tdownstream_intron\n')

        for gene_name in dict_poison_5_junction:
            if gene_name not in dict_poison_3_junction:
                continue

            non_coding_upstream_intron = list(dict_poison_5_junction[gene_name])
            non_coding_upstream_intron.sort()
            non_coding_downstream_intron = list(dict_poison_3_junction[gene_name])
            non_coding_downstream_intron.sort()

            for up_end, up_start, up_chr, up_coding_exon, coding_strand in non_coding_upstream_intron:
                for down_start, down_end, down_chr, down_coding_exon, coding_strand in non_coding_downstream_intron:
                    if (up_chr != down_chr or
                        int(up_end+1) == int(down_end-1) or
                        int(up_start) == int(down_start) or
                        int(up_end) > int(down_start)):
                        continue
                    if abs(up_end - down_start) <= max_distance:

                        upstream_intron = f"{up_chr}:{up_start}-{up_end-1}"
                        downstream_intron = f"{down_chr}:{down_start}-{down_end-1}"

                        record_key = (up_chr, up_end, down_start, gene_name, coding_strand,
                                    upstream_intron, downstream_intron)

                        if record_key not in seen_records:
                            seen_records.add(record_key)

                            record = {
                                '#chr': up_chr,
                                'start': up_end,
                                'end': down_start,
                                'name': gene_name,
                                'length': abs(up_end - down_start),
                                'strand': coding_strand,
                                'upstream_intron': upstream_intron,
                                'downstream_intron': downstream_intron
                            }
                            poison_records.append(record)

                            bed.write(
                                f'{up_chr}\t{up_end}\t{down_start}\t{gene_name}\t'
                                f'{abs(up_end - down_start)}\t{coding_strand}\t'
                                f'{upstream_intron}\t{downstream_intron}\n'
                            )

    return pd.DataFrame(poison_records)


def write_potential_cryptic_sites(output_path, dict_poison_5_junction, dict_poison_3_junction):
    """Write potential cryptic splice sites (single-sided junctions) to BED."""
    seen_records = set()
    cryptic_records = []

    with open(output_path, 'w') as bed:
        bed.write('#chr\tstart\tend\tname\tstrand\tjunction_type\tcryptic_junction\n')

        all_genes = set(dict_poison_5_junction.keys()) | set(dict_poison_3_junction.keys())

        for gene_name in all_genes:
            has_5 = gene_name in dict_poison_5_junction
            has_3 = gene_name in dict_poison_3_junction

            # Only keep genes with one-sided events
            if has_5 and has_3:
                continue

            if has_5:
                for non_code_stop, non_code_start, non_code_chr, coding_exon, strand in dict_poison_5_junction[gene_name]:
                    record_key = (non_code_chr, non_code_start, non_code_stop, gene_name, '5prime')
                    if record_key in seen_records:
                        continue
                    seen_records.add(record_key)

                    cryptic_junction = f"{non_code_chr}:{non_code_start}-{non_code_stop}"
                    record = {
                        '#chr': non_code_chr,
                        'start': non_code_start,
                        'end': non_code_stop,
                        'name': gene_name,
                        'strand': strand,
                        'junction_type': '5prime',
                        'cryptic_junction': cryptic_junction
                    }
                    cryptic_records.append(record)
                    bed.write(
                        f"{non_code_chr}\t{non_code_start}\t{non_code_stop}\t{gene_name}"
                        f"\t{strand}\t5prime\t{cryptic_junction}\n"
                    )

            if has_3:
                for non_code_start, non_code_stop, non_code_chr, coding_exon, strand in dict_poison_3_junction[gene_name]:
                    record_key = (non_code_chr, non_code_start, non_code_stop, gene_name, '3prime')
                    if record_key in seen_records:
                        continue
                    seen_records.add(record_key)

                    cryptic_junction = f"{non_code_chr}:{non_code_start}-{non_code_stop}"
                    record = {
                        '#chr': non_code_chr,
                        'start': non_code_start,
                        'end': non_code_stop,
                        'name': gene_name,
                        'strand': strand,
                        'junction_type': '3prime',
                        'cryptic_junction': cryptic_junction
                    }
                    cryptic_records.append(record)
                    bed.write(
                        f"{non_code_chr}\t{non_code_start}\t{non_code_stop}\t{gene_name}"
                        f"\t{strand}\t3prime\t{cryptic_junction}\n"
                    )

    return pd.DataFrame(cryptic_records)


def run_step1(args, pe_dir):
    """Execute Step 1: Poison exon identification and cryptic splice site detection."""
    print("\n" + "="*80, file=sys.stderr)
    print("STEP 1: Poison Exon Identification and Cryptic Splice Site Detection", file=sys.stderr)
    print("="*80 + "\n", file=sys.stderr)

    # Get sample name from output prefix
    sample_name = os.path.basename(args.output_prefix)

    # Process GTF file in two ways
    dict_cds = process_gtf_to_dict(args.gtf)
    cds_df = process_gtf_to_df(args.gtf)

    # Process Leafcutter2 file
    dict_leafcutter2 = process_leafcutter(args.leafcutter2_annotations)

    # Diagnostic logging: check overlap between GTF and leafcutter2 gene identifiers
    gtf_genes = set(dict_cds.keys())
    leafcutter_genes = set(dict_leafcutter2.keys())
    overlapping_genes = gtf_genes & leafcutter_genes

    print(f"[Step 1] GTF genes with CDS: {len(gtf_genes)}", file=sys.stderr)
    print(f"[Step 1] Leafcutter2 genes with non-coding introns: {len(leafcutter_genes)}", file=sys.stderr)
    print(f"[Step 1] Overlapping genes: {len(overlapping_genes)}", file=sys.stderr)

    if overlapping_genes:
        sample_genes = list(overlapping_genes)[:5]
        print(f"[Step 1] Sample overlapping genes: {', '.join(sample_genes)}", file=sys.stderr)
    else:
        print("[Step 1] WARNING: No overlapping genes found between GTF and leafcutter2!", file=sys.stderr)
        print("[Step 1] Sample GTF genes: {}".format(', '.join(list(gtf_genes)[:3])), file=sys.stderr)
        print("[Step 1] Sample leafcutter2 genes: {}".format(', '.join(list(leafcutter_genes)[:3])), file=sys.stderr)

    # Detect poison exons
    dict_poison_5_junction, dict_poison_3_junction = detect_poison_exons(
        dict_cds, dict_leafcutter2, args.max_distance
    )

    # Define output paths in PE directory
    all_output = os.path.join(pe_dir, f"{sample_name}_all.bed")
    filtered_output = os.path.join(pe_dir, f"{sample_name}_filtered.bed")
    cryptic_output = os.path.join(pe_dir, f"{sample_name}_potential_cryptic_ss_all.bed")

    # Write all poison exons
    print(f"[Step 1] Writing all poison exons to: {all_output}", file=sys.stderr)
    poison_df = write_poison_exons(
        all_output,
        dict_poison_5_junction,
        dict_poison_3_junction,
        args.max_distance
    )

    # Write potential cryptic splice sites
    print(f"[Step 1] Writing potential cryptic splice site candidates to: {cryptic_output}", file=sys.stderr)
    cryptic_df = write_potential_cryptic_sites(
        cryptic_output,
        dict_poison_5_junction,
        dict_poison_3_junction
    )

    # Filter and write non-overlapping poison exons
    if poison_df.empty:
        print("[Step 1] No poison exons detected. This could mean:", file=sys.stderr)
        print("  - No suitable gene candidates found", file=sys.stderr)
        print("  - Gene identifiers don't match between GTF and leafcutter annotations", file=sys.stderr)
        print("  - No non-coding introns adjacent to CDS regions", file=sys.stderr)
        # Create empty filtered file with header
        with open(filtered_output, 'w') as f:
            f.write('#chr\tstart\tend\tname\tlength\tstrand\tupstream_intron\tdownstream_intron\n')
        filtered_df = poison_df
    else:
        print("[Step 1] Filtering poison exons that overlap with CDS regions...", file=sys.stderr)
        filtered_df = remove_overlapping_exons(poison_df, cds_df)

        if not filtered_df.empty:
            filtered_df = filtered_df.drop_duplicates()
            filtered_df.to_csv(filtered_output, sep='\t', index=False)
            print(f"[Step 1] Successfully filtered poison exons. Saved to {filtered_output}", file=sys.stderr)
            print(f"[Step 1] Retained {len(filtered_df)} out of {len(poison_df)} original regions", file=sys.stderr)
        else:
            print("[Step 1] No poison exons remained after CDS overlap filtering", file=sys.stderr)
            # Create empty filtered file with header
            with open(filtered_output, 'w') as f:
                f.write('#chr\tstart\tend\tname\tlength\tstrand\tupstream_intron\tdownstream_intron\n')

    print("[Step 1] Done!\n", file=sys.stderr)

    return all_output, filtered_output, cryptic_output


# ============================================================================
# STEP 2: complete_poison_junction.py functionality
# ============================================================================

def parse_leafcutter_counts(input_file, output_file):
    """Parses Leafcutter2 cluster_ratios file into structured format.

    This function handles the Leafcutter2 format: chr:start:end:clu_X_strand:annotation_type
    Coordinates are kept as-is to match the junction classifications format.
    The extra annotation_type field (e.g., ":IN") is automatically handled by the parsing logic.
    """
    print(f"[Step 2] Parsing Leafcutter2 cluster_ratios from: {input_file}", file=sys.stderr)

    if not Path(input_file).exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    Path(output_file).parent.mkdir(parents=True, exist_ok=True)

    with gzip.open(input_file, 'rt') as infile, \
         gzip.open(output_file, 'wt') as outfile:

        header = infile.readline().strip().split()
        fieldnames = ['original_string', '#chr', 'start', 'stop', 'cluster', 'strand', 'intron_coord'] + header
        outfile.write('\t'.join(fieldnames) + '\n')

        for line in infile:
            row = line.strip().split()
            original_string = row[0]
            intron_info = original_string.split(':')

            chr_name = intron_info[0]
            start = intron_info[1]
            stop = intron_info[2]  # Keep original coordinate from leafcutter2
            cluster = f'clu_{intron_info[3].split("_")[1]}'
            strand = intron_info[3].split('_')[2]
            intron_coord = f'{chr_name}:{start}-{stop}'

            parsed_row = [original_string, chr_name, start, stop, cluster, strand, intron_coord] + row[1:]
            outfile.write('\t'.join(parsed_row) + '\n')

    return output_file


def process_leafcutter_data(parsed_counts_file):
    """Creates cluster and junction mappings from parsed LeafCutter data."""
    dict_leafcutter_by_cluster = {}
    dict_leafcutter_by_junction = {}

    with gzip.open(parsed_counts_file, 'rt') as f:
        next(f)  # Skip header
        for line in f:
            fields = line.strip().split('\t')
            cluster = fields[4]
            intron_coord = fields[6]

            if cluster not in dict_leafcutter_by_cluster:
                dict_leafcutter_by_cluster[cluster] = set()
            dict_leafcutter_by_cluster[cluster].add(intron_coord)

            if intron_coord not in dict_leafcutter_by_junction:
                dict_leafcutter_by_junction[intron_coord] = []
            dict_leafcutter_by_junction[intron_coord].append(cluster)

    return dict_leafcutter_by_cluster, dict_leafcutter_by_junction


def process_leafcutter2_annotations(leafcutter2_path):
    """Process Leafcutter2 annotations and create mapping of junctions to coding status."""
    dict_junction_coding_status = {}

    print(f"[Step 2] Processing Leafcutter2 annotations: {leafcutter2_path}", file=sys.stderr)

    for line in open(leafcutter2_path):
        if line.strip().startswith("Gene_name"):
            continue

        fields = line.strip().split("\t")
        intron_coord = fields[1]
        coding = fields[4]

        dict_junction_coding_status[intron_coord] = (coding == "True")

    return dict_junction_coding_status


def process_poison_exons(poison_exon_path, dict_leafcutter_by_junction):
    """Identifies potential coding junctions from poison exon BED file."""
    dict_potential_coding_junction_by_cluster = {}
    dict_potential_coding_junction_with_surroundig_introns = {}

    for line in open(poison_exon_path):
        if line.strip().startswith("#"):
            continue

        columns = line.strip().split("\t")
        upstream_intron = columns[6]
        downstream_intron = columns[7]

        if upstream_intron not in dict_leafcutter_by_junction or \
           downstream_intron not in dict_leafcutter_by_junction:
            continue

        chrom = upstream_intron.split(":")[0]
        upstream_start = upstream_intron.split(":")[1].split("-")[0]
        upstream_end = upstream_intron.split(":")[1].split("-")[1]
        downstream_start = downstream_intron.split(":")[1].split("-")[0]
        downstream_end = downstream_intron.split(":")[1].split("-")[1]

        potential_coding_junction = f"{chrom}:{upstream_start}-{downstream_end}"

        if int(upstream_end) == int(downstream_end):
            continue

        if dict_leafcutter_by_junction[upstream_intron] == dict_leafcutter_by_junction[downstream_intron]:
            cluster = dict_leafcutter_by_junction[upstream_intron][0]

            if cluster not in dict_potential_coding_junction_by_cluster:
                dict_potential_coding_junction_by_cluster[cluster] = set()
            dict_potential_coding_junction_by_cluster[cluster].add(potential_coding_junction)

            if potential_coding_junction not in dict_potential_coding_junction_with_surroundig_introns:
                dict_potential_coding_junction_with_surroundig_introns[potential_coding_junction] = set()
            dict_potential_coding_junction_with_surroundig_introns[potential_coding_junction].add(
                (upstream_intron, downstream_intron)
            )

    return dict_potential_coding_junction_by_cluster, dict_potential_coding_junction_with_surroundig_introns


def find_complete_pe_clusters(dict_leafcutter_by_cluster, dict_potential_coding_junction_by_cluster,
                            dict_junction_coding_status):
    """Identifies clusters containing both poison exons and coding junctions."""
    dict_complete_pe_by_cluster = {}

    for cluster in dict_leafcutter_by_cluster:
        if cluster not in dict_potential_coding_junction_by_cluster:
            continue

        potential_coding_junction = dict_potential_coding_junction_by_cluster[cluster]
        cluster_junctions = dict_leafcutter_by_cluster[cluster]

        for potential_junction in potential_coding_junction:
            for known_junction in cluster_junctions:
                if potential_junction == known_junction:
                    # Check if this junction is annotated as coding
                    if known_junction in dict_junction_coding_status and dict_junction_coding_status[known_junction]:
                        if cluster not in dict_complete_pe_by_cluster:
                            dict_complete_pe_by_cluster[cluster] = []
                        dict_complete_pe_by_cluster[cluster].append(
                            (cluster_junctions, potential_coding_junction)
                        )
                    break

    return dict_complete_pe_by_cluster


def get_pe_dictionary(dict_complete_pe_by_cluster, dict_potential_coding_junction_with_surroundig_introns):
    """Creates mapping of coding junctions to surrounding introns."""
    dict_PE = {}

    for cluster, junctions in dict_complete_pe_by_cluster.items():
        all_junction = junctions[0][0]
        coding_junction = next(iter(junctions[0][1]))

        if coding_junction not in dict_potential_coding_junction_with_surroundig_introns:
            continue

        surrounding_introns = dict_potential_coding_junction_with_surroundig_introns[coding_junction]
        dict_PE[coding_junction] = surrounding_introns

    return dict_PE


def save_step2_outputs(dict_complete_pe_by_cluster, dict_PE, pe_dir, sample_name):
    """Save Step 2 analysis results to output files."""
    complete_pe_path = os.path.join(pe_dir, f"{sample_name}_complete_PE.tsv")
    introns_path = os.path.join(pe_dir, f"{sample_name}_coding_noncoding_introns.tsv")

    with open(complete_pe_path, 'w') as f:
        f.write("cluster\tjunctions\tcoding_junction\n")
        for cluster, data_list in dict_complete_pe_by_cluster.items():
            for junctions, coding_junction in data_list:
                junction_str = ','.join(list(junctions))
                coding_str = next(iter(coding_junction))
                f.write(f"{cluster}\t{junction_str}\t{coding_str}\n")

    with open(introns_path, 'w') as f:
        f.write("coding_junction\tupstream_intron\tdownstream_intron\n")
        for coding_junction, intron_pairs in dict_PE.items():
            for upstream_intron, downstream_intron in intron_pairs:
                f.write(f"{coding_junction}\t{upstream_intron}\t{downstream_intron}\n")

    print(f"[Step 2] Complete poison exon data saved to: {complete_pe_path}", file=sys.stderr)
    print(f"[Step 2] Introns data saved to: {introns_path}", file=sys.stderr)

    return complete_pe_path, introns_path


def run_step2(args, all_bed, filtered_bed, pe_dir):
    """Execute Step 2: Complete poison junction analysis for both all and filtered BED files."""
    print("\n" + "="*80, file=sys.stderr)
    print("STEP 2: Complete Poison Junction Analysis", file=sys.stderr)
    print("="*80 + "\n", file=sys.stderr)

    sample_name = os.path.basename(args.output_prefix)

    # Parse LeafCutter counts (only need to do once)
    parsed_counts_file = os.path.join(pe_dir, f"{sample_name}_parsed_counts.gz")
    parsed_file = parse_leafcutter_counts(args.leafcutter_counts, parsed_counts_file)
    print(f"[Step 2] Parsed counts saved to: {parsed_counts_file}", file=sys.stderr)

    # Process parsed LeafCutter data (only need to do once)
    dict_leafcutter_by_cluster, dict_leafcutter_by_junction = process_leafcutter_data(parsed_file)

    # Process Leafcutter2 annotations (only need to do once)
    dict_junction_coding_status = process_leafcutter2_annotations(args.leafcutter2_annotations)

    # Initialize output file paths
    introns_file_all = None
    introns_file_filtered = None

    # Process "all" poison exons
    if all_bed and os.path.exists(all_bed):
        print(f"[Step 2] Processing all poison exons from: {all_bed}", file=sys.stderr)

        dict_potential_coding_junction_by_cluster, dict_potential_coding_junction_with_surroundig_introns = \
            process_poison_exons(all_bed, dict_leafcutter_by_junction)

        dict_complete_pe_by_cluster = find_complete_pe_clusters(
            dict_leafcutter_by_cluster,
            dict_potential_coding_junction_by_cluster,
            dict_junction_coding_status
        )

        dict_PE = get_pe_dictionary(
            dict_complete_pe_by_cluster,
            dict_potential_coding_junction_with_surroundig_introns
        )

        # Save outputs with "_all" suffix
        complete_pe_path_all = os.path.join(pe_dir, f"{sample_name}_complete_PE_all.tsv")
        introns_file_all = os.path.join(pe_dir, f"{sample_name}_coding_noncoding_introns_all.tsv")

        with open(complete_pe_path_all, 'w') as f:
            f.write("cluster\tjunctions\tcoding_junction\n")
            for cluster, data_list in dict_complete_pe_by_cluster.items():
                for junctions, coding_junction in data_list:
                    junction_str = ','.join(list(junctions))
                    coding_str = next(iter(coding_junction))
                    f.write(f"{cluster}\t{junction_str}\t{coding_str}\n")

        with open(introns_file_all, 'w') as f:
            f.write("coding_junction\tupstream_intron\tdownstream_intron\n")
            for coding_junction, intron_pairs in dict_PE.items():
                for upstream_intron, downstream_intron in intron_pairs:
                    f.write(f"{coding_junction}\t{upstream_intron}\t{downstream_intron}\n")

        print(f"[Step 2] Complete PE (all) saved to: {complete_pe_path_all}", file=sys.stderr)
        print(f"[Step 2] Introns (all) saved to: {introns_file_all}", file=sys.stderr)

    # Process "filtered" poison exons
    if filtered_bed and os.path.exists(filtered_bed):
        print(f"[Step 2] Processing filtered poison exons from: {filtered_bed}", file=sys.stderr)

        dict_potential_coding_junction_by_cluster, dict_potential_coding_junction_with_surroundig_introns = \
            process_poison_exons(filtered_bed, dict_leafcutter_by_junction)

        dict_complete_pe_by_cluster = find_complete_pe_clusters(
            dict_leafcutter_by_cluster,
            dict_potential_coding_junction_by_cluster,
            dict_junction_coding_status
        )

        dict_PE = get_pe_dictionary(
            dict_complete_pe_by_cluster,
            dict_potential_coding_junction_with_surroundig_introns
        )

        # Save outputs with "_filtered" suffix
        complete_pe_path_filtered = os.path.join(pe_dir, f"{sample_name}_complete_PE_filtered.tsv")
        introns_file_filtered = os.path.join(pe_dir, f"{sample_name}_coding_noncoding_introns_filtered.tsv")

        with open(complete_pe_path_filtered, 'w') as f:
            f.write("cluster\tjunctions\tcoding_junction\n")
            for cluster, data_list in dict_complete_pe_by_cluster.items():
                for junctions, coding_junction in data_list:
                    junction_str = ','.join(list(junctions))
                    coding_str = next(iter(coding_junction))
                    f.write(f"{cluster}\t{junction_str}\t{coding_str}\n")

        with open(introns_file_filtered, 'w') as f:
            f.write("coding_junction\tupstream_intron\tdownstream_intron\n")
            for coding_junction, intron_pairs in dict_PE.items():
                for upstream_intron, downstream_intron in intron_pairs:
                    f.write(f"{coding_junction}\t{upstream_intron}\t{downstream_intron}\n")

        print(f"[Step 2] Complete PE (filtered) saved to: {complete_pe_path_filtered}", file=sys.stderr)
        print(f"[Step 2] Introns (filtered) saved to: {introns_file_filtered}", file=sys.stderr)

    print("[Step 2] Done!\n", file=sys.stderr)

    return introns_file_all, introns_file_filtered


# ============================================================================
# STEP 3: BED_color_convert.py functionality
# ============================================================================

def convert_coords_to_bed(coord_str):
    """Convert chromosome coordinate string to BED format components."""
    chrom, pos = coord_str.split(':')
    start, end = map(int, pos.split('-'))
    return chrom, start, end


def run_step3(args, introns_file_all, introns_file_filtered, cryptic_bed_file, colored_bed_dir):
    """Execute Step 3: BED color conversion for all three file types."""
    print("\n" + "="*80, file=sys.stderr)
    print("STEP 3: BED Color Conversion", file=sys.stderr)
    print("="*80 + "\n", file=sys.stderr)

    sample_name = os.path.basename(args.output_prefix)

    # Colors in RGB format
    #BLUE = "0,0,255"      # for coding junction
    ORANGE = "255,165,0"  # for introns
    TEAL = "153,255,204"  # for cryptic splice sites

    BLUE = "0,0,0"      # really black
    ORANGE = "128,128,128" # really grey

    output_files = {}

    # Process introns file for "all" poison exons
    if introns_file_all and os.path.exists(introns_file_all):
        output_bed_all = os.path.join(colored_bed_dir, f"{sample_name}_colored_junctions_all.bed")
        print(f"[Step 3] Creating colored BED for all poison exons: {output_bed_all}", file=sys.stderr)

        try:
            with open(introns_file_all, 'r') as f, open(output_bed_all, 'w') as out_f:
                header = next(f).strip().split('\t')

                for line in f:
                    fields = line.strip().split('\t')
                    if len(fields) != 3:
                        continue
                    coding_junction, upstream_intron, downstream_intron = fields

                    # Process coding junction (blue)
                    chrom, start, end = convert_coords_to_bed(coding_junction)
                    out_f.write(f"{chrom}\t{start}\t{end}\tcoding_junction\t1000\t.\t{start}\t{end}\t{BLUE}\n")

                    # Process upstream intron (orange)
                    chrom, start, end = convert_coords_to_bed(upstream_intron)
                    out_f.write(f"{chrom}\t{start}\t{end}\tupstream_intron\t1000\t.\t{start}\t{end}\t{ORANGE}\n")

                    # Process downstream intron (orange)
                    chrom, start, end = convert_coords_to_bed(downstream_intron)
                    out_f.write(f"{chrom}\t{start}\t{end}\tdownstream_intron\t1000\t.\t{start}\t{end}\t{ORANGE}\n")

            output_files['all'] = output_bed_all
            print(f"[Step 3] Colored BED (all) saved to: {output_bed_all}", file=sys.stderr)

        except FileNotFoundError:
            print(f"[Step 3] Warning: File {introns_file_all} not found, skipping.", file=sys.stderr)

    # Process introns file for "filtered" poison exons
    if introns_file_filtered and os.path.exists(introns_file_filtered):
        output_bed_filtered = os.path.join(colored_bed_dir, f"{sample_name}_colored_junctions_filtered.bed")
        print(f"[Step 3] Creating colored BED for filtered poison exons: {output_bed_filtered}", file=sys.stderr)

        try:
            with open(introns_file_filtered, 'r') as f, open(output_bed_filtered, 'w') as out_f:
                header = next(f).strip().split('\t')

                for line in f:
                    fields = line.strip().split('\t')
                    if len(fields) != 3:
                        continue
                    coding_junction, upstream_intron, downstream_intron = fields

                    # Process coding junction (blue)
                    chrom, start, end = convert_coords_to_bed(coding_junction)
                    out_f.write(f"{chrom}\t{start}\t{end}\tcoding_junction\t1000\t.\t{start}\t{end}\t{BLUE}\n")

                    # Process upstream intron (orange)
                    chrom, start, end = convert_coords_to_bed(upstream_intron)
                    out_f.write(f"{chrom}\t{start}\t{end}\tupstream_intron\t1000\t.\t{start}\t{end}\t{ORANGE}\n")

                    # Process downstream intron (orange)
                    chrom, start, end = convert_coords_to_bed(downstream_intron)
                    out_f.write(f"{chrom}\t{start}\t{end}\tdownstream_intron\t1000\t.\t{start}\t{end}\t{ORANGE}\n")

            output_files['filtered'] = output_bed_filtered
            print(f"[Step 3] Colored BED (filtered) saved to: {output_bed_filtered}", file=sys.stderr)

        except FileNotFoundError:
            print(f"[Step 3] Warning: File {introns_file_filtered} not found, skipping.", file=sys.stderr)

    # Process cryptic splice sites
    if cryptic_bed_file and os.path.exists(cryptic_bed_file):
        output_bed_cryptic = os.path.join(colored_bed_dir, f"{sample_name}_colored_cryptic_ss.bed")
        print(f"[Step 3] Creating colored BED for cryptic splice sites: {output_bed_cryptic}", file=sys.stderr)

        try:
            with open(cryptic_bed_file, 'r') as f, open(output_bed_cryptic, 'w') as out_f:
                header = next(f).strip().split('\t')

                for line in f:
                    fields = line.strip().split('\t')
                    if len(fields) < 7:
                        continue

                    chrom = fields[0].lstrip('#')
                    start = fields[1]
                    end = fields[2]
                    gene_name = fields[3]
                    strand = fields[4]
                    junction_type = fields[5]
                    cryptic_junction = fields[6]

                    # Write cryptic splice site with teal color
                    out_f.write(f"{chrom}\t{start}\t{end}\t{gene_name}_{junction_type}\t1000\t{strand}\t{start}\t{end}\t{TEAL}\n")

            output_files['cryptic'] = output_bed_cryptic
            print(f"[Step 3] Colored BED (cryptic) saved to: {output_bed_cryptic}", file=sys.stderr)

        except FileNotFoundError:
            print(f"[Step 3] Warning: File {cryptic_bed_file} not found, skipping.", file=sys.stderr)

    print("[Step 3] Done!\n", file=sys.stderr)

    return output_files


# ============================================================================
# Main pipeline execution
# ============================================================================

def main():
    args = parse_args()

    # Setup organized directory structure
    pe_dir, colored_bed_dir = setup_directories(args.output_prefix)

    print("\n" + "="*80, file=sys.stderr)
    print("INTEGRATED POISON EXON DETECTION PIPELINE (Leafcutter2)", file=sys.stderr)
    print("="*80, file=sys.stderr)
    print(f"GTF file: {args.gtf}", file=sys.stderr)
    print(f"Leafcutter2 annotations: {args.leafcutter2_annotations}", file=sys.stderr)
    print(f"Leafcutter2 cluster_ratios: {args.leafcutter_counts}", file=sys.stderr)
    print(f"Output prefix: {args.output_prefix}", file=sys.stderr)
    print(f"PE directory: {pe_dir}", file=sys.stderr)
    print(f"ColoredBed directory: {colored_bed_dir}", file=sys.stderr)
    print(f"Max distance: {args.max_distance}", file=sys.stderr)
    print("="*80 + "\n", file=sys.stderr)

    try:
        # Step 1: Identify poison exons and cryptic splice sites
        all_bed, filtered_bed, cryptic_bed = run_step1(args, pe_dir)

        # Step 2: Complete poison junction analysis for both all and filtered
        introns_file_all, introns_file_filtered = run_step2(args, all_bed, filtered_bed, pe_dir)

        # Step 3: Convert all three file types to colored BED format
        colored_bed_files = run_step3(args, introns_file_all, introns_file_filtered, cryptic_bed, colored_bed_dir)

        # Summary
        print("\n" + "="*80, file=sys.stderr)
        print("PIPELINE COMPLETE", file=sys.stderr)
        print("="*80, file=sys.stderr)
        print("\nOutput files generated:", file=sys.stderr)
        print(f"\nPE Directory ({pe_dir}):", file=sys.stderr)
        print(f"  1. All poison exons: {os.path.basename(all_bed)}", file=sys.stderr)
        print(f"  2. Filtered poison exons: {os.path.basename(filtered_bed)}", file=sys.stderr)
        print(f"  3. Potential cryptic splice sites: {os.path.basename(cryptic_bed)}", file=sys.stderr)
        if introns_file_all:
            print(f"  4. Complete PE clusters (all): {os.path.basename(args.output_prefix)}_complete_PE_all.tsv", file=sys.stderr)
            print(f"  5. Coding/noncoding introns (all): {os.path.basename(introns_file_all)}", file=sys.stderr)
        if introns_file_filtered:
            print(f"  6. Complete PE clusters (filtered): {os.path.basename(args.output_prefix)}_complete_PE_filtered.tsv", file=sys.stderr)
            print(f"  7. Coding/noncoding introns (filtered): {os.path.basename(introns_file_filtered)}", file=sys.stderr)
        print(f"  8. Parsed LeafCutter2 counts: {os.path.basename(args.output_prefix)}_parsed_counts.gz", file=sys.stderr)

        print(f"\nColoredBed Directory ({colored_bed_dir}):", file=sys.stderr)
        if 'all' in colored_bed_files:
            print(f"  9. Colored BED (all): {os.path.basename(colored_bed_files['all'])}", file=sys.stderr)
        if 'filtered' in colored_bed_files:
            print(f"  10. Colored BED (filtered): {os.path.basename(colored_bed_files['filtered'])}", file=sys.stderr)
        if 'cryptic' in colored_bed_files:
            print(f"  11. Colored BED (cryptic): {os.path.basename(colored_bed_files['cryptic'])}", file=sys.stderr)
        print("="*80 + "\n", file=sys.stderr)

    except Exception as e:
        print(f"\n[Pipeline] Error: {str(e)}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
