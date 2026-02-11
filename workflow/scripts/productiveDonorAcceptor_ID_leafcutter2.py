#!/usr/bin/env python3
"""
altDonorAcceptor_ID_leafcutter2.py
Identify Alt Donor and Alt Acceptor Junctions - Leafcutter2 Version

Detects alternative coding junctions where:
- altDonor: 5' donor site is alternative (not at CDS boundary), 3' acceptor site matches CDS boundary
- altAcceptor: 5' donor site matches CDS boundary, 3' acceptor site is alternative (not at CDS boundary)

This identifies alternative splicing events between two PRODUCTIVE (coding) isoforms.

This version is optimized for Leafcutter2 cluster_ratios.gz format.

Author: Dylan Stermer
Date: February 2026
"""

import argparse
import sys
import pandas as pd
import gzip
from collections import defaultdict
import os


def parse_args(args=None):
    parser = argparse.ArgumentParser(
        description="Identify alt donor and acceptor junctions (coding alternative splicing)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('-g', '--gtf', required=True,
                        help='Input GTF file path (e.g., GRCh38 Gencode)')
    parser.add_argument('-l', '--leafcutter2-annotations', required=True,
                        help='Leafcutter2 junction classification file path')
    parser.add_argument('-c', '--leafcutter-counts', required=True,
                        help='Path to Leafcutter2 cluster_ratios file (gzipped)')
    parser.add_argument('-o', '--output-prefix', required=True,
                        help='Output prefix for files (base directory and sample name)')
    parser.add_argument('--chromosome', default=None,
                        help='Filter for a specific chromosome (e.g., chr1) for faster testing')

    return parser.parse_args(args)


def setup_directories(output_prefix):
    """Create output directory structure."""
    base_dir = os.path.dirname(output_prefix)
    if not base_dir:
        base_dir = '.'

    productive_donor_acceptor_dir = os.path.join(base_dir, 'ProductiveDonorAcceptor')
    os.makedirs(productive_donor_acceptor_dir, exist_ok=True)

    return productive_donor_acceptor_dir


def process_gtf_to_dict(gtf_path, chromosome_filter=None):
    """Process GTF file and extract CDS regions."""
    dict_cds = {}

    print(f"[Processing] GTF file: {gtf_path}", file=sys.stderr)
    if chromosome_filter:
        print(f"[Processing] Filtering for chromosome: {chromosome_filter}", file=sys.stderr)
    for line in open(gtf_path):
        if line[0] == "#":
            continue
        chr, source, ann_type, start, stop, dot, strand, frame, info = line.split("\t")
        if ann_type != "CDS":
            continue
        if chromosome_filter and chr != chromosome_filter:
            continue
        gene_name = info.split('gene_name "')[1].split('"')[0]
        if gene_name not in dict_cds:
            dict_cds[gene_name] = set()
        dict_cds[gene_name].add((chr, int(start), int(stop), strand))

    print(f"[Processing] Found CDS regions for {len(dict_cds)} genes", file=sys.stderr)
    return dict_cds


def parse_leafcutter_counts(counts_path, chromosome_filter=None):
    """
    Parse Leafcutter2 cluster_ratios file to get junction-to-cluster mappings.

    Handles Leafcutter2 format: chr:start:end:clu_X_strand:annotation_type
    Coordinates are kept as-is to match the junction classifications format.

    Returns:
        dict_leafcutter_by_junction: {intron_coord: [cluster_ids]}
        dict_leafcutter_by_cluster: {cluster_id: set(intron_coords)}
    """
    dict_leafcutter_by_junction = defaultdict(list)
    dict_leafcutter_by_cluster = defaultdict(set)

    print(f"[Processing] Leafcutter2 cluster_ratios file: {counts_path}", file=sys.stderr)
    if chromosome_filter:
        print(f"[Processing] Filtering for chromosome: {chromosome_filter}", file=sys.stderr)

    opener = gzip.open if counts_path.endswith('.gz') else open

    with opener(counts_path, 'rt') as f:
        # Skip header
        next(f)

        for line in f:
            parts = line.strip().split()
            if not parts:
                continue

            junction_string = parts[0]

            # Parse format: chr:start:end:clu_X_strand or chr:start:end:clu_X_strand:annotation_type
            intron_info = junction_string.split(':')
            if len(intron_info) < 4:  # Allow for extra annotation field in leafcutter2
                continue

            chr_name = intron_info[0]
            if chromosome_filter and chr_name != chromosome_filter:
                continue
            start = intron_info[1]
            stop = intron_info[2]  # Keep original coordinate from leafcutter2
            cluster_info = intron_info[3]  # clu_X_strand
            # intron_info[4] is annotation_type (e.g., "IN", "UP") - ignored

            intron_coord = f"{chr_name}:{start}-{stop}"

            # Extract cluster ID
            cluster_id = cluster_info  # Keep full cluster ID with strand

            dict_leafcutter_by_junction[intron_coord].append(cluster_id)
            dict_leafcutter_by_cluster[cluster_id].add(intron_coord)

    print(f"[Processing] Found {len(dict_leafcutter_by_junction)} junctions in {len(dict_leafcutter_by_cluster)} clusters", file=sys.stderr)
    return dict_leafcutter_by_junction, dict_leafcutter_by_cluster


def process_leafcutter2_coding(leafcutter_path, chromosome_filter=None):
    """
    Process Leafcutter2 annotations for coding, non-UTR junctions.
    Returns:
        dict_coding: {gene_name: set((chr, start, stop, strand, intron_coord))}
        dict_junction_coding_status: {intron_coord: True} for coding junctions
    """
    dict_coding = {}
    dict_junction_coding_status = {}

    print(f"[Processing] Leafcutter2 coding junctions: {leafcutter_path}", file=sys.stderr)
    if chromosome_filter:
        print(f"[Processing] Filtering for chromosome: {chromosome_filter}", file=sys.stderr)
    for line in open(leafcutter_path):
        if line.strip().startswith("Gene_name"):
            continue
        gene_name, intron_cord, strand, annotation, coding, UTR, GencodePC = line.strip().split("\t")

        chrom, positions = intron_cord.split(':')
        if chromosome_filter and chrom != chromosome_filter:
            continue
        start, stop = positions.split('-')

        # Store coding status for all junctions
        if coding == "True":
            dict_junction_coding_status[intron_cord] = True

        # Filter for coding, non-UTR junctions
        if coding == "True" and UTR == "False":
            if gene_name not in dict_coding:
                dict_coding[gene_name] = set()

            # Apply +1 adjustment to stop coordinate for both strands
            stop = int(stop) + 1
            dict_coding[gene_name].add((chrom, int(start), stop, strand, intron_cord))

    print(f"[Processing] Found coding junctions for {len(dict_coding)} genes", file=sys.stderr)
    print(f"[Processing] Coding status recorded for {len(dict_junction_coding_status)} junctions", file=sys.stderr)
    return dict_coding, dict_junction_coding_status


def detect_productive_donors_acceptors(dict_cds, dict_coding,
                                        dict_leafcutter_by_junction, dict_junction_coding_status):
    """
    Detect alt donor and acceptor junctions.

    LOGIC (coordinate-based, works for both strands):
    - altDonor: coding junction where acceptor matches CDS boundary, donor does NOT
    - altAcceptor: coding junction where donor matches CDS boundary, acceptor does NOT

    For each productive junction, find other CODING junctions that share the canonical splice site.

    COORDINATE-BASED:
    - On + strand: start=donor, stop=acceptor
    - On - strand: stop=donor, start=acceptor

    VALIDATION:
    - Junction must exist in LeafCutter counts
    - Junction must be coding
    """
    productive_donors = []
    productive_acceptors = []

    print("[Processing] Identifying alt donor and acceptor junctions...", file=sys.stderr)

    for gene_name in dict_coding:
        if gene_name not in dict_cds:
            continue

        coding_exons = dict_cds[gene_name]
        coding_junctions = dict_coding[gene_name]

        # Build sets of CDS boundaries
        cds_starts = set()
        cds_stops = set()
        for cds_chr, cds_start, cds_stop, cds_strand in coding_exons:
            cds_starts.add((cds_chr, cds_start))
            cds_stops.add((cds_chr, cds_stop))

        for c_chr, c_start, c_stop, c_strand, c_coord in coding_junctions:
            # VALIDATION: Junction must exist in LeafCutter counts
            if c_coord not in dict_leafcutter_by_junction:
                continue

            # Get CDS boundaries for this chromosome
            chr_cds_starts = {start for chr, start in cds_starts if chr == c_chr}
            chr_cds_stops = {stop for chr, stop in cds_stops if chr == c_chr}

            # COORDINATE-BASED DETECTION
            # Case 1: c_stop matches CDS start (acceptor at CDS boundary)
            matches_cds_start_at_stop = c_stop in chr_cds_starts
            not_matches_cds_stop_at_start = c_start not in chr_cds_stops

            # Case 2: c_start matches CDS stop (donor at CDS boundary)
            matches_cds_stop_at_start = c_start in chr_cds_stops
            not_matches_cds_start_at_stop = c_stop not in chr_cds_starts

            # STRAND-DEPENDENT LABELING
            if c_strand == "+":
                # Positive strand: start=donor(5'), stop=acceptor(3')

                # altDonor: acceptor matches CDS, donor is alternative
                if matches_cds_start_at_stop and not_matches_cds_stop_at_start:
                    junction_type = 'altDonor'
                    donor_coord = c_start
                    acceptor_coord = c_stop

                    # Find ALL other coding junctions sharing acceptor (c_stop)
                    spanning_junctions = []
                    for other_chr, other_start, other_stop, other_strand, other_coord in coding_junctions:
                        if other_chr == c_chr and other_strand == c_strand:
                            # On + strand, acceptor is at stop position
                            if other_stop == c_stop and other_coord != c_coord:
                                if other_coord in dict_junction_coding_status:
                                    spanning_junctions.append(other_coord)

                    productive_donors.append({
                        'chr': c_chr,
                        'start': c_start,
                        'stop': c_stop,
                        'gene_name': gene_name,
                        'strand': c_strand,
                        'junction_type': junction_type,
                        'intron_coord': c_coord,
                        'spanning_productive': ','.join(spanning_junctions) if spanning_junctions else 'None',
                        'donor_coord': donor_coord,
                        'acceptor_coord': acceptor_coord
                    })

                # altAcceptor: donor matches CDS, acceptor is alternative
                elif matches_cds_stop_at_start and not_matches_cds_start_at_stop:
                    junction_type = 'altAcceptor'
                    donor_coord = c_start
                    acceptor_coord = c_stop

                    # Find ALL other coding junctions sharing donor (c_start)
                    spanning_junctions = []
                    for other_chr, other_start, other_stop, other_strand, other_coord in coding_junctions:
                        if other_chr == c_chr and other_strand == c_strand:
                            # On + strand, donor is at start position
                            if other_start == c_start and other_coord != c_coord:
                                if other_coord in dict_junction_coding_status:
                                    spanning_junctions.append(other_coord)

                    productive_acceptors.append({
                        'chr': c_chr,
                        'start': c_start,
                        'stop': c_stop,
                        'gene_name': gene_name,
                        'strand': c_strand,
                        'junction_type': junction_type,
                        'intron_coord': c_coord,
                        'spanning_productive': ','.join(spanning_junctions) if spanning_junctions else 'None',
                        'donor_coord': donor_coord,
                        'acceptor_coord': acceptor_coord
                    })

            else:  # Negative strand
                # Negative strand: stop=donor(5'), start=acceptor(3')

                # altAcceptor: donor (c_stop) matches CDS start, acceptor is alternative
                if matches_cds_start_at_stop and not_matches_cds_stop_at_start:
                    junction_type = 'altAcceptor'
                    donor_coord = c_stop
                    acceptor_coord = c_start

                    # Find ALL other coding junctions sharing donor (c_stop on - strand)
                    spanning_junctions = []
                    for other_chr, other_start, other_stop, other_strand, other_coord in coding_junctions:
                        if other_chr == c_chr and other_strand == c_strand:
                            # On - strand, donor is at stop position
                            if other_stop == c_stop and other_coord != c_coord:
                                if other_coord in dict_junction_coding_status:
                                    spanning_junctions.append(other_coord)

                    productive_acceptors.append({
                        'chr': c_chr,
                        'start': c_start,
                        'stop': c_stop,
                        'gene_name': gene_name,
                        'strand': c_strand,
                        'junction_type': junction_type,
                        'intron_coord': c_coord,
                        'spanning_productive': ','.join(spanning_junctions) if spanning_junctions else 'None',
                        'donor_coord': donor_coord,
                        'acceptor_coord': acceptor_coord
                    })

                # altDonor: acceptor (c_start) matches CDS stop, donor is alternative
                elif matches_cds_stop_at_start and not_matches_cds_start_at_stop:
                    junction_type = 'altDonor'
                    donor_coord = c_stop
                    acceptor_coord = c_start

                    # Find ALL other coding junctions sharing acceptor (c_start on - strand)
                    spanning_junctions = []
                    for other_chr, other_start, other_stop, other_strand, other_coord in coding_junctions:
                        if other_chr == c_chr and other_strand == c_strand:
                            # On - strand, acceptor is at start position
                            if other_start == c_start and other_coord != c_coord:
                                if other_coord in dict_junction_coding_status:
                                    spanning_junctions.append(other_coord)

                    productive_donors.append({
                        'chr': c_chr,
                        'start': c_start,
                        'stop': c_stop,
                        'gene_name': gene_name,
                        'strand': c_strand,
                        'junction_type': junction_type,
                        'intron_coord': c_coord,
                        'spanning_productive': ','.join(spanning_junctions) if spanning_junctions else 'None',
                        'donor_coord': donor_coord,
                        'acceptor_coord': acceptor_coord
                    })

    print(f"[Processing] Found {len(productive_donors)} altDonors", file=sys.stderr)
    print(f"[Processing] Found {len(productive_acceptors)} altAcceptors", file=sys.stderr)

    return productive_donors, productive_acceptors


def find_spanning_productive_by_boundary(productive_donors, productive_acceptors, dict_coding, dict_junction_coding_status):
    """
    Post-processing: For each alt junction, find ALL spanning productive junctions that share a splice site boundary.
    This supplements the inline matching done during detection.
    """
    print("[Processing] Post-processing: finding spanning productive junctions by shared boundaries...", file=sys.stderr)

    # Build a lookup by coordinate for all coding junctions
    coding_by_start = defaultdict(list)  # {(chr, start, strand): [junction_coords]}
    coding_by_stop = defaultdict(list)   # {(chr, stop, strand): [junction_coords]}

    for gene_name, junctions in dict_coding.items():
        for chr, start, stop, strand, junction_coord in junctions:
            # Only include if it's truly coding
            if junction_coord in dict_junction_coding_status:
                coding_by_start[(chr, start, strand)].append(junction_coord)
                coding_by_stop[(chr, stop, strand)].append(junction_coord)

    donors_updated = 0
    acceptors_updated = 0

    # Update altDonors with spanning productive junctions
    for entry in productive_donors:
        chr = entry['chr']
        strand = entry['strand']
        acceptor_coord = entry['acceptor_coord']
        current_junction = entry['intron_coord']

        # For altDonor: donor is alternative, so find junctions sharing the acceptor
        if strand == "+":
            # Acceptor is at stop position on + strand
            associated = coding_by_stop.get((chr, acceptor_coord, strand), [])
        else:
            # Acceptor is at start position on - strand
            associated = coding_by_start.get((chr, acceptor_coord, strand), [])

        # Filter out the current junction
        associated = [j for j in associated if j != current_junction]

        if associated:
            entry['spanning_productive'] = ','.join(associated)
            donors_updated += 1

    # Update altAcceptors with spanning productive junctions
    for entry in productive_acceptors:
        chr = entry['chr']
        strand = entry['strand']
        donor_coord = entry['donor_coord']
        current_junction = entry['intron_coord']

        # For altAcceptor: acceptor is alternative, so find junctions sharing the donor
        if strand == "+":
            # Donor is at start position on + strand
            associated = coding_by_start.get((chr, donor_coord, strand), [])
        else:
            # Donor is at stop position on - strand
            associated = coding_by_stop.get((chr, donor_coord, strand), [])

        # Filter out the current junction
        associated = [j for j in associated if j != current_junction]

        if associated:
            entry['spanning_productive'] = ','.join(associated)
            acceptors_updated += 1

    print(f"[Processing] Updated {donors_updated} altDonors and {acceptors_updated} altAcceptors with spanning productive junctions", file=sys.stderr)

    return productive_donors, productive_acceptors


def calculate_frame_differences(productive_donors, productive_acceptors):
    """
    Calculate the difference in nucleotides between alt junctions and their spanning productive junctions.
    Check if the difference is a multiple of 3 (frame-preserving).

    For altDonor: difference is in donor position (acceptor is shared)
    For altAcceptor: difference is in acceptor position (donor is shared)
    """
    print("[Processing] Calculating frame differences...", file=sys.stderr)

    for entry in productive_donors:
        if entry['spanning_productive'] == 'None':
            entry['frame_diffs'] = 'None'
            entry['frame_preserved'] = 'None'
            continue

        strand = entry['strand']
        # For altDonor, the donor is the alternative site
        if strand == "+":
            alt_donor = entry['start']  # donor is at start on + strand
        else:
            alt_donor = entry['stop']   # donor is at stop on - strand

        frame_diffs = []
        frame_preserved = []

        for spanning_coord in entry['spanning_productive'].split(','):
            span_chr, span_positions = spanning_coord.split(':')
            span_start, span_stop = map(int, span_positions.split('-'))

            if strand == "+":
                span_donor = span_start
            else:
                span_donor = span_stop + 1  # +1 to match the adjustment

            diff = abs(alt_donor - span_donor)
            frame_diffs.append(str(diff))
            frame_preserved.append('Y' if diff % 3 == 0 else 'N')

        entry['frame_diffs'] = ','.join(frame_diffs)
        entry['frame_preserved'] = ','.join(frame_preserved)

    for entry in productive_acceptors:
        if entry['spanning_productive'] == 'None':
            entry['frame_diffs'] = 'None'
            entry['frame_preserved'] = 'None'
            continue

        strand = entry['strand']
        # For altAcceptor, the acceptor is the alternative site
        if strand == "+":
            alt_acceptor = entry['stop']   # acceptor is at stop on + strand
        else:
            alt_acceptor = entry['start']  # acceptor is at start on - strand

        frame_diffs = []
        frame_preserved = []

        for spanning_coord in entry['spanning_productive'].split(','):
            span_chr, span_positions = spanning_coord.split(':')
            span_start, span_stop = map(int, span_positions.split('-'))

            if strand == "+":
                span_acceptor = span_stop + 1  # +1 to match the adjustment
            else:
                span_acceptor = span_start

            diff = abs(alt_acceptor - span_acceptor)
            frame_diffs.append(str(diff))
            frame_preserved.append('Y' if diff % 3 == 0 else 'N')

        entry['frame_diffs'] = ','.join(frame_diffs)
        entry['frame_preserved'] = ','.join(frame_preserved)

    # Calculate statistics
    donor_with_spanning = [d for d in productive_donors if d['spanning_productive'] != 'None']
    acceptor_with_spanning = [a for a in productive_acceptors if a['spanning_productive'] != 'None']

    # Count frame-preserved pairs
    donor_frame_preserved = sum(1 for d in donor_with_spanning if 'Y' in d['frame_preserved'])
    acceptor_frame_preserved = sum(1 for a in acceptor_with_spanning if 'Y' in a['frame_preserved'])

    print(f"[Processing] altDonors with at least one frame-preserving spanning: {donor_frame_preserved}/{len(donor_with_spanning)}", file=sys.stderr)
    print(f"[Processing] altAcceptors with at least one frame-preserving spanning: {acceptor_frame_preserved}/{len(acceptor_with_spanning)}", file=sys.stderr)

    return productive_donors, productive_acceptors


def write_bed_output(productive_donors, productive_acceptors, output_prefix, output_dir):
    """Write BED format output file."""
    sample_name = os.path.basename(output_prefix)
    output_bed = os.path.join(output_dir, f"{sample_name}_altDonorAcceptor.bed")

    print(f"[Output] Writing BED file: {output_bed}", file=sys.stderr)

    with open(output_bed, 'w') as f:
        # Write header
        f.write("#chr\tstart\tend\tname\tstrand\tjunction_type\tintron_coord\tspanning_productive\tframe_diffs\tframe_preserved\tdonor_coord\tacceptor_coord\n")

        # Write altDonors
        for entry in productive_donors:
            # Subtract 1 from stop coordinate (reversing the +1 from parsing)
            stop_coord = entry['stop'] - 1
            f.write(f"{entry['chr']}\t{entry['start']}\t{stop_coord}\t"
                   f"{entry['gene_name']}\t{entry['strand']}\t{entry['junction_type']}\t"
                   f"{entry['intron_coord']}\t{entry['spanning_productive']}\t"
                   f"{entry['frame_diffs']}\t{entry['frame_preserved']}\t"
                   f"{entry['donor_coord']}\t{entry['acceptor_coord']}\n")

        # Write altAcceptors
        for entry in productive_acceptors:
            # Subtract 1 from stop coordinate (reversing the +1 from parsing)
            stop_coord = entry['stop'] - 1
            f.write(f"{entry['chr']}\t{entry['start']}\t{stop_coord}\t"
                   f"{entry['gene_name']}\t{entry['strand']}\t{entry['junction_type']}\t"
                   f"{entry['intron_coord']}\t{entry['spanning_productive']}\t"
                   f"{entry['frame_diffs']}\t{entry['frame_preserved']}\t"
                   f"{entry['donor_coord']}\t{entry['acceptor_coord']}\n")

    print(f"[Output] BED file written successfully", file=sys.stderr)
    return output_bed


def write_colored_bed(productive_donors, productive_acceptors, output_prefix, output_dir):
    """Write colored BED9 format for visualization."""
    sample_name = os.path.basename(output_prefix)
    output_bed = os.path.join(output_dir, f"{sample_name}_colored_altDonorAcceptor.bed")

    print(f"[Output] Writing colored BED file: {output_bed}", file=sys.stderr)

    # Define colors (different from poison to distinguish)
    BLUE = "0,0,255"           # altDonor
    CYAN = "0,200,200"         # altAcceptor
    ORANGE = "255,165,0"       # associated junctions

    with open(output_bed, 'w') as f:
        # Write altDonors in blue
        for entry in productive_donors:
            # Subtract 1 from stop coordinate (reversing the +1 from parsing)
            stop_coord = entry['stop'] - 1
            f.write(f"{entry['chr']}\t{entry['start']}\t{stop_coord}\t"
                   f"{entry['gene_name']}_{entry['junction_type']}\t1000\t{entry['strand']}\t"
                   f"{entry['start']}\t{stop_coord}\t{BLUE}\n")

            # Write associated coding junctions in orange
            if entry['spanning_productive'] != 'None':
                for assoc_junction in entry['spanning_productive'].split(','):
                    assoc_chr, assoc_positions = assoc_junction.split(':')
                    assoc_start, assoc_stop = assoc_positions.split('-')
                    f.write(f"{assoc_chr}\t{assoc_start}\t{assoc_stop}\t"
                           f"{entry['gene_name']}_spanning_productive\t1000\t{entry['strand']}\t"
                           f"{assoc_start}\t{assoc_stop}\t{ORANGE}\n")

        # Write altAcceptors in cyan
        for entry in productive_acceptors:
            # Subtract 1 from stop coordinate (reversing the +1 from parsing)
            stop_coord = entry['stop'] - 1
            f.write(f"{entry['chr']}\t{entry['start']}\t{stop_coord}\t"
                   f"{entry['gene_name']}_{entry['junction_type']}\t1000\t{entry['strand']}\t"
                   f"{entry['start']}\t{stop_coord}\t{CYAN}\n")

            # Write associated coding junctions in orange
            if entry['spanning_productive'] != 'None':
                for assoc_junction in entry['spanning_productive'].split(','):
                    assoc_chr, assoc_positions = assoc_junction.split(':')
                    assoc_start, assoc_stop = assoc_positions.split('-')
                    f.write(f"{assoc_chr}\t{assoc_start}\t{assoc_stop}\t"
                           f"{entry['gene_name']}_spanning_productive\t1000\t{entry['strand']}\t"
                           f"{assoc_start}\t{assoc_stop}\t{ORANGE}\n")

    print(f"[Output] Colored BED file written successfully", file=sys.stderr)
    return output_bed


def write_summary_stats(productive_donors, productive_acceptors, output_prefix, output_dir):
    """Write summary statistics."""
    sample_name = os.path.basename(output_prefix)
    output_stats = os.path.join(output_dir, f"{sample_name}_altDonorAcceptor_stats.txt")

    print(f"[Output] Writing summary statistics: {output_stats}", file=sys.stderr)

    # Calculate statistics
    total_donors = len(productive_donors)
    total_acceptors = len(productive_acceptors)

    donors_with_associated = sum(1 for d in productive_donors if d['spanning_productive'] != 'None')
    acceptors_with_associated = sum(1 for a in productive_acceptors if a['spanning_productive'] != 'None')

    # Count associated junctions per productive junction
    donor_assoc_counts = [len(d['spanning_productive'].split(',')) for d in productive_donors if d['spanning_productive'] != 'None']
    acceptor_assoc_counts = [len(a['spanning_productive'].split(',')) for a in productive_acceptors if a['spanning_productive'] != 'None']

    # Frame preservation statistics
    donors_with_spanning = [d for d in productive_donors if d['spanning_productive'] != 'None']
    acceptors_with_spanning = [a for a in productive_acceptors if a['spanning_productive'] != 'None']

    # Count junctions with at least one frame-preserving pair
    donor_any_frame_preserved = sum(1 for d in donors_with_spanning if 'Y' in d.get('frame_preserved', ''))
    acceptor_any_frame_preserved = sum(1 for a in acceptors_with_spanning if 'Y' in a.get('frame_preserved', ''))

    # Count all individual pairs
    donor_total_pairs = 0
    donor_frame_preserved_pairs = 0
    for d in donors_with_spanning:
        if d.get('frame_preserved') and d['frame_preserved'] != 'None':
            pairs = d['frame_preserved'].split(',')
            donor_total_pairs += len(pairs)
            donor_frame_preserved_pairs += pairs.count('Y')

    acceptor_total_pairs = 0
    acceptor_frame_preserved_pairs = 0
    for a in acceptors_with_spanning:
        if a.get('frame_preserved') and a['frame_preserved'] != 'None':
            pairs = a['frame_preserved'].split(',')
            acceptor_total_pairs += len(pairs)
            acceptor_frame_preserved_pairs += pairs.count('Y')

    with open(output_stats, 'w') as f:
        f.write("Alt Donor/Acceptor Junction Detection Summary\n")
        f.write("=" * 60 + "\n\n")

        f.write(f"Total altDonors identified: {total_donors}\n")
        f.write(f"Total altAcceptors identified: {total_acceptors}\n")
        f.write(f"Total alt donor/acceptor junctions: {total_donors + total_acceptors}\n\n")

        f.write(f"altDonors with spanning productive junctions: {donors_with_associated}\n")
        f.write(f"altAcceptors with spanning productive junctions: {acceptors_with_associated}\n\n")

        if donor_assoc_counts:
            f.write(f"Spanning productive per altDonor (mean): {sum(donor_assoc_counts)/len(donor_assoc_counts):.2f}\n")
            f.write(f"Spanning productive per altDonor (max): {max(donor_assoc_counts)}\n\n")

        if acceptor_assoc_counts:
            f.write(f"Spanning productive per altAcceptor (mean): {sum(acceptor_assoc_counts)/len(acceptor_assoc_counts):.2f}\n")
            f.write(f"Spanning productive per altAcceptor (max): {max(acceptor_assoc_counts)}\n\n")

        f.write("=" * 60 + "\n")
        f.write("FRAME PRESERVATION ANALYSIS\n")
        f.write("=" * 60 + "\n\n")

        f.write(f"altDonors with >= 1 frame-preserving spanning: {donor_any_frame_preserved}/{len(donors_with_spanning)}")
        if len(donors_with_spanning) > 0:
            f.write(f" ({100*donor_any_frame_preserved/len(donors_with_spanning):.1f}%)\n")
        else:
            f.write("\n")

        f.write(f"altAcceptors with >= 1 frame-preserving spanning: {acceptor_any_frame_preserved}/{len(acceptors_with_spanning)}")
        if len(acceptors_with_spanning) > 0:
            f.write(f" ({100*acceptor_any_frame_preserved/len(acceptors_with_spanning):.1f}%)\n\n")
        else:
            f.write("\n\n")

        f.write(f"Total altDonor-spanning pairs: {donor_total_pairs}\n")
        f.write(f"  Frame-preserving (diff % 3 == 0): {donor_frame_preserved_pairs}")
        if donor_total_pairs > 0:
            f.write(f" ({100*donor_frame_preserved_pairs/donor_total_pairs:.1f}%)\n")
        else:
            f.write("\n")
        f.write(f"  Frame-shifting (diff % 3 != 0): {donor_total_pairs - donor_frame_preserved_pairs}")
        if donor_total_pairs > 0:
            f.write(f" ({100*(donor_total_pairs - donor_frame_preserved_pairs)/donor_total_pairs:.1f}%)\n\n")
        else:
            f.write("\n\n")

        f.write(f"Total altAcceptor-spanning pairs: {acceptor_total_pairs}\n")
        f.write(f"  Frame-preserving (diff % 3 == 0): {acceptor_frame_preserved_pairs}")
        if acceptor_total_pairs > 0:
            f.write(f" ({100*acceptor_frame_preserved_pairs/acceptor_total_pairs:.1f}%)\n")
        else:
            f.write("\n")
        f.write(f"  Frame-shifting (diff % 3 != 0): {acceptor_total_pairs - acceptor_frame_preserved_pairs}")
        if acceptor_total_pairs > 0:
            f.write(f" ({100*(acceptor_total_pairs - acceptor_frame_preserved_pairs)/acceptor_total_pairs:.1f}%)\n")
        else:
            f.write("\n")

    print(f"[Output] Summary statistics written successfully", file=sys.stderr)
    return output_stats


def main():
    args = parse_args()

    # Setup directories
    output_dir = setup_directories(args.output_prefix)

    print("\n" + "="*80, file=sys.stderr)
    print("ALT DONOR/ACCEPTOR JUNCTION DETECTION", file=sys.stderr)
    print("="*80, file=sys.stderr)
    print(f"GTF file: {args.gtf}", file=sys.stderr)
    print(f"Leafcutter2 annotations: {args.leafcutter2_annotations}", file=sys.stderr)
    print(f"Leafcutter2 cluster_ratios: {args.leafcutter_counts}", file=sys.stderr)
    print(f"Output prefix: {args.output_prefix}", file=sys.stderr)
    print(f"Output directory: {output_dir}", file=sys.stderr)
    if args.chromosome:
        print(f"Chromosome filter: {args.chromosome}", file=sys.stderr)
    print("="*80 + "\n", file=sys.stderr)

    try:
        # Step 1: Process GTF for CDS regions
        dict_cds = process_gtf_to_dict(args.gtf, args.chromosome)

        # Step 2: Parse LeafCutter counts file
        dict_leafcutter_by_junction, dict_leafcutter_by_cluster = parse_leafcutter_counts(args.leafcutter_counts, args.chromosome)

        # Step 3: Process Leafcutter2 for coding junctions
        dict_coding, dict_junction_coding_status = process_leafcutter2_coding(args.leafcutter2_annotations, args.chromosome)

        # Step 4: Detect alt donors and acceptors
        productive_donors, productive_acceptors = detect_productive_donors_acceptors(
            dict_cds, dict_coding,
            dict_leafcutter_by_junction, dict_junction_coding_status
        )

        # Step 5: Post-process to find all associated coding junctions sharing boundaries
        productive_donors, productive_acceptors = find_spanning_productive_by_boundary(
            productive_donors, productive_acceptors, dict_coding, dict_junction_coding_status
        )

        # Step 6: Calculate frame differences
        productive_donors, productive_acceptors = calculate_frame_differences(
            productive_donors, productive_acceptors
        )

        # Step 7: Write outputs
        bed_file = write_bed_output(productive_donors, productive_acceptors, args.output_prefix, output_dir)
        colored_bed_file = write_colored_bed(productive_donors, productive_acceptors, args.output_prefix, output_dir)
        stats_file = write_summary_stats(productive_donors, productive_acceptors, args.output_prefix, output_dir)

        # Summary
        print("\n" + "="*80, file=sys.stderr)
        print("DETECTION COMPLETE", file=sys.stderr)
        print("="*80, file=sys.stderr)
        print("\nOutput files generated:", file=sys.stderr)
        print(f"  1. BED file: {os.path.basename(bed_file)}", file=sys.stderr)
        print(f"  2. Colored BED file: {os.path.basename(colored_bed_file)}", file=sys.stderr)
        print(f"  3. Summary statistics: {os.path.basename(stats_file)}", file=sys.stderr)
        print("="*80 + "\n", file=sys.stderr)

    except Exception as e:
        print(f"\n[Error] {str(e)}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
