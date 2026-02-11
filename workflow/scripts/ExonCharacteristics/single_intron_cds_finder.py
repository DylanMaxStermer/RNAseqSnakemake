#!/usr/bin/env python3
"""
Single-Intron CDS Finder - Identify CDS exons adjacent to single-intron Leafcutter clusters

This script identifies single-intron Leafcutter clusters and determines which CDS exons
are adjacent to those introns using GTF annotation. It correctly handles coordinate
conversion between BED-style (0-based, half-open) and GTF (1-based, inclusive) formats.

Author: Dylan Stermer
Date: 2026-02-04
"""

import argparse
import sys
import os
from collections import defaultdict
import re


def parse_leafcutter_clusters(leafcutter_file, min_counts, chromosome_filter='All'):
    """
    Parse Leafcutter2 clustering output and extract single-junction clusters.

    Leafcutter format: Chromosome:Strand Junction1 Junction2 ...
    Junction format: start:end:count (BED-style: 0-based, half-open)

    Args:
        leafcutter_file: Path to Leafcutter2 clustering output
        min_counts: Minimum junction read count threshold
        chromosome_filter: Specific chromosome to analyze (default: 'All')

    Returns:
        Dictionary: {intron_coord: (chr, start, end, strand, counts)}
        where intron_coord is formatted as "chr:start-end"
    """
    introns = {}
    total_clusters = 0
    single_junction_clusters = 0
    filtered_by_counts = 0
    filtered_by_chr = 0

    print(f"\nParsing Leafcutter clusters from: {leafcutter_file}")
    print(f"  Minimum read count threshold: {min_counts}")
    if chromosome_filter != 'All':
        print(f"  Filtering to chromosome: {chromosome_filter}")

    with open(leafcutter_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            total_clusters += 1

            # Parse format: Chromosome:Strand Junction1 Junction2 ...
            # Example: chr1:+ 3365091:3436686:9 3436566:3436686:7
            tokens = line.split()
            if len(tokens) < 1:
                print(f"Warning: Empty line or malformed: {line[:100]}", file=sys.stderr)
                continue

            chr_strand = tokens[0]
            junctions = tokens[1:] if len(tokens) > 1 else []

            # Parse chromosome and strand
            if ':' not in chr_strand:
                print(f"Warning: Cannot parse chr:strand: {chr_strand}", file=sys.stderr)
                continue

            chr_name, strand = chr_strand.rsplit(':', 1)

            # Filter by chromosome if specified
            if chromosome_filter != 'All' and chr_name != chromosome_filter:
                filtered_by_chr += 1
                continue

            # Only keep single-junction clusters
            if len(junctions) == 1:
                single_junction_clusters += 1

                # Parse junction: start:end:count
                junction_parts = junctions[0].split(':')
                if len(junction_parts) != 3:
                    print(f"Warning: Malformed junction: {junctions[0]}", file=sys.stderr)
                    continue

                try:
                    start = int(junction_parts[0])
                    end = int(junction_parts[1])
                    counts = int(junction_parts[2])
                except ValueError:
                    print(f"Warning: Non-integer values in junction: {junctions[0]}", file=sys.stderr)
                    continue

                # Filter by read count
                if counts >= min_counts:
                    intron_coord = f"{chr_name}:{start}-{end}"
                    introns[intron_coord] = (chr_name, start, end, strand, counts)
                else:
                    filtered_by_counts += 1

    print(f"\n  Cluster parsing summary:")
    print(f"    Total clusters: {total_clusters}")
    print(f"    Single-junction clusters: {single_junction_clusters}")
    if chromosome_filter != 'All':
        print(f"    Filtered by chromosome: {filtered_by_chr}")
    print(f"    Filtered by read count (<{min_counts}): {filtered_by_counts}")
    print(f"    Final introns retained: {len(introns)}")

    return introns


def parse_gtf_attribute(attr_string, key):
    """
    Extract a specific attribute value from GTF attribute string.

    Args:
        attr_string: GTF attributes (9th column)
        key: Attribute key to extract (e.g., 'gene_name')

    Returns:
        Attribute value or None if not found
    """
    # Match pattern: key "value"; or key "value"
    pattern = f'{key} "([^"]+)"'
    match = re.search(pattern, attr_string)
    if match:
        return match.group(1)
    return None


def parse_gtf_cds(gtf_file, chromosome_filter='All'):
    """
    Parse GTF file and return CDS features with transcript structure.

    Returns both individual CDS and transcript-level first/last exon info.
    """
    from collections import defaultdict

    cds_count = 0
    filtered_count = 0
    transcript_cds = defaultdict(list)  # transcript_id -> list of (chr, start, end, strand, gene_name, gene_id)

    print(f"\nParsing GTF file: {gtf_file}")
    if chromosome_filter != 'All':
        print(f"  Filtering to chromosome: {chromosome_filter}")

    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            line = line.strip()
            if not line:
                continue

            fields = line.split('\t')
            if len(fields) < 9:
                continue

            if fields[2] != 'CDS':
                continue

            chr_name = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attributes = fields[8]

            if chromosome_filter != 'All' and chr_name != chromosome_filter:
                filtered_count += 1
                continue

            cds_count += 1
            if cds_count % 50000 == 0:
                print(f"  Processing CDS: {cds_count}...", end='\r')

            gene_name = parse_gtf_attribute(attributes, 'gene_name')
            transcript_id = parse_gtf_attribute(attributes, 'transcript_id')
            gene_id = parse_gtf_attribute(attributes, 'gene_id')

            transcript_cds[transcript_id].append((chr_name, start, end, strand, gene_name, gene_id))

    print(f"  Total CDS processed: {cds_count}")

    # Determine first/last CDS per transcript
    first_last_cds = {}  # (chr, start, end, strand) -> set of ('first', 'last')

    for transcript_id, cds_list in transcript_cds.items():
        if not cds_list:
            continue

        strand = cds_list[0][3]

        if strand == '+':
            first_cds = min(cds_list, key=lambda x: x[1])  # Lowest start
            last_cds = max(cds_list, key=lambda x: x[2])   # Highest end
        else:  # strand == '-'
            first_cds = max(cds_list, key=lambda x: x[2])  # Highest end (5' on - strand)
            last_cds = min(cds_list, key=lambda x: x[1])   # Lowest start (3' on - strand)

        first_key = (first_cds[0], first_cds[1], first_cds[2], first_cds[3])
        last_key = (last_cds[0], last_cds[1], last_cds[2], last_cds[3])

        if first_key not in first_last_cds:
            first_last_cds[first_key] = set()
        first_last_cds[first_key].add('first')

        if last_key not in first_last_cds:
            first_last_cds[last_key] = set()
        first_last_cds[last_key].add('last')

    # Yield all CDS with terminal info
    for transcript_id, cds_list in transcript_cds.items():
        for cds in cds_list:
            chr_name, start, end, strand, gene_name, gene_id = cds
            cds_key = (chr_name, start, end, strand)
            terminal_status = first_last_cds.get(cds_key, set())
            yield (chr_name, start, end, strand, gene_name, transcript_id, gene_id, terminal_status)


def match_cds_to_introns(gtf_file, introns, chromosome_filter='All'):
    """
    Match CDS boundaries to intron boundaries with correct coordinate conversion.

    COORDINATE SYSTEMS:
    - Leafcutter introns: 0-based, half-open [start, end)
      Example: [100, 200) means bases at indices 100-199 (0-based)
               which are bases 101-200 in 1-based coordinates

    - GTF CDS: 1-based, inclusive [start, end]
      Example: [95, 100] means bases 95-100 (1-based)
               which are bases at indices 94-99 (0-based)

    MATCHING LOGIC:

    On + strand (5'→3' left to right):
    - 5' splice site (donor): CDS ends, intron begins
      → CDS_end (1-based) == intron_start (0-based)
      Example: CDS ends at 100 (1-based), intron starts at [100, ...) (0-based)

    - 3' splice site (acceptor): intron ends, CDS begins
      → CDS_start (1-based) == intron_end (0-based) + 1
      Example: Intron ends at [..., 200) (0-based), CDS starts at 201 (1-based)

    On - strand (5'→3' right to left, so 5' and 3' are flipped genomically):
    - 5' splice site: genomically downstream (higher coord), CDS starts after intron
      → CDS_start (1-based) == intron_end (0-based) + 1

    - 3' splice site: genomically upstream (lower coord), CDS ends before intron
      → CDS_end (1-based) == intron_start (0-based)

    Args:
        gtf_file: Path to GTF file
        introns: Dictionary of introns from parse_leafcutter_clusters()
        chromosome_filter: Specific chromosome to analyze (default: 'All')

    Returns:
        Tuple: (five_prime_dict, three_prime_dict)
        Each dict maps CDS_ID (by coordinates) to:
        {
            'cds_info': {chr, start, end, strand, gene_names, transcript_ids},
            'introns': [(intron_coord, counts), ...]
        }
    """
    # Use coordinate-based CDS ID to deduplicate
    five_prime_matches = defaultdict(lambda: {'cds_info': None, 'introns': set(), 'transcripts': set(), 'genes': set(), 'gene_ids': set(), 'terminal_status': set()})
    three_prime_matches = defaultdict(lambda: {'cds_info': None, 'introns': set(), 'transcripts': set(), 'genes': set(), 'gene_ids': set(), 'terminal_status': set()})

    # Organize introns by chromosome and strand for faster lookup
    introns_by_chr_strand = defaultdict(list)
    for intron_coord, (chr_name, start, end, strand, counts) in introns.items():
        introns_by_chr_strand[(chr_name, strand)].append((intron_coord, start, end, counts))

    print("\nMatching CDS to intron boundaries...")
    cds_processed = 0

    # Stream through GTF file
    for cds in parse_gtf_cds(gtf_file, chromosome_filter):
        chr_name, start, end, strand, gene_name, transcript_id, gene_id, terminal_status = cds

        cds_processed += 1
        if cds_processed % 50000 == 0:
            print(f"  Matching progress: {cds_processed} CDS features...", end='\r')

        # Create unique CDS identifier based on COORDINATES only (not transcript)
        cds_id = f"{chr_name}:{start}-{end}:{strand}"

        # Get relevant introns (same chromosome and strand)
        relevant_introns = introns_by_chr_strand.get((chr_name, strand), [])

        # Check each intron for boundary matches
        for intron_coord, i_start, i_end, i_counts in relevant_introns:
            matched = False
            boundary_type = None

            if strand == '+':
                # 5' boundary (donor site): CDS ends right before intron starts
                if end == i_start:
                    matched = True
                    boundary_type = '5prime'

                # 3' boundary (acceptor site): CDS starts right after intron ends
                if start == i_end + 1:
                    matched = True
                    boundary_type = '3prime'

            elif strand == '-':
                # 5' boundary (donor site, genomically downstream): CDS starts after intron
                if start == i_end + 1:
                    matched = True
                    boundary_type = '5prime'

                # 3' boundary (acceptor site, genomically upstream): CDS ends before intron
                if end == i_start:
                    matched = True
                    boundary_type = '3prime'

            if matched:
                matches_dict = five_prime_matches if boundary_type == '5prime' else three_prime_matches

                # Store CDS info (will be same for all transcripts with same coordinates)
                if matches_dict[cds_id]['cds_info'] is None:
                    matches_dict[cds_id]['cds_info'] = {
                        'chr': chr_name,
                        'start': start,
                        'end': end,
                        'strand': strand
                    }

                # Accumulate transcripts, genes, gene_ids, introns, and terminal status
                matches_dict[cds_id]['transcripts'].add(transcript_id or 'Unknown')
                matches_dict[cds_id]['genes'].add(gene_name or 'Unknown')
                matches_dict[cds_id]['gene_ids'].add(gene_id or 'Unknown')
                matches_dict[cds_id]['introns'].add((intron_coord, i_counts))
                if terminal_status:
                    matches_dict[cds_id]['terminal_status'].update(terminal_status)

    # Convert sets to sorted lists for consistent output
    for cds_id in five_prime_matches:
        five_prime_matches[cds_id]['transcripts'] = sorted(five_prime_matches[cds_id]['transcripts'])
        five_prime_matches[cds_id]['genes'] = sorted(five_prime_matches[cds_id]['genes'])
        five_prime_matches[cds_id]['gene_ids'] = sorted(five_prime_matches[cds_id]['gene_ids'])
        five_prime_matches[cds_id]['introns'] = sorted(five_prime_matches[cds_id]['introns'])

    for cds_id in three_prime_matches:
        three_prime_matches[cds_id]['transcripts'] = sorted(three_prime_matches[cds_id]['transcripts'])
        three_prime_matches[cds_id]['genes'] = sorted(three_prime_matches[cds_id]['genes'])
        three_prime_matches[cds_id]['gene_ids'] = sorted(three_prime_matches[cds_id]['gene_ids'])
        three_prime_matches[cds_id]['introns'] = sorted(three_prime_matches[cds_id]['introns'])

    print(f"  Total CDS features processed: {cds_processed}")
    print(f"\n  Matching summary (deduplicated by CDS coordinates):")
    print(f"    Unique CDS with 5' boundary matches: {len(five_prime_matches)}")
    print(f"    Unique CDS with 3' boundary matches: {len(three_prime_matches)}")

    return dict(five_prime_matches), dict(three_prime_matches)


def find_constitutive_exons(five_prime_dict, three_prime_dict):
    """
    Find CDS that appear in both 5' and 3' dictionaries.

    These CDS have single introns flanking both boundaries, indicating they are
    "constitutive cassette exons" - exons that are flanked by single introns on both sides.

    Args:
        five_prime_dict: Dictionary of 5' boundary matches
        three_prime_dict: Dictionary of 3' boundary matches

    Returns:
        Dictionary of constitutive exons with both intron coordinates
    """
    constitutive = {}

    for cds_id in five_prime_dict:
        if cds_id in three_prime_dict:
            data_5p = five_prime_dict[cds_id]
            data_3p = three_prime_dict[cds_id]

            constitutive[cds_id] = {
                'cds_info': data_5p['cds_info'],
                'genes': sorted(set(data_5p['genes']) | set(data_3p['genes'])),
                'gene_ids': sorted(set(data_5p['gene_ids']) | set(data_3p['gene_ids'])),
                'transcripts': sorted(set(data_5p['transcripts']) | set(data_3p['transcripts'])),
                'five_prime_introns': data_5p['introns'],
                'three_prime_introns': data_3p['introns']
            }

    print(f"    Constitutive exons (CDS in both dictionaries): {len(constitutive)}")

    return constitutive


def classify_all_exons(five_prime_dict, three_prime_dict):
    """
    Classify ALL exons by their position based on junction presence.

    Position classifications:
    - 'both_junctions': Has both 5' and 3' single-intron junctions
    - 'only_5prime': Has only 5' single-intron junction (may/may not be last exon)
    - 'only_3prime': Has only 3' single-intron junction (may/may not be first exon)

    Note: "only" means we only found one junction type in our filtered data.
    The missing junction may be filtered out (multi-junction cluster, low reads)
    or may not exist (true terminal exon).

    Args:
        five_prime_dict: Dictionary of 5' boundary matches
        three_prime_dict: Dictionary of 3' boundary matches

    Returns:
        Dictionary with all CDS classified by position
    """
    all_exons = {}

    # Get all unique CDS IDs
    all_cds_ids = set(five_prime_dict.keys()) | set(three_prime_dict.keys())

    for cds_id in all_cds_ids:
        in_5p = cds_id in five_prime_dict
        in_3p = cds_id in three_prime_dict

        if in_5p and in_3p:
            # Internal exon: has junctions on both sides
            data_5p = five_prime_dict[cds_id]
            data_3p = three_prime_dict[cds_id]

            all_exons[cds_id] = {
                'cds_info': data_5p['cds_info'],
                'genes': sorted(set(data_5p['genes']) | set(data_3p['genes'])),
                'gene_ids': sorted(set(data_5p['gene_ids']) | set(data_3p['gene_ids'])),
                'transcripts': sorted(set(data_5p['transcripts']) | set(data_3p['transcripts'])),
                'five_prime_introns': data_5p['introns'],
                'three_prime_introns': data_3p['introns'],
                'position': 'internal'
            }

        elif in_5p:
            # Only 5' junction (donor) - include ONLY if it's actually a FIRST exon
            data_5p = five_prime_dict[cds_id]
            if 'first' in data_5p.get('terminal_status', set()):
                all_exons[cds_id] = {
                    'cds_info': data_5p['cds_info'],
                    'genes': data_5p['genes'],
                    'gene_ids': data_5p['gene_ids'],
                    'transcripts': data_5p['transcripts'],
                    'five_prime_introns': data_5p['introns'],
                    'three_prime_introns': [],
                    'position': 'only_5prime'
                }

        elif in_3p:
            # Only 3' junction (acceptor) - include ONLY if it's actually a LAST exon
            data_3p = three_prime_dict[cds_id]
            if 'last' in data_3p.get('terminal_status', set()):
                all_exons[cds_id] = {
                    'cds_info': data_3p['cds_info'],
                    'genes': data_3p['genes'],
                    'gene_ids': data_3p['gene_ids'],
                    'transcripts': data_3p['transcripts'],
                    'five_prime_introns': [],
                    'three_prime_introns': data_3p['introns'],
                    'position': 'only_3prime'
                }

    # Count by position
    position_counts = {'internal': 0, 'only_3prime': 0, 'only_5prime': 0}
    for data in all_exons.values():
        position_counts[data['position']] += 1

    print(f"\n  Exon position classification:")
    print(f"    Internal (both junctions): {position_counts['internal']}")
    print(f"    Only 3' junction: {position_counts['only_3prime']}")
    print(f"    Only 5' junction: {position_counts['only_5prime']}")
    print(f"    Total: {len(all_exons)}")

    return all_exons


def write_bed12(output_file, cds_dict, description, boundary_type=None, exon_position=None):
    """
    Write CDS entries to BED12 format file.

    BED12 format uses 0-based, half-open coordinates [start, end)

    Args:
        output_file: Output file path
        cds_dict: Dictionary of CDS matches (deduplicated by coordinates)
        description: Description for output message
        boundary_type: '5prime', '3prime', or 'constitutive'
        exon_position: Optional position label
    """
    print(f"\n  Writing {description}")
    print(f"    Output: {output_file}")

    entries_written = 0
    strand_counts = {'+': 0, '-': 0}

    with open(output_file, 'w') as f:
        # Write BED header comments
        f.write(f"# {description}\n")
        f.write("# Format: BED12 + exon_position (column 13)\n")
        f.write("# Colors: Blue=internal, Red=only_5prime, Green=only_3prime\n")
        f.write("# Columns: chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts, exon_position\n")

        # Sort entries by chromosome and start position
        sorted_entries = []
        for cds_id, data in cds_dict.items():
            cds_info = data['cds_info']
            sorted_entries.append((cds_info['chr'], cds_info['start'], cds_id, data))

        sorted_entries.sort(key=lambda x: (x[0], x[1]))

        for _, _, cds_id, data in sorted_entries:
            cds_info = data['cds_info']
            chr_name = cds_info['chr']
            start = cds_info['start']
            end = cds_info['end']
            strand = cds_info['strand']

            genes = ','.join(data['genes'])
            gene_ids = ','.join(data.get('gene_ids', ['Unknown']))

            if exon_position is None:
                position = data.get('position', 'unknown')
            else:
                position = exon_position

            if boundary_type == 'constitutive':
                intron_5p_list = [f"{intron}(n={count})" for intron, count in data['five_prime_introns']]
                intron_3p_list = [f"{intron}(n={count})" for intron, count in data['three_prime_introns']]
                intron_info = f"5p:{','.join(intron_5p_list)}|3p:{','.join(intron_3p_list)}"
            else:
                intron_list = [f"{intron}(n={count})" for intron, count in data['introns']]
                if boundary_type == '5prime':
                    intron_info = f"5p:{','.join(intron_list)}"
                elif boundary_type == '3prime':
                    intron_info = f"3p:{','.join(intron_list)}"
                else:
                    intron_info = ','.join(intron_list)

            bed_start = start - 1
            bed_end = end

            if boundary_type == 'constitutive':
                all_counts = [count for _, count in data['five_prime_introns']] + \
                           [count for _, count in data['three_prime_introns']]
                score = min(all_counts) if all_counts else 0
            else:
                all_counts = [count for _, count in data['introns']]
                score = min(all_counts) if all_counts else 0
            score = min(score, 1000)

            if position == 'internal':
                itemRgb = "0,0,255"  # Blue
            elif position == 'only_5prime':
                itemRgb = "255,0,0"  # Red
            elif position == 'only_3prime':
                itemRgb = "0,255,0"  # Green
            else:
                itemRgb = "0,0,0"

            name = f"{genes}|{gene_ids}|{intron_info}"
            bed_line = f"{chr_name}\t{bed_start}\t{bed_end}\t{name}\t{score}\t{strand}\t{bed_start}\t{bed_end}\t{itemRgb}\t1\t{bed_end - bed_start}\t0\t{position}\n"
            f.write(bed_line)
            entries_written += 1
            strand_counts[strand] = strand_counts.get(strand, 0) + 1

    print(f"    Entries written: {entries_written}")
    print(f"      + strand: {strand_counts.get('+', 0)}")
    print(f"      - strand: {strand_counts.get('-', 0)}")


def write_classified_exons_bed(output_file, all_exons_dict, description):
    """
    Write all classified exons to BED12 format with position labels.

    Args:
        output_file: Output file path
        all_exons_dict: Dictionary of all exons with position classification
        description: Description for output message
    """
    print(f"\n  Writing {description}")
    print(f"    Output: {output_file}")

    entries_written = 0
    strand_counts = {'+': 0, '-': 0}
    position_counts = {'internal': 0, 'only_3prime': 0, 'only_5prime': 0}

    with open(output_file, 'w') as f:
        # Write BED header comments
        f.write(f"# {description}\n")
        f.write("# Format: BED12 + exon_position (column 13)\n")
        f.write("# Coordinates: 0-based, half-open (converted from GTF 1-based inclusive)\n")
        f.write("# Name field format: Gene1,Gene2|Transcript1,Transcript2|5p:intron1(n=count1)|3p:intron2\n")
        f.write("# Column 13 (exon_position): only_3prime, only_5prime, both_junctions\n")
        f.write("#   only_3prime: Has single-intron 3' junction (may/may not be first exon)\n")
        f.write("#   only_5prime: Has single-intron 5' junction (may/may not be last exon)\n")
        f.write("#   both_junctions: Has single-intron junctions on both sides (cassette exon)\n")
        f.write("# Columns: chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts, exon_position\n")

        # Sort entries by chromosome and start position
        sorted_entries = []
        for cds_id, data in all_exons_dict.items():
            cds_info = data['cds_info']
            sorted_entries.append((cds_info['chr'], cds_info['start'], cds_id, data))

        sorted_entries.sort(key=lambda x: (x[0], x[1]))

        for _, _, cds_id, data in sorted_entries:
            cds_info = data['cds_info']
            chr_name = cds_info['chr']
            start = cds_info['start']
            end = cds_info['end']
            strand = cds_info['strand']
            position = data['position']

            # Format gene names (comma-separated if multiple)
            genes = ','.join(data['genes'])
            gene_ids = ','.join(data.get('gene_ids', ['Unknown']))

            # Format intron information based on what's available
            intron_parts = []

            if data['five_prime_introns']:
                intron_5p_list = [f"{intron}(n={count})" for intron, count in data['five_prime_introns']]
                intron_parts.append(f"5p:{','.join(intron_5p_list)}")

            if data['three_prime_introns']:
                intron_3p_list = [f"{intron}(n={count})" for intron, count in data['three_prime_introns']]
                intron_parts.append(f"3p:{','.join(intron_3p_list)}")

            intron_info = '|'.join(intron_parts)

            # Calculate numeric score for IGV compatibility
            all_counts = []
            if data['five_prime_introns']:
                all_counts.extend([count for _, count in data['five_prime_introns']])
            if data['three_prime_introns']:
                all_counts.extend([count for _, count in data['three_prime_introns']])
            score = min(all_counts) if all_counts else 0
            score = min(score, 1000)  # BED score must be 0-1000

            # Convert GTF (1-based inclusive) to BED (0-based half-open)
            bed_start = start - 1
            bed_end = end

            # Color by position
            if position == 'internal':
                itemRgb = "0,0,255"  # Blue
            elif position == 'only_5prime':
                itemRgb = "255,0,0"  # Red
            elif position == 'only_3prime':
                itemRgb = "0,255,0"  # Green
            else:
                itemRgb = "0,0,0"

            name = f"{genes}|{gene_ids}|{intron_info}"
            bed_line = f"{chr_name}\t{bed_start}\t{bed_end}\t{name}\t{score}\t{strand}\t{bed_start}\t{bed_end}\t{itemRgb}\t1\t{bed_end - bed_start}\t0\t{position}\n"
            f.write(bed_line)
            entries_written += 1
            strand_counts[strand] = strand_counts.get(strand, 0) + 1
            position_counts[position] = position_counts.get(position, 0) + 1

    print(f"    Entries written: {entries_written}")
    print(f"      + strand: {strand_counts.get('+', 0)}")
    print(f"      - strand: {strand_counts.get('-', 0)}")
    print(f"      Internal: {position_counts.get('internal', 0)}")
    print(f"      Only_3p: {position_counts.get('only_3prime', 0)}")
    print(f"      Only_5p: {position_counts.get('only_5prime', 0)}")


def main():
    parser = argparse.ArgumentParser(
        description='Identify CDS exons adjacent to single-intron Leafcutter clusters',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
COORDINATE HANDLING:
  This script correctly handles coordinate conversion between:
  - Leafcutter (BED-style): 0-based, half-open [start, end)
  - GTF annotation: 1-based, inclusive [start, end]

  For a CDS ending at position 100 (1-based) adjacent to an intron
  starting at [100, ...) (0-based), the match is: CDS_end == intron_start

EXAMPLES:
  # Basic usage with default threshold
  %(prog)s \\
    --leafcutter /path/to/leafcutter2_clusters \\
    --gtf /path/to/annotation.gtf

  # Test on chr1 only with custom threshold
  %(prog)s \\
    --leafcutter /path/to/leafcutter2_clusters \\
    --gtf /path/to/annotation.gtf \\
    --chromosome chr1 \\
    --min_counts 50

  # Full genome analysis with custom output directory
  %(prog)s \\
    --leafcutter /path/to/leafcutter2_clusters \\
    --gtf /path/to/annotation.gtf \\
    --min_counts 100 \\
    --output_dir results/single_intron_analysis/
        """
    )

    parser.add_argument('--leafcutter', '-l',
                        required=True,
                        help='Leafcutter2 clustering output file')
    parser.add_argument('--gtf', '-g',
                        required=True,
                        help='GTF annotation file')
    parser.add_argument('--min_counts', '-m', type=int, default=100,
                        help='Minimum junction read count threshold (default: 100)')
    parser.add_argument('--chromosome', '-c', type=str, default='All',
                        help='Chromosome to analyze (default: All). Use chr1, chr2, etc. for testing')
    parser.add_argument('--output_dir', '-o', type=str, default='.',
                        help='Output directory for BED files (default: current directory)')

    args = parser.parse_args()

    # Validate input files
    if not os.path.exists(args.leafcutter):
        print(f"Error: Leafcutter file not found: {args.leafcutter}", file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(args.gtf):
        print(f"Error: GTF file not found: {args.gtf}", file=sys.stderr)
        sys.exit(1)

    # Create output directory if needed
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
        print(f"\nCreated output directory: {args.output_dir}")

    print("=" * 80)
    print("  Single-Intron CDS Finder - Leafcutter + GTF Analysis")
    print("=" * 80)
    if args.chromosome != 'All':
        print(f"  Analyzing chromosome: {args.chromosome}")
        print("=" * 80)

    # Step 1: Parse Leafcutter clusters
    introns = parse_leafcutter_clusters(
        args.leafcutter,
        args.min_counts,
        args.chromosome
    )

    if not introns:
        print("\nNo single-intron clusters found meeting criteria.", file=sys.stderr)
        print("Try lowering --min_counts threshold or checking --chromosome filter.", file=sys.stderr)
        sys.exit(0)

    # Step 2: Match CDS to introns
    five_prime_dict, three_prime_dict = match_cds_to_introns(
        args.gtf,
        introns,
        args.chromosome
    )

    # Step 3: Find constitutive exons
    constitutive_dict = find_constitutive_exons(five_prime_dict, three_prime_dict)

    # Step 4: Classify all exons by position
    all_exons_dict = classify_all_exons(five_prime_dict, three_prime_dict)

    # Step 5: Write output BED12 files
    print("\nGenerating output files...")

    # Add chromosome suffix to filenames if filtering
    chr_suffix = f"_{args.chromosome}" if args.chromosome != 'All' else ""

    output_5p = os.path.join(args.output_dir, f"cds_5prime_adjacent_introns{chr_suffix}.bed")
    output_3p = os.path.join(args.output_dir, f"cds_3prime_adjacent_introns{chr_suffix}.bed")
    output_const = os.path.join(args.output_dir, f"constitutive_exons{chr_suffix}.bed")
    output_all = os.path.join(args.output_dir, f"all_exons_classified{chr_suffix}.bed")
    output_terminal = os.path.join(args.output_dir, f"terminal_exons{chr_suffix}.bed")

    write_bed12(output_5p, five_prime_dict,
                "CDS with 5' boundary adjacent to single introns", '5prime', exon_position='only_5prime')
    write_bed12(output_3p, three_prime_dict,
                "CDS with 3' boundary adjacent to single introns", '3prime', exon_position='only_3prime')
    write_bed12(output_const, constitutive_dict,
                "Constitutive cassette exons (flanked by single introns)", 'constitutive', exon_position='internal')
    write_classified_exons_bed(output_all, all_exons_dict,
                              "All exons classified by position (first/internal/last)")

    # Terminal exons only (first and last)
    terminal_exons_dict = {k: v for k, v in all_exons_dict.items()
                          if v.get('position') in ['only_5prime', 'only_3prime']}
    write_classified_exons_bed(output_terminal, terminal_exons_dict,
                              "Terminal exons only (first and last exons)")

    print("\n" + "=" * 80)
    print("  ANALYSIS COMPLETE")
    print("=" * 80)
    print(f"\nOutput files:")
    print(f"  {output_5p}  (only_5prime: first exons)")
    print(f"  {output_3p}  (only_3prime: last exons)")
    print(f"  {output_const}  (internal: cassette exons)")
    print(f"  {output_terminal}  (first + last exons only)")
    print(f"  {output_all}  (ALL: internal + first + last)")
    print()


if __name__ == '__main__':
    main()
