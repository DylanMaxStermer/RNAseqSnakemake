#!/usr/bin/env python3
"""
Annotate constitutive exons with CDS reading frame (phase) from a reference GTF.

Matches each exon in exons_std.tsv to a GTF CDS feature using exact coordinate
lookup and records the phase (0, 1, or 2) from the GTF frame column.

Coordinate conversion:
  exons_std.tsv uses BED (0-based half-open): start, end
  GTF CDS uses 1-based inclusive:             gtf_start = bed_start + 1, gtf_end = bed_end

Author: Dylan Stermer
"""

import argparse
import os
import sys
from collections import defaultdict


def load_cds_phase(gtf_file):
    """
    Stream the GTF and build a lookup dict of CDS phase values.

    Returns:
        dict: {(chrom, start_1based, end_1based, strand): [phase_int, ...]}
    """
    cds_phase = defaultdict(list)
    cds_count = 0

    print(f"Loading GTF: {gtf_file}")
    with open(gtf_file) as fh:
        for line_num, line in enumerate(fh, 1):
            if line_num % 1_000_000 == 0:
                print(f"  Processed {line_num:,} lines...", end="\r")
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            if fields[2] != "CDS":
                continue

            chrom  = fields[0]
            start  = int(fields[3])   # 1-based inclusive
            end    = int(fields[4])   # 1-based inclusive
            strand = fields[6]
            frame  = fields[7]        # '0', '1', '2'

            if frame not in ("0", "1", "2"):
                continue

            key = (chrom, start, end, strand)
            cds_phase[key].append(int(frame))
            cds_count += 1

    print(f"\n  CDS features loaded: {cds_count:,}")
    print(f"  Unique CDS coordinate keys: {len(cds_phase):,}")
    return cds_phase


def annotate_exons(exon_file, cds_phase, output_file):
    """
    Read exons_std.tsv, match each row to a CDS phase, and write output.
    """
    # Read input, preserving comment header lines
    comments = []
    header = None
    rows = []

    with open(exon_file) as fh:
        for line in fh:
            if line.startswith("#"):
                comments.append(line.rstrip("\n"))
                continue
            if header is None:
                header = line.rstrip("\n").split("\t")
                continue
            rows.append(line.rstrip("\n").split("\t"))

    if header is None:
        print("Error: no header row found in input file.", file=sys.stderr)
        sys.exit(1)

    print(f"Read {len(rows)} exon rows")

    # Column indices
    try:
        i_chrom  = header.index("chrom")
        i_start  = header.index("start")
        i_end    = header.index("end")
        i_strand = header.index("strand")
    except ValueError as e:
        print(f"Error: required column missing: {e}", file=sys.stderr)
        sys.exit(1)

    matched = 0
    unmatched = 0
    phase_values = []

    for row in rows:
        chrom  = row[i_chrom]
        start  = int(row[i_start])   # BED 0-based
        end    = int(row[i_end])     # BED half-open
        strand = row[i_strand]

        # Convert to GTF 1-based inclusive
        gtf_start = start + 1
        gtf_end   = end

        key = (chrom, gtf_start, gtf_end, strand)
        phases = cds_phase.get(key)

        if phases:
            phase = str(min(phases))
            matched += 1
        else:
            phase = "NA"
            unmatched += 1

        phase_values.append(phase)

    print(f"  Matched:   {matched:,}")
    print(f"  Unmatched: {unmatched:,}")

    # Write output
    os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)

    out_header = header + ["phase"]

    with open(output_file, "w") as fh:
        for c in comments:
            fh.write(c + "\n")
        # Update total columns comment if present
        fh.write("\t".join(out_header) + "\n")
        for row, phase in zip(rows, phase_values):
            fh.write("\t".join(row + [phase]) + "\n")

    print(f"Wrote {len(rows)} rows to {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Annotate constitutive exons with CDS phase from GTF"
    )
    parser.add_argument(
        "-e", "--exons", required=True,
        help="Input exons_std.tsv file",
    )
    parser.add_argument(
        "-g", "--gtf", required=True,
        help="Reference GTF file",
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output file (exons_std_phased.tsv)",
    )
    args = parser.parse_args()

    for path, label in [(args.exons, "exons"), (args.gtf, "GTF")]:
        if not os.path.exists(path):
            print(f"Error: {label} file not found: {path}", file=sys.stderr)
            sys.exit(1)

    cds_phase = load_cds_phase(args.gtf)
    annotate_exons(args.exons, cds_phase, args.output)


if __name__ == "__main__":
    main()
