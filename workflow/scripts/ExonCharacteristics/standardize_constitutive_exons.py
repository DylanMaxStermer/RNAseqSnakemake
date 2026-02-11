#!/usr/bin/env python3
"""
Standardize constitutive_exons.bed into the common exons_std.tsv format.

Input BED12 + exon_position columns:
  chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd,
  itemRgb, blockCount, blockSizes, blockStarts, exon_position

Name field format:
  NOC2L|ENSG00000188976.11|5p:chr1:951238-951999(n=1008)|3p:chr1:952139-952411(n=967)

Output standard columns (first 10) + junction columns:
  chrom, start, end, name, strand, length, gene_name, gene_type, tx_type,
  exon_position, 5primeJunc, 3primeJunc

Author: Dylan Stermer
"""

import argparse
import os
import re
import sys


INPUT_COLS = [
    "chrom", "chromStart", "chromEnd", "name", "score", "strand",
    "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes",
    "blockStarts", "exon_position",
]


def parse_name_field(name):
    """
    Parse the composite name field into its components.

    Format: GENE|GENE_ID|5p:coord(n=N)|3p:coord(n=N)

    Returns:
        gene_name (str), five_prime_junc (str), three_prime_junc (str)
        Missing fields are returned as 'NA'.
    """
    parts = name.split("|")

    gene_name = parts[0] if parts else "NA"

    five_prime_junc = "NA"
    three_prime_junc = "NA"

    for part in parts:
        if part.startswith("5p:"):
            # Strip prefix and trailing (n=...)
            coord = re.sub(r"\(n=\d+\)$", "", part[3:])
            five_prime_junc = coord
        elif part.startswith("3p:"):
            coord = re.sub(r"\(n=\d+\)$", "", part[3:])
            three_prime_junc = coord

    return gene_name, five_prime_junc, three_prime_junc


def standardize(input_file, output_file):
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    rows = []
    with open(input_file) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < len(INPUT_COLS):
                # pad if exon_position column is missing
                fields += ["NA"] * (len(INPUT_COLS) - len(fields))
            row = dict(zip(INPUT_COLS, fields))
            rows.append(row)

    print(f"Read {len(rows)} data rows from {input_file}")

    out_rows = []
    for row in rows:
        chrom = row["chrom"]
        start = int(row["chromStart"])
        end = int(row["chromEnd"])
        length = end - start
        strand = row["strand"]
        exon_position = row["exon_position"]

        gene_name, five_junc, three_junc = parse_name_field(row["name"])

        out_rows.append({
            "chrom": chrom,
            "start": start,
            "end": end,
            "name": gene_name,
            "strand": strand,
            "length": length,
            "gene_name": gene_name,
            "gene_type": "protein_coding",
            "tx_type": "protein_coding",
            "exon_position": exon_position,
            "5primeJunc": five_junc,
            "3primeJunc": three_junc,
        })

    out_cols = [
        "chrom", "start", "end", "name", "strand", "length",
        "gene_name", "gene_type", "tx_type", "exon_position",
        "5primeJunc", "3primeJunc",
    ]

    with open(output_file, "w") as fh:
        fh.write("# Standardized exon annotation file\n")
        fh.write("# Standard columns (first 10): chrom, start, end, name, strand, length, "
                 "gene_name, gene_type, tx_type, exon_position\n")
        fh.write("# Coordinates: 0-based half-open (BED format)\n")
        fh.write("# Original format: constitutive_bed\n")
        fh.write(f"# Total columns: {len(out_cols)}\n")
        fh.write("\t".join(out_cols) + "\n")
        for row in out_rows:
            fh.write("\t".join(str(row[c]) for c in out_cols) + "\n")

    print(f"Wrote {len(out_rows)} rows to {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Standardize constitutive_exons.bed to exons_std.tsv format"
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Input constitutive_exons.bed file",
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output standardized TSV file",
    )
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"Error: input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    standardize(args.input, args.output)


if __name__ == "__main__":
    main()
