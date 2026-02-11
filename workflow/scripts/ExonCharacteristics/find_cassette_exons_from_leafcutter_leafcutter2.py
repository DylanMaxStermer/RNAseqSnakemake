#!/usr/bin/env python3
"""
Cassette Exon Identification from LeafCutter2 Data

Identifies potential cassette exons from LeafCutter2 junction clusters.
Does not require a skipping junction to be present.

Usage:
    python find_cassette_exons_leafcutter2.py \
        --leafcutter-counts counts.tsv \
        --leafcutter-annotation annotation.tsv \
        --output-dir output/
"""

import argparse
import os
import re
import sys
from pathlib import Path

import pandas as pd
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(
        description="Identify cassette exons from LeafCutter2 junction data",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Input files
    parser.add_argument(
        "--leafcutter-counts", "-c",
        required=True,
        help="Path to LeafCutter2 junction counts file (space-delimited, junction IDs as row names)"
    )
    parser.add_argument(
        "--leafcutter-annotation", "-a",
        required=True,
        help="Path to junction annotation file (TSV with Gene_name, Intron_coord, Strand, etc.)"
    )

    # Output
    parser.add_argument(
        "--output-dir", "-o",
        default="output",
        help="Output directory for results"
    )

    # Filtering parameters
    parser.add_argument(
        "--min-cassette-length",
        type=int,
        default=30,
        help="Minimum cassette exon length in bp"
    )
    parser.add_argument(
        "--max-cassette-length",
        type=int,
        default=500,
        help="Maximum cassette exon length in bp"
    )
    parser.add_argument(
        "--min-junction-reads",
        type=int,
        default=5,
        help="Minimum total reads across all samples for a junction"
    )

    # Annotation filtering flags (each filters independently)
    parser.add_argument(
        "--keep-annot",
        type=lambda x: x.lower() in ('true', '1', 'yes'),
        default=True,
        help="Require Annot=True (keep only annotated junctions)"
    )
    parser.add_argument(
        "--keep-coding",
        type=lambda x: x.lower() in ('true', '1', 'yes'),
        default=True,
        help="Require Coding=True (keep only coding junctions)"
    )
    parser.add_argument(
        "--exclude-utr",
        type=lambda x: x.lower() in ('true', '1', 'yes'),
        default=True,
        help="Exclude UTR junctions (keep only UTR=False)"
    )
    parser.add_argument(
        "--keep-gencodepc",
        type=lambda x: x.lower() in ('true', '1', 'yes'),
        default=True,
        help="Require GencodePC=True (keep only GENCODE protein-coding)"
    )

    return parser.parse_args()


def load_and_parse_counts(counts_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load LeafCutter2 counts and parse junction IDs to extract cluster info.

    Junction ID format: chr:start:stop:clu_N_strand:TYPE
    Example: chr1:961729:961825:clu_3_+:PR

    Note: LeafCutter2 uses actual stop position (no +1 offset like old LeafCutter).
    """

    print("Loading and parsing junction counts...")

    # Load counts - space delimited with row names
    # LeafCutter2 has a header row with "chrom" as first column
    counts = pd.read_csv(counts_path, sep=r'\s+', index_col=0)
    print(f"  Loaded counts: {counts.shape[0]} junctions x {counts.shape[1]} samples")

    # Parse junction IDs
    junction_ids = counts.index.tolist()

    parsed_rows = []
    for jid in junction_ids:
        # Format: chr:start:stop:clu_N_strand:TYPE
        parts = jid.split(":")
        if len(parts) != 5:
            print(f"  Warning: skipping malformed junction ID: {jid}")
            continue

        chrom = parts[0]
        start = int(parts[1])
        stop = int(parts[2])  # LeafCutter2 uses actual stop, no +1
        clu_strand = parts[3]
        junction_type = parts[4]

        # Parse clu_N_strand (e.g., "clu_1_-")
        match = re.match(r"(clu_\d+)_([+-])$", clu_strand)
        if not match:
            print(f"  Warning: skipping malformed cluster/strand: {clu_strand}")
            continue

        clu = match.group(1)
        strand = match.group(2)

        # Create Intron_coord for joining with annotation (no coordinate adjustment needed)
        intron_coord = f"{chrom}:{start}-{stop}"

        parsed_rows.append({
            "original_rowname": jid,
            "clu": clu,
            "chrom": chrom,
            "start": start,
            "stop": stop,
            "strand": strand,
            "junction_type": junction_type,
            "Intron_coord": intron_coord
        })

    parsed_df = pd.DataFrame(parsed_rows)
    print(f"  Parsed {len(parsed_df)} junction coordinates")

    return counts, parsed_df


def load_annotation(annotation_path: str) -> pd.DataFrame:
    """Load junction classification/annotation file."""

    print("Loading junction annotation...")
    annotation = pd.read_csv(annotation_path, sep="\t")
    print(f"  Loaded annotation: {annotation.shape[0]} junctions")
    print(f"  Columns: {list(annotation.columns)}")

    return annotation


def merge_counts_with_annotation(
    parsed_counts: pd.DataFrame,
    annotation: pd.DataFrame
) -> pd.DataFrame:
    """Join parsed count data with annotation on Intron_coord."""

    print("Merging counts with annotation...")

    merged = parsed_counts.merge(
        annotation,
        on="Intron_coord",
        how="left"
    )

    # Check for unmatched junctions
    unmatched = merged["Gene_name"].isna().sum()
    if unmatched > 0:
        print(f"  Warning: {unmatched} junctions have no annotation match")

    matched = (~merged["Gene_name"].isna()).sum()
    print(f"  Matched {matched} junctions with annotation")

    return merged


def filter_by_annotation(
    df: pd.DataFrame,
    keep_annot: bool,
    keep_coding: bool,
    exclude_utr: bool,
    keep_gencodepc: bool
) -> pd.DataFrame:
    """
    Filter junctions based on annotation flags.

    Logic mirrors the R code:
        filter(Coding == TRUE, UTR == FALSE, Annot == TRUE)

    Each flag enables that filter. If flag is True, apply the filter.
    """

    initial_count = len(df)
    mask = pd.Series([True] * len(df), index=df.index)

    filters_applied = []

    if keep_annot and "Annot" in df.columns:
        mask &= df["Annot"].fillna(False).astype(bool)
        filters_applied.append("Annot=True")

    if keep_coding and "Coding" in df.columns:
        mask &= df["Coding"].fillna(False).astype(bool)
        filters_applied.append("Coding=True")

    if exclude_utr and "UTR" in df.columns:
        mask &= ~df["UTR"].fillna(True).astype(bool)
        filters_applied.append("UTR=False (excluded)")

    if keep_gencodepc and "GencodePC" in df.columns:
        mask &= df["GencodePC"].fillna(False).astype(bool)
        filters_applied.append("GencodePC=True")

    filtered = df[mask].copy()

    print(f"Annotation filtering: {initial_count} -> {len(filtered)} junctions")
    if filters_applied:
        print(f"  Filters applied: {', '.join(filters_applied)}")
    else:
        print("  No annotation filters applied")

    return filtered


def find_cassette_exons(
    counts: pd.DataFrame,
    cluster_coords: pd.DataFrame,
    min_cassette_length: int,
    max_cassette_length: int,
    min_junction_reads: int
) -> pd.DataFrame:
    """Identify cassette exons from junction data."""

    print("\nStep 1: Calculating junction read counts...")

    # Calculate total reads per junction
    junction_reads = counts.sum(axis=1).to_dict()

    # Also create a TYPE-agnostic lookup for skip junctions
    # Key format: chr:start:stop:clu_N_strand (without :TYPE)
    junction_reads_no_type = {}
    for jid, reads in junction_reads.items():
        # Remove the :TYPE suffix
        parts = jid.rsplit(":", 1)
        if len(parts) == 2:
            base_id = parts[0]
            # If we see this base_id multiple times (different TYPEs), sum the reads
            junction_reads_no_type[base_id] = junction_reads_no_type.get(base_id, 0) + reads

    cluster_coords["total_reads"] = cluster_coords["original_rowname"].map(junction_reads)

    # Filter by minimum reads
    junctions = cluster_coords[
        cluster_coords["total_reads"] >= min_junction_reads
    ].copy()

    # Also require annotation match (Gene_name not null)
    junctions = junctions[junctions["Gene_name"].notna()].copy()

    print(f"  Found {len(junctions)} junctions with >= {min_junction_reads} reads and annotation")
    plus_count = (junctions["strand"] == "+").sum()
    minus_count = (junctions["strand"] == "-").sum()
    print(f"  Plus strand: {plus_count} | Minus strand: {minus_count}")

    if len(junctions) == 0:
        return pd.DataFrame()

    print("\nStep 2: Finding junction pairs in clusters...")

    # Prepare for self-join
    j1 = junctions[["clu", "chrom", "strand", "Gene_name", "start", "stop",
                    "original_rowname", "total_reads"]].copy()
    j1.columns = ["clu", "chrom", "Strand", "Gene_name", "j1_start", "j1_stop",
                  "junction_1", "reads_1"]

    j2 = junctions[["clu", "chrom", "strand", "start", "stop",
                    "original_rowname", "total_reads"]].copy()
    j2.columns = ["clu", "chrom", "Strand", "j2_start", "j2_stop",
                  "junction_2", "reads_2"]

    # Self-join on cluster, chromosome, strand
    cassettes = j1.merge(j2, on=["clu", "chrom", "Strand"], how="inner")
    cassettes = cassettes[cassettes["junction_1"] != cassettes["junction_2"]].copy()

    print(f"  Found {len(cassettes)} potential junction pairs")

    if len(cassettes) == 0:
        return pd.DataFrame()

    print("\nStep 3: Filtering for non-overlapping junction pairs...")

    # Keep only non-overlapping pairs (true cassette exons)
    # Plus strand: upstream junction ends before downstream junction starts
    # Minus strand: downstream junction ends before upstream junction starts
    def is_j1_upstream(row):
        if row["Strand"] == "+":
            return row["j1_stop"] < row["j2_start"]
        elif row["Strand"] == "-":
            return row["j2_stop"] < row["j1_start"]
        return False

    cassettes["is_j1_upstream"] = cassettes.apply(is_j1_upstream, axis=1)
    cassettes = cassettes[cassettes["is_j1_upstream"]].copy()

    # Rename for clarity
    cassettes["upstream_start"] = cassettes["j1_start"]
    cassettes["upstream_stop"] = cassettes["j1_stop"]
    cassettes["downstream_start"] = cassettes["j2_start"]
    cassettes["downstream_stop"] = cassettes["j2_stop"]
    cassettes["junction_upstream"] = cassettes["junction_1"]
    cassettes["junction_downstream"] = cassettes["junction_2"]
    cassettes["reads_upstream"] = cassettes["reads_1"]
    cassettes["reads_downstream"] = cassettes["reads_2"]

    print(f"  Kept {len(cassettes)} non-overlapping pairs")

    if len(cassettes) == 0:
        return pd.DataFrame()

    print("\nStep 4: Calculating cassette exon coordinates...")

    # Calculate cassette coordinates (the exon between the two introns)
    def calc_cassette_coords(row):
        if row["Strand"] == "+":
            return row["upstream_stop"], row["downstream_start"]
        else:
            return row["downstream_stop"], row["upstream_start"]

    coords = cassettes.apply(calc_cassette_coords, axis=1, result_type="expand")
    cassettes["cassette_start"] = coords[0]
    cassettes["cassette_stop"] = coords[1]
    cassettes["cassette_length"] = cassettes["cassette_stop"] - cassettes["cassette_start"]
    cassettes["cassette_coord"] = (
        cassettes["chrom"] + ":" +
        cassettes["cassette_start"].astype(str) + "-" +
        cassettes["cassette_stop"].astype(str)
    )

    # Filter by length
    cassettes = cassettes[
        (cassettes["cassette_length"] >= min_cassette_length) &
        (cassettes["cassette_length"] <= max_cassette_length)
    ].copy()

    print(f"  Found {len(cassettes)} cassettes with valid lengths ({min_cassette_length}-{max_cassette_length} bp)")
    plus_count = (cassettes["Strand"] == "+").sum()
    minus_count = (cassettes["Strand"] == "-").sum()
    print(f"  Plus strand: {plus_count} | Minus strand: {minus_count}")

    if len(cassettes) == 0:
        return pd.DataFrame()

    print("\nStep 5: Checking for skipping junctions...")

    # Calculate skip junction coordinates
    def calc_skip_coords(row):
        if row["Strand"] == "+":
            return row["upstream_start"], row["downstream_stop"]
        else:
            return row["downstream_start"], row["upstream_stop"]

    skip_coords = cassettes.apply(calc_skip_coords, axis=1, result_type="expand")
    cassettes["skip_start"] = skip_coords[0]
    cassettes["skip_stop"] = skip_coords[1]

    # Extract cluster number for expected junction ID
    cassettes["clu_num"] = cassettes["clu"].str.extract(r"clu_(\d+)")[0]

    # Construct expected skipping junction ID (without TYPE suffix)
    # Format: chr:start:stop:clu_N_strand
    cassettes["skip_junction_base"] = (
        cassettes["chrom"] + ":" +
        cassettes["skip_start"].astype(str) + ":" +
        cassettes["skip_stop"].astype(str) + ":" +
        "clu_" + cassettes["clu_num"] + "_" + cassettes["Strand"]
    )

    # Look up reads for skipping junction (TYPE-agnostic)
    cassettes["reads_skip"] = cassettes["skip_junction_base"].map(
        lambda x: junction_reads_no_type.get(x, 0)
    )
    cassettes["skip_detected"] = cassettes["reads_skip"] >= min_junction_reads

    skip_count = cassettes["skip_detected"].sum()
    skip_pct = 100 * cassettes["skip_detected"].mean()
    print(f"  Cassettes WITH skipping junction: {skip_count} ({skip_pct:.1f}%)")

    plus_skip = ((cassettes["skip_detected"]) & (cassettes["Strand"] == "+")).sum()
    minus_skip = ((cassettes["skip_detected"]) & (cassettes["Strand"] == "-")).sum()
    print(f"    Plus: {plus_skip} | Minus: {minus_skip}")

    no_skip_count = (~cassettes["skip_detected"]).sum()
    no_skip_pct = 100 * (~cassettes["skip_detected"]).mean()
    print(f"  Cassettes WITHOUT skipping junction: {no_skip_count} ({no_skip_pct:.1f}%)")

    plus_no_skip = ((~cassettes["skip_detected"]) & (cassettes["Strand"] == "+")).sum()
    minus_no_skip = ((~cassettes["skip_detected"]) & (cassettes["Strand"] == "-")).sum()
    print(f"    Plus: {plus_no_skip} | Minus: {minus_no_skip}")

    # Create cassette ID
    cassettes["cassette_id"] = (
        cassettes["Gene_name"] + "_" +
        cassettes["chrom"] + ":" +
        cassettes["cassette_start"].astype(str) + "-" +
        cassettes["cassette_stop"].astype(str) + "_" +
        cassettes["Strand"]
    )

    # Sort by position
    cassettes = cassettes.sort_values(["chrom", "cassette_start"]).reset_index(drop=True)

    return cassettes


def export_results(cassettes: pd.DataFrame, output_dir: str) -> pd.DataFrame:
    """Export results to various file formats."""

    os.makedirs(output_dir, exist_ok=True)

    print("\nExporting results...")

    # Main results table
    cassette_table = cassettes[[
        "cassette_id", "Gene_name", "cassette_coord", "chrom",
        "cassette_start", "cassette_stop", "cassette_length", "Strand",
        "clu", "junction_upstream", "junction_downstream",
        "skip_junction_base", "reads_upstream", "reads_downstream",
        "reads_skip", "skip_detected"
    ]].copy()

    cassette_table.columns = [
        "cassette_id", "gene_name", "cassette_coord", "chrom",
        "cassette_start", "cassette_stop", "cassette_length", "strand",
        "cluster", "upstream_junction", "downstream_junction",
        "skipping_junction", "reads_upstream", "reads_downstream",
        "reads_skip", "skip_detected"
    ]

    output_path = os.path.join(output_dir, "cassetteExons.tsv")
    cassette_table.to_csv(output_path, sep="\t", index=False)
    print(f"  Saved: {output_path}")

    # Filtered version with only cassettes that have skipping junctions
    cassette_table_with_skipping = cassette_table[cassette_table["skip_detected"]].copy()
    output_path = os.path.join(output_dir, "cassetteExons_withSkipping.tsv")
    cassette_table_with_skipping.to_csv(output_path, sep="\t", index=False)
    print(f"  Saved: {output_path} ({len(cassette_table_with_skipping)} cassettes)")

    # BED file of cassette exon coordinates
    bed_cassettes = pd.DataFrame({
        "chrom": cassettes["chrom"],
        "start": cassettes["cassette_start"],
        "stop": cassettes["cassette_stop"],
        "name": cassettes["cassette_id"],
        "score": cassettes["reads_upstream"] + cassettes["reads_downstream"],
        "strand": cassettes["Strand"],
        "thickStart": cassettes["cassette_start"],
        "thickEnd": cassettes["cassette_stop"],
        "itemRgb": cassettes["skip_detected"].apply(
            lambda x: "255,0,0" if x else "150,150,150"
        )
    })

    output_path = os.path.join(output_dir, "cassetteExons.bed")
    bed_cassettes.to_csv(output_path, sep="\t", index=False, header=False)
    print(f"  Saved: {output_path}")

    # BED file with all junction coordinates
    bed_upstream = pd.DataFrame({
        "chrom": cassettes["chrom"],
        "start": cassettes["upstream_start"],
        "stop": cassettes["upstream_stop"],
        "name": "upstream_" + cassettes["cassette_id"],
        "score": cassettes["reads_upstream"],
        "strand": cassettes["Strand"],
        "thickStart": cassettes["upstream_start"],
        "thickEnd": cassettes["upstream_stop"],
        "itemRgb": "0,200,0"
    })

    bed_downstream = pd.DataFrame({
        "chrom": cassettes["chrom"],
        "start": cassettes["downstream_start"],
        "stop": cassettes["downstream_stop"],
        "name": "downstream_" + cassettes["cassette_id"],
        "score": cassettes["reads_downstream"],
        "strand": cassettes["Strand"],
        "thickStart": cassettes["downstream_start"],
        "thickEnd": cassettes["downstream_stop"],
        "itemRgb": "0,0,200"
    })

    cassettes_with_skip = cassettes[cassettes["skip_detected"]]
    bed_skip = pd.DataFrame({
        "chrom": cassettes_with_skip["chrom"],
        "start": cassettes_with_skip["skip_start"],
        "stop": cassettes_with_skip["skip_stop"],
        "name": "skip_" + cassettes_with_skip["cassette_id"],
        "score": cassettes_with_skip["reads_skip"],
        "strand": cassettes_with_skip["Strand"],
        "thickStart": cassettes_with_skip["skip_start"],
        "thickEnd": cassettes_with_skip["skip_stop"],
        "itemRgb": "255,150,0"
    })

    bed_junctions = pd.concat([bed_upstream, bed_downstream, bed_skip], ignore_index=True)
    bed_junctions = bed_junctions.sort_values(["chrom", "start"]).reset_index(drop=True)

    output_path = os.path.join(output_dir, "cassetteExons_junctions.bed")
    bed_junctions.to_csv(output_path, sep="\t", index=False, header=False)
    print(f"  Saved: {output_path}")

    # Per-gene summary
    gene_summary = cassette_table.groupby("gene_name").agg(
        n_cassettes_total=("cassette_id", "count"),
        n_cassettes_with_skipping=("skip_detected", "sum"),
        n_cassettes_without_skipping=("skip_detected", lambda x: (~x).sum()),
        median_cassette_length=("cassette_length", "median")
    ).reset_index().sort_values("n_cassettes_total", ascending=False)

    output_path = os.path.join(output_dir, "cassetteExons_per_gene.tsv")
    gene_summary.to_csv(output_path, sep="\t", index=False)
    print(f"  Saved: {output_path}")

    # Summary statistics
    print("\n=== SUMMARY ===")
    print(f"Total cassette exons identified: {len(cassette_table)}")
    print(f"  Plus strand: {(cassette_table['strand'] == '+').sum()}")
    print(f"  Minus strand: {(cassette_table['strand'] == '-').sum()}")

    skip_count = cassette_table["skip_detected"].sum()
    skip_pct = 100 * cassette_table["skip_detected"].mean()
    print(f"\nWith detectable skipping: {skip_count} ({skip_pct:.1f}%)")

    no_skip_count = (~cassette_table["skip_detected"]).sum()
    no_skip_pct = 100 * (~cassette_table["skip_detected"]).mean()
    print(f"Without detectable skipping: {no_skip_count} ({no_skip_pct:.1f}%)")

    print(f"\nGenes containing cassette exons: {cassette_table['gene_name'].nunique()}")
    print(f"Median cassette length: {cassette_table['cassette_length'].median():.0f} bp")
    print(f"Cassette length range: {cassette_table['cassette_length'].min():.0f} - "
          f"{cassette_table['cassette_length'].max():.0f} bp")

    print("\nTop 10 genes by cassette count:")
    top_genes = gene_summary.head(10)[["gene_name", "n_cassettes_total", "n_cassettes_with_skipping"]]
    print(top_genes.to_string(index=False))

    return cassette_table


def main():
    print("========================================")
    print("   CASSETTE EXON IDENTIFICATION")
    print("   (LeafCutter2 Format)")
    print("========================================\n")

    args = parse_args()

    # Step 1: Load and parse counts (extracts cluster info)
    counts, parsed_counts = load_and_parse_counts(args.leafcutter_counts)

    # Step 2: Load annotation
    annotation = load_annotation(args.leafcutter_annotation)

    # Step 3: Merge counts with annotation
    cluster_coords = merge_counts_with_annotation(parsed_counts, annotation)

    # Step 4: Filter by annotation flags
    cluster_coords_filtered = filter_by_annotation(
        cluster_coords,
        keep_annot=args.keep_annot,
        keep_coding=args.keep_coding,
        exclude_utr=args.exclude_utr,
        keep_gencodepc=args.keep_gencodepc
    )

    # Step 5: Find cassette exons
    cassettes = find_cassette_exons(
        counts=counts,
        cluster_coords=cluster_coords_filtered,
        min_cassette_length=args.min_cassette_length,
        max_cassette_length=args.max_cassette_length,
        min_junction_reads=args.min_junction_reads
    )

    if len(cassettes) == 0:
        print("\nNo cassette exons found with current parameters.")
        sys.exit(0)

    # Step 6: Export results
    results = export_results(cassettes, args.output_dir)

    print("\n========================================")
    print("   COMPLETE")
    print("========================================")
    print(f"\nOutput files saved to: {args.output_dir}/")
    print("  - cassetteExons.tsv (all cassettes)")
    print("  - cassetteExons_withSkipping.tsv (only cassettes with skipping)")
    print("  - cassetteExons.bed (cassette coordinates)")
    print("  - cassetteExons_junctions.bed (all junctions)")
    print("  - cassetteExons_per_gene.tsv (gene summary)")


if __name__ == "__main__":
    main()
