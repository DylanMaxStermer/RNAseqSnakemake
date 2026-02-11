# Splice Site Score Heatmap Generator

## Overview

This tool integrates splice site analysis data with reference splice site scores to create visualizations showing how different scoring metrics (Balance, Splice Score, LUC7) correlate with cryptic splice site usage in both exon lengthening and shortening events.

## What It Does

1. **Extracts 9-mer sequences** from canonical and cryptic splice sites (3 exon bases + 6 intron bases)
2. **Matches** these 9-mers against a reference database of scored splice sites
3. **Creates heatmaps** showing score patterns across ranked events
4. **Optionally clusters** events by canonical site metrics instead of by rank

## Input Files Required

### 1. Reference Score File
**File:** `conor_handedness.tsv`
- **Format:** Tab-separated values
- **Size:** 262,144 scored 9-mer sequences
- **Key Columns:**
  - `SE.5SS.bed.seq`: 9-mer splice site sequence (e.g., "CAGGTGAGG")
  - `Balance`: Handedness metric (-3 to 3)
  - `SE.5SS.bed.score`: Splice site prediction score
  - `LUC7_score`: LUC7 protein binding score

### 2. Analysis Files
**Files:**
- `splice_site_analysis_longer.tsv` (exon lengthening events, ΔSJ > 0)
- `splice_site_analysis_shorter.tsv` (exon shortening events, ΔSJ < 0)

**Format:** Tab-separated values from splice site grid analysis

**Key Columns:**
- `canonical_sequence`: Canonical splice site with format "EXONBASES|INTRONBASES" (e.g., "GTTCCCACAG|GTATTCACAT")
- `cryptic_sequence`: Cryptic/alternative splice site in same format
- `cryptic_length`: Size difference between canonical and cryptic
- `rank`: Event ranking by |ΔSJ| (1 = largest)

## 9-mer Extraction Rule

The script extracts a 9-base sequence from the longer sequences:

```
Input:  "GTTCCCACAG|GTATTCACAT"
         └─exon─┘  └─intron─┘

Extract: Last 3 exon bases + First 6 intron bases
         "CAG" + "GTATTC" = "CAGGTATTC"
```

This 9-mer is then matched against the reference file to retrieve scores.

## Command-Line Flags

### Required Arguments

| Flag | Description | Example |
|------|-------------|---------|
| `--reference` | Path to reference score file (conor_handedness.tsv) | `--reference conor_handedness.tsv` |
| `--longer` | Path to longer events analysis file | `--longer splice_site_analysis_longer.tsv` |
| `--shorter` | Path to shorter events analysis file | `--shorter splice_site_analysis_shorter.tsv` |
| `--output` | Output directory for results | `--output heatmap_results/` |

### Optional Arguments

| Flag | Default | Choices | Description |
|------|---------|---------|-------------|
| `--cluster` | `none` | `none`, `balance`, `score`, `luc7` | Clustering method based on canonical site metrics |
| `--scores` | All three | Any subset | Score columns to visualize (for future use) |

### Clustering Options Explained

- **`none`** (default): Events displayed in original rank order (by |ΔSJ|, largest to smallest)
- **`balance`**: Events reordered by canonical Balance score (highest to lowest)
- **`score`**: Events reordered by canonical Splice Score (highest to lowest)
- **`luc7`**: Events reordered by canonical LUC7 score (highest to lowest)

> **Note:** When clustering is enabled, the rank row shows original rank values in the new clustered order, allowing you to see how clustering changes the arrangement.

## Output Files

The script generates 4 files in the specified output directory:

### Heatmap Images

1. **`heatmap_longer.png`** - Visualization for exon lengthening events (ΔSJ > 0)
2. **`heatmap_shorter.png`** - Visualization for exon shortening events (ΔSJ < 0)

**Heatmap Layout:**
```
┌─────────────────────────────────────┐
│  Balance - Canonical                │  ← RdBu_r, -3 to 3
├─────────────────────────────────────┤
│  Balance - Cryptic                  │  ← RdBu_r, -3 to 3
├─────────────────────────────────────┤
│  Splice Score - Canonical           │  ← viridis, data range
├─────────────────────────────────────┤
│  Splice Score - Cryptic             │  ← viridis, data range
├─────────────────────────────────────┤
│  LUC7 - Canonical                   │  ← plasma, data range
├─────────────────────────────────────┤
│  LUC7 - Cryptic                     │  ← plasma, data range
├─────────────────────────────────────┤
│  Rank - Canonical                   │  ← Blues, 1 to max
├─────────────────────────────────────┤
│  Rank - Cryptic                     │  ← Blues, 1 to max
└─────────────────────────────────────┘
       Events (left to right) →
```

**Color Schemes:**
- **Balance:** RdBu_r (red-white-blue), range -3 to 3, discrete ticks
- **Splice Score:** viridis (yellow-green-blue), data-driven range
- **LUC7 Score:** plasma (purple-pink-yellow), data-driven range
- **Rank:** Blues (light → dark blue), continuous gradient

### Data Tables

3. **`scored_longer.tsv`** - Original longer events data with matched scores
4. **`scored_shorter.tsv`** - Original shorter events data with matched scores

**New columns added:**
- `canonical_9mer` - Extracted 9-mer from canonical sequence
- `cryptic_9mer` - Extracted 9-mer from cryptic sequence
- `canonical_Balance`, `canonical_score`, `canonical_LUC7` - Scores for canonical site
- `cryptic_Balance`, `cryptic_score`, `cryptic_LUC7` - Scores for cryptic site

## Usage Examples

### Basic Usage (Default - Rank Order)

```bash
python workflow/scripts/create_score_heatmap.py \
  --reference test_output_filtering/conor_handedness.tsv \
  --longer test_output_filtering/splice_grid_output/test2/splice_site_analysis_longer.tsv \
  --shorter test_output_filtering/splice_grid_output/test2/splice_site_analysis_shorter.tsv \
  --output test_output_filtering/splice_grid_output/test3/heatmap_analysis/
```

### Cluster by Canonical Balance

```bash
python workflow/scripts/create_score_heatmap.py \
  --reference test_output_filtering/conor_handedness.tsv \
  --longer test_output_filtering/splice_grid_output/test2/splice_site_analysis_longer.tsv \
  --shorter test_output_filtering/splice_grid_output/test2/splice_site_analysis_shorter.tsv \
  --output test_output_filtering/splice_grid_output/test3/heatmap_balance/ \
  --cluster balance
```

### Cluster by Canonical Splice Score

```bash
python workflow/scripts/create_score_heatmap.py \
  --reference test_output_filtering/conor_handedness.tsv \
  --longer test_output_filtering/splice_grid_output/test2/splice_site_analysis_longer.tsv \
  --shorter test_output_filtering/splice_grid_output/test2/splice_site_analysis_shorter.tsv \
  --output test_output_filtering/splice_grid_output/test3/heatmap_score/ \
  --cluster score
```

### Cluster by Canonical LUC7 Score

```bash
python workflow/scripts/create_score_heatmap.py \
  --reference test_output_filtering/conor_handedness.tsv \
  --longer test_output_filtering/splice_grid_output/test2/splice_site_analysis_longer.tsv \
  --shorter test_output_filtering/splice_grid_output/test2/splice_site_analysis_shorter.tsv \
  --output test_output_filtering/splice_grid_output/test3/heatmap_luc7/ \
  --cluster luc7
```

## Expected Output

When you run the script, you'll see:

```
============================================================
SPLICE SITE SCORE HEATMAP GENERATION
============================================================

Loading reference scores: conor_handedness.tsv
  Reference sequences: 262144
  Columns: [list of columns]

------------------------------------------------------------
PROCESSING LONGER EVENTS (ΔSJ > 0)
------------------------------------------------------------
  Loading analysis file: splice_site_analysis_longer.tsv
    Events: 232
  Extracting 9-mers...
    Canonical 9-mers extracted: 232/232
    Cryptic 9-mers extracted: 232/232
  Matching canonical 9-mers to reference...
    Canonical matched: 230/232
  Matching cryptic 9-mers to reference...
    Cryptic matched: 228/232
  Clustering by canonical LUC7 score
  Creating heatmap: heatmap_longer.png
    Saved heatmap with 232 events
  Saved scored data: scored_longer.tsv

[Similar output for shorter events]

============================================================
SUMMARY
============================================================
Output directory: heatmap_analysis/
Longer events: 232
Shorter events: 519

Generated files:
  - heatmap_longer.png
  - heatmap_shorter.png
  - scored_longer.tsv
  - scored_shorter.tsv
============================================================
```

## Interpreting Results

### Unmatched Sequences

Some 9-mers may not exist in the reference database. These are filled with:
- **Balance:** 0 (neutral)
- **Splice Score:** Mean score from matched events
- **LUC7:** 0 (neutral)

### Clustering vs Rank Order

- **No clustering:** Events ordered by |ΔSJ| magnitude (original rank)
- **With clustering:** Events reordered by selected canonical metric
  - Allows identification of score-based patterns
  - Rank row shows where high-|ΔSJ| events fall in the new order
  - Useful for finding correlations between scores and cryptic site usage

## Troubleshooting

### No 9-mers extracted
- **Cause:** Sequences too short or missing '|' delimiter
- **Solution:** Check that input sequences follow format "EXONBASES|INTRONBASES" with at least 3 exon and 6 intron bases

### Low match rate
- **Cause:** 9-mers not present in reference database
- **Solution:** This is expected for some sequences; script handles gracefully with neutral values

### Empty output
- **Cause:** Input files have no events or wrong format
- **Solution:** Verify input files exist and have correct columns

## Related Files

- **Script:** `workflow/scripts/create_score_heatmap.py`
- **Reference Data:** `test_output_filtering/conor_handedness.tsv`
- **Analysis Data:** `test_output_filtering/splice_grid_output/test2/splice_site_analysis_*.tsv`

## Citation

If using this analysis, please cite the reference database source and the splice site analysis pipeline.

---

**Last Updated:** 2026-02-09
**Script Version:** 1.0
