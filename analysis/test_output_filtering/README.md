# Splice Site Analysis Pipeline

This directory contains tools and outputs for analyzing alternative splice site usage, including distribution analysis, grid visualizations, and score-based CDF comparisons.

---

## Table of Contents

1. [Pipeline Overview](#pipeline-overview)
2. [Input Files Required](#input-files-required)
3. [Scripts and Usage](#scripts-and-usage)
4. [Step-by-Step Workflow](#step-by-step-workflow)
5. [Output Files](#output-files)
6. [Dependencies](#dependencies)

---

## Pipeline Overview

The analysis pipeline consists of multiple stages, starting from junction detection:

```
Stage -1: Junction Detection (from Leafcutter2)
├─ Input: Leafcutter2 outputs (cluster_ratios, junction_classifications), GTF
├─ Script: productiveDonorAcceptor_ID_leafcutter2.py
└─ Output: BED files (altDonorAcceptor.bed) → Input to Stage 0

Stage 0: Flanking CDS Annotation
├─ Input: BED files (from Stage -1), Reference GTF
├─ Scripts: flanking_CDS_productiveDonorAcceptor.py
│           flanking_CDS_poisonDonorAcceptor.py
└─ Output: productive_da_flanking, poison_da_flanking → Input to Stages 1-2

Stage 1: Distribution Analysis
├─ Input: productive_da_flanking, poison_da_flanking
├─ Script: splice_site_distribution.py
└─ Output: Histograms and statistics → splice_region_analysis/

Stage 2: Grid Visualization & Sequence Analysis
├─ Input: productive_da_flanking, Reference.fa
├─ Script: splice_site_grid_analysis.py
└─ Output: Grid plots + TSV files → splice_grid_output/

Stage 3: Score-Based CDF Analysis
├─ Input: splice_site_analysis_*.tsv, conor_handedness.tsv
├─ Script: create_score_cdf.py
└─ Output: CDF plots + binned data → connor_handedness/
```

---

## Input Files Required

### 1. Productive and Poison Flanking Files

**Location:**
- `productive/ProductiveDonorAcceptor/productive_da_flanking`
- `poison/poison_da_flanking.tsv`

**Description:**
- Contains splice junction information with flanking CDS annotations
- Includes columns: `chrom`, `start`, `end`, `name`, `strand`, `junction_type`, `intron_coord`, `spanning_productive`, `donor_CDS`, `acceptor_CDS`, `shared_CDS`, `alternative_CDS`, etc.

**How they're created:**
- **Stage 0 (Required):** Created by FlankingCDS scripts (see below)
- Input to FlankingCDS scripts: BED files from main Snakemake pipeline
  - `all_chr_altDonorAcceptor.bed` (for productive events)
  - `poisonDonorAcceptor.bed` (for poison events)
- These BED files are generated from: BAM files → junction detection → filtering
- See [Stage 0: FlankingCDS Scripts](#stage-0-flankingcds-scripts-prerequisite) for detailed usage

### 2. Reference Genome FASTA

**Location:**
- Project-specific (e.g., `/project/yangili1/genomic_references/hg38.fa`)

**Description:**
- Standard reference genome in FASTA format
- Must be indexed (`.fai` file should exist)
- Used to extract sequences around splice sites

### 3. Splice Site Score Reference

**Location:**
- `test_output_filtering/conor_handedness.tsv`

**Description:**
- Reference database with 262,144 scored 9-mer splice site sequences
- Contains columns:
  - `SE.5SS.bed.seq`: 9-mer sequence
  - `Balance`: Handedness score (-3 to 3)
  - `SE.5SS.bed.score`: MaxEnt splice score
  - `LUC7_score`: LUC7 binding score

**Source:**
- Provided by collaborator (Conor)
- Pre-computed scores for all possible 9-mer donor splice sites

---

## Scripts and Usage

### Stage -1: Junction Detection from Leafcutter2

This stage identifies alternative donor/acceptor junctions between productive (coding) isoforms using Leafcutter2 output.

---

#### Script -1: `productiveDonorAcceptor_ID_leafcutter2.py`

**Purpose:** Identify alternative donor and acceptor junctions from Leafcutter2 data

**Location:** `/project/yangili1/dylan_stermer/GitHubTesting/PoisonID/scripts/productiveDonorAcceptor_ID_leafcutter2.py`

**Usage:**
```bash
python /project/yangili1/dylan_stermer/GitHubTesting/PoisonID/scripts/productiveDonorAcceptor_ID_leafcutter2.py \
  -g /project/yangili1/genomic_references/gencode.v38.annotation.gtf \
  -l /path/to/leafcutter2/junction_classifications.txt.gz \
  -c /path/to/leafcutter2/cluster_ratios.txt.gz \
  -o test_output_filtering/productive/all_chr \
  --chromosome chr1
```

**Arguments:**
- `-g, --gtf`: Reference GTF file (required)
- `-l, --leafcutter2-annotations`: Leafcutter2 junction classification file (required)
- `-c, --leafcutter-counts`: Leafcutter2 cluster_ratios file, gzipped (required)
- `-o, --output-prefix`: Output prefix for files (base directory and sample name) (required)
- `--chromosome`: Filter for specific chromosome (optional, for faster testing)

**What it does:**
1. **Loads CDS regions** from GTF file
2. **Parses Leafcutter2 junction classifications** to identify coding junctions
3. **Identifies alternative splicing events** where:
   - **altDonor**: Donor site (5' splice site) is alternative, acceptor site matches CDS boundary
   - **altAcceptor**: Donor site matches CDS boundary, acceptor site (3' splice site) is alternative
4. **Finds spanning productive junctions** that share one splice site
5. **Calculates frame preservation** between alternative and spanning productive junctions

**Input File Requirements:**

**Leafcutter2 junction_classifications.txt.gz:**
- Format: `chr:start:end:clu_X_strand annotation_type`
- Annotation types: `annotated`, `cryptic_threeprime`, `cryptic_fiveprime`, `unannotated`
- Generated by Leafcutter2 pipeline from BAM files

**Leafcutter2 cluster_ratios.txt.gz:**
- Contains junction usage counts per cluster
- Format: columns for junction coordinates, rows for samples
- Junction format: `chr:start:end:clu_X_strand:annotation`

**Outputs:**

Creates in `{output_dir}/ProductiveDonorAcceptor/`:

1. **`{sample_name}_altDonorAcceptor.bed`** - Main output BED file
   - Format: Standard BED + custom columns
   - Columns:
     - `chr`, `start`, `end`: Junction coordinates (0-based)
     - `name`: Gene name
     - `strand`: + or -
     - `junction_type`: altDonor or altAcceptor
     - `intron_coord`: Junction coordinate string (chr:start-end)
     - `spanning_productive`: Comma-separated productive junctions sharing one boundary
     - `frame_diffs`: Frame differences for each spanning productive junction
     - `frame_preserved`: Boolean flags for frame preservation
     - `donor_coord`: Donor position
     - `acceptor_coord`: Acceptor position

2. **`{sample_name}_colored_altDonorAcceptor.bed`** - BED file for IGV visualization
   - Same junctions plus associated spanning productive junctions as separate entries
   - Color-coded: altDonor = blue, altAcceptor = red, spanning productive = gray

3. **`{sample_name}_altDonorAcceptor_stats.txt`** - Summary statistics
   - Total counts of altDonor and altAcceptor events
   - Spanning productive junction statistics
   - Frame preservation analysis

**Example Output:**
```
# After running:
test_output_filtering/productive/ProductiveDonorAcceptor/
├── all_chr_altDonorAcceptor.bed          ← Input to flanking_CDS_productiveDonorAcceptor.py
├── all_chr_colored_altDonorAcceptor.bed
└── all_chr_altDonorAcceptor_stats.txt
```

**Key Concepts:**
- **Productive junction**: Junction between two coding exons (both exons in-frame)
- **altDonor**: Alternative donor site creates a longer or shorter exon at the 5' end
- **altAcceptor**: Alternative acceptor site creates a longer or shorter exon at the 3' end
- **Spanning productive**: The canonical productive junction that spans across the alternative splice site
- **Frame preservation**: Whether the alternative junction maintains the reading frame

---

### Stage 0: FlankingCDS Scripts

These scripts annotate splice junctions with flanking CDS information, creating the input files needed for downstream analysis.

---

#### Script 0A: `flanking_CDS_productiveDonorAcceptor.py`

**Purpose:** Identify flanking CDS features for productive alternative donor/acceptor junctions

**Location:** `workflow/scripts/FlankingCDS/flanking_CDS_productiveDonorAcceptor.py`

**Usage:**
```bash
python workflow/scripts/FlankingCDS/flanking_CDS_productiveDonorAcceptor.py \
  -g /project/yangili1/genomic_references/gencode.v38.annotation.gtf \
  -p test_output_filtering/productive/ProductiveDonorAcceptor/all_chr_altDonorAcceptor.bed \
  -o test_output_filtering/productive/ProductiveDonorAcceptor/productive_da_flanking \
  -f test_output_filtering/productive/ProductiveDonorAcceptor/flanking_cds.tsv \
  --tolerance 0 \
  --filter-mode closest
```

**Arguments:**
- `-g, --gtf`: Reference GTF file with gene annotations (required)
- `-p, --productive-bed`: Productive alternative donor/acceptor BED file (required)
- `-o, --output`: Output file for junctions with flanking CDS info (required)
- `-f, --flanking-output`: Output file for flanking CDS features only (required)
- `--tolerance`: Coordinate matching tolerance in bp (default: 0)
- `--chromosome`: Process only this chromosome (optional, for testing)
- `--filter-mode`: Spanning productive junction filtering mode (default: closest)
  - `closest`: Keep only the closest spanning productive junction
  - `threshold`: Keep junctions within distance threshold (500 bp)
  - `none`: Keep all spanning productive junctions

**Input BED File Requirements:**
- Must contain columns: `chrom`, `start`, `end`, `name`, `strand`, `junction_type`, `spanning_productive`
- `junction_type` should be `altDonor` or `altAcceptor`
- `spanning_productive`: Comma-separated list of productive junction coordinates

**Outputs:**
1. **Main output** (`productive_da_flanking`):
   - All input columns plus:
   - `donor_CDS`: CDS region at donor splice site (format: chr:start-end)
   - `acceptor_CDS`: CDS region at acceptor splice site
   - `shared_CDS`: Which CDS is shared (unchanged) between spanning productive and alternative junction
   - `alternative_CDS`: Which CDS is alternative (changed)
   - For altDonor: acceptor_CDS is shared, donor_CDS is alternative
   - For altAcceptor: donor_CDS is shared, acceptor_CDS is alternative

2. **Flanking CDS features** (`flanking_cds.tsv`):
   - Subset of events with valid flanking CDS annotations

**How it works:**
1. Loads CDS features from GTF file
2. For each alternative donor/acceptor junction:
   - Identifies the spanning productive junction(s)
   - Finds CDS at the donor position (5' splice site)
   - Finds CDS at the acceptor position (3' splice site)
   - Annotates which CDS is shared vs alternative based on junction type

---

#### Script 0B: `flanking_CDS_poisonDonorAcceptor.py`

**Purpose:** Identify flanking CDS features for poison donor/acceptor junctions

**Location:** `workflow/scripts/FlankingCDS/flanking_CDS_poisonDonorAcceptor.py`

**Usage:**
```bash
python workflow/scripts/FlankingCDS/flanking_CDS_poisonDonorAcceptor.py \
  -g /project/yangili1/genomic_references/gencode.v38.annotation.gtf \
  -p path/to/poisonDonorAcceptor.bed \
  -o test_output_filtering/poison/poison_da_flanking.tsv \
  -f test_output_filtering/poison/flanking_cds.tsv \
  --tolerance 0 \
  --filter-mode closest
```

**Arguments:**
- `-g, --gtf`: Reference GTF file with gene annotations (required)
- `-p, --poison-bed`: Poison donor/acceptor BED file (required)
- `-o, --output`: Output file for junctions with flanking CDS info (required)
- `-f, --flanking-output`: Output file for flanking CDS features only (required)
- `--tolerance`: Coordinate matching tolerance in bp (default: 0)
- `--chromosome`: Process only this chromosome (optional, for testing)
- `--filter-mode`: Productive junction filtering mode (default: closest)

**Input BED File Requirements:**
- Must contain columns: `chrom`, `start`, `end`, `name`, `strand`, `junction_type`, `productive_junctions`
- `junction_type` should be `poisonDonor` or `poisonAcceptor`
- `productive_junctions`: Comma-separated list of associated productive junction coordinates

**Outputs:**
1. **Main output** (`poison_da_flanking.tsv`):
   - All input columns plus:
   - `donor_CDS`: CDS region at donor splice site
   - `acceptor_CDS`: CDS region at acceptor splice site
   - `shared_CDS`: Which CDS is shared between productive and poison junction
   - `alternative_CDS`: Which CDS is alternative
   - For poisonDonor: acceptor_CDS is shared, donor_CDS is alternative
   - For poisonAcceptor: donor_CDS is shared, acceptor_CDS is alternative

2. **Flanking CDS features** (`flanking_cds.tsv`):
   - Subset of events with valid flanking CDS annotations

**Notes:**
- Both scripts use the same core logic but expect different input BED formats
- Productive script expects `spanning_productive` column
- Poison script expects `productive_junctions` column
- GTF coordinates are 1-based, BED coordinates are 0-based (scripts handle conversion)

---

### Script 1: `splice_site_distribution.py`

**Purpose:** Analyze distribution of splice site distances (exon lengthening vs shortening)

**Location:** `workflow/scripts/splice_site_distribution.py`

**Usage:**
```bash
python workflow/scripts/splice_site_distribution.py \
  --poison test_output_filtering/poison/poison_da_flanking.tsv \
  --productive test_output_filtering/productive/ProductiveDonorAcceptor/productive_da_flanking \
  --output test_output_filtering/splice_region_analysis/ \
  --chromosome All \
  --intron-range 200 \
  --exon-range 200
```

**Arguments:**
- `--poison`: Path to poison events file (required)
- `--productive`: Path to productive events file (required)
- `--output`: Output directory (required)
- `--chromosome`: Filter to specific chromosome (default: All)
- `--intron-range`: X-axis limit for shortening events in bp (default: 200)
- `--exon-range`: X-axis limit for elongation events in bp (default: 200)

**Outputs:**
Creates in `splice_region_analysis/`:
- `donor/` - Donor splice site analysis
  - `counts/` - Count histograms
  - `proportion/` - Proportion histograms
- `acceptor/` - Acceptor splice site analysis
  - `counts/` - Count histograms
  - `proportion/` - Proportion histograms

---

### Script 2: `splice_site_grid_analysis.py`

**Purpose:** Create grid visualizations of ranked splice events and extract sequences

**Location:** `workflow/scripts/splice_site_grid_analysis.py`

**Usage (Single View - Longer Events):**
```bash
python workflow/scripts/splice_site_grid_analysis.py \
  --productive test_output_filtering/productive/ProductiveDonorAcceptor/productive_da_flanking \
  --reference /project/yangili1/genomic_references/hg38.fa \
  --output test_output_filtering/splice_grid_output/longer_only/ \
  --junction-type altDonor \
  --exon-change longer \
  --sequence-range 15 \
  --chromosome All
```

**Usage (Single View - Shorter Events):**
```bash
python workflow/scripts/splice_site_grid_analysis.py \
  --productive test_output_filtering/productive/ProductiveDonorAcceptor/productive_da_flanking \
  --reference /project/yangili1/genomic_references/hg38.fa \
  --output test_output_filtering/splice_grid_output/shorter_only/ \
  --junction-type altDonor \
  --exon-change shorter \
  --sequence-range 15 \
  --chromosome All
```

**Usage (Merged View - Both Longer and Shorter):**
```bash
python workflow/scripts/splice_site_grid_analysis.py \
  --productive test_output_filtering/productive/ProductiveDonorAcceptor/productive_da_flanking \
  --reference /project/yangili1/genomic_references/hg38.fa \
  --output test_output_filtering/splice_grid_output/merged/ \
  --junction-type altDonor \
  --merge \
  --sequence-range 15 \
  --chromosome All
```

**Arguments:**
- `--productive`: Path to productive_da_flanking file (required)
- `--reference`: Path to reference genome FASTA (required)
- `--output`: Output directory (required)
- `--junction-type`: `altDonor` or `altAcceptor` (default: altDonor)
- `--exon-change`: `longer` or `shorter` (default: longer) - ignored if --merge is used
- `--merge`: Plot both longer and shorter on same plot (optional)
- `--sequence-range`: bp flanking splice site for padding (default: 15)
- `--chromosome`: Filter to specific chromosome (default: All)
- `--IntronExon-range`: Max absolute cryptic_length to include (optional)

**Outputs:**
Creates in output directory:
- `splice_site_grid.png` - Grid visualization (or `splice_site_grid_merged.png` if --merge)
- `splice_site_analysis.tsv` - All events with sequences
- `splice_site_analysis_longer.tsv` - Only exon lengthening events (cryptic_length > 0)
- `splice_site_analysis_shorter.tsv` - Only exon shortening events (cryptic_length < 0)

**Key Output Columns (TSV files):**
- `cryptic_length`: Distance between cryptic and canonical splice site (bp)
- `rank`: Event ranking by |cryptic_length| (1 = largest)
- `canonical_sequence`: Canonical splice site sequence (format: "EXONBASES|INTRONBASES")
- `cryptic_sequence`: Cryptic splice site sequence (format: "EXONBASES|INTRONBASES")
- `junction_type`: altDonor or altAcceptor
- Genomic coordinates and gene information

---

### Script 3: `create_score_cdf.py`

**Purpose:** Generate binned CDF plots comparing canonical vs cryptic splice site scores

**Location:** `workflow/scripts/create_score_cdf.py`

**Usage (Quantile-Based Binning - Default):**
```bash
python workflow/scripts/create_score_cdf.py \
  --reference test_output_filtering/conor_handedness.tsv \
  --longer test_output_filtering/splice_grid_output/test2/splice_site_analysis_longer.tsv \
  --shorter test_output_filtering/splice_grid_output/test2/splice_site_analysis_shorter.tsv \
  --output test_output_filtering/connor_handedness/cdf_quantile/ \
  --bin 4
```

**Usage (Rank-Based Binning):**
```bash
python workflow/scripts/create_score_cdf.py \
  --reference test_output_filtering/conor_handedness.tsv \
  --longer test_output_filtering/splice_grid_output/test2/splice_site_analysis_longer.tsv \
  --shorter test_output_filtering/splice_grid_output/test2/splice_site_analysis_shorter.tsv \
  --output test_output_filtering/connor_handedness/cdf_rank/ \
  --bin 4 \
  --no-quantile
```

**Arguments:**
- `--reference`: Path to conor_handedness.tsv reference scores (required)
- `--longer`: Path to splice_site_analysis_longer.tsv (required)
- `--shorter`: Path to splice_site_analysis_shorter.tsv (required)
- `--output`: Output directory (required)
- `--bin`: Number of bins (default: 4)
- `--quantile` / `--no-quantile`: Binning method (default: quantile-based)
- `--scores`: Score columns to visualize (default: Balance, SE.5SS.bed.score, LUC7_score)

**Binning Methods:**

| Method | Flag | Behavior | N= per Bin | Use Case |
|--------|------|----------|------------|----------|
| **Quantile** | `--quantile` (default) | Bins by cryptic_length VALUES | ~Equal | Events with same length stay together |
| **Rank** | `--no-quantile` | Bins by RANK intervals | Exactly equal | Top 25%, second 25%, etc. |

**Key Difference:**
- **Quantile binning**: All events with the same `cryptic_length` stay in the same bin
- **Rank binning**: Events with the same `cryptic_length` can be split across bins

**Outputs:**
Creates in output directory:
- `cdf_plots.png` - 2×3 grid of CDF plots
  - Columns: Shortening events | Lengthening events
  - Rows: Balance | LUC7 Score | Splice Score
  - Line styles: Solid = Canonical, Dashed = Cryptic
  - Colors: Blue (Bin 1), Red (Bin 2), Green (Bin 3), Purple (Bin 4)
- `binned_longer.tsv` - Lengthening events with bin assignments
- `binned_shorter.tsv` - Shortening events with bin assignments

**Legend Format:**
```
Bin 1 (2 to 15 nt, n=58) - Canonical
Bin 1 (2 to 15 nt, n=58) - Cryptic
```
- Shows bin number, length range, event count, and site type

---

## Step-by-Step Workflow

### Complete Analysis Example (All Stages)

#### Stage -1: Identify Alternative Donor/Acceptor Junctions (from Leafcutter2)
```bash
# Identify productive alternative donor/acceptor junctions
python /project/yangili1/dylan_stermer/GitHubTesting/PoisonID/scripts/productiveDonorAcceptor_ID_leafcutter2.py \
  -g /project/yangili1/genomic_references/gencode.v38.annotation.gtf \
  -l /path/to/leafcutter2/junction_classifications.txt.gz \
  -c /path/to/leafcutter2/cluster_ratios.txt.gz \
  -o test_output_filtering/productive/all_chr

# Output: test_output_filtering/productive/ProductiveDonorAcceptor/all_chr_altDonorAcceptor.bed
```

#### Stage 0: Annotate Flanking CDS Features
```bash
# Productive events
python workflow/scripts/FlankingCDS/flanking_CDS_productiveDonorAcceptor.py \
  -g /project/yangili1/genomic_references/gencode.v38.annotation.gtf \
  -p test_output_filtering/productive/ProductiveDonorAcceptor/all_chr_altDonorAcceptor.bed \
  -o test_output_filtering/productive/ProductiveDonorAcceptor/productive_da_flanking \
  -f test_output_filtering/productive/ProductiveDonorAcceptor/flanking_cds.tsv \
  --filter-mode closest

# Poison events (if available)
python workflow/scripts/FlankingCDS/flanking_CDS_poisonDonorAcceptor.py \
  -g /project/yangili1/genomic_references/gencode.v38.annotation.gtf \
  -p /path/to/poisonDonorAcceptor.bed \
  -o test_output_filtering/poison/poison_da_flanking.tsv \
  -f test_output_filtering/poison/flanking_cds.tsv \
  --filter-mode closest

# Outputs:
#   - productive_da_flanking (required for Stages 1-2)
#   - poison_da_flanking.tsv (required for Stage 1)
```

#### Stage 1: Distribution Analysis (Optional)
```bash
python workflow/scripts/splice_site_distribution.py \
  --poison test_output_filtering/poison/poison_da_flanking.tsv \
  --productive test_output_filtering/productive/ProductiveDonorAcceptor/productive_da_flanking \
  --output test_output_filtering/splice_region_analysis/
```

#### Stage 2: Grid Analysis and Sequence Extraction
```bash
# Create merged view with both longer and shorter events
python workflow/scripts/splice_site_grid_analysis.py \
  --productive test_output_filtering/productive/ProductiveDonorAcceptor/productive_da_flanking \
  --reference /project/yangili1/genomic_references/hg38.fa \
  --output test_output_filtering/splice_grid_output/analysis_run1/ \
  --junction-type altDonor \
  --merge \
  --sequence-range 15
```

This creates:
- `splice_site_analysis_longer.tsv`
- `splice_site_analysis_shorter.tsv`
- Grid visualization PNG files

#### Stage 3: Score-Based CDF Analysis
```bash
python workflow/scripts/create_score_cdf.py \
  --reference test_output_filtering/conor_handedness.tsv \
  --longer test_output_filtering/splice_grid_output/analysis_run1/splice_site_analysis_longer.tsv \
  --shorter test_output_filtering/splice_grid_output/analysis_run1/splice_site_analysis_shorter.tsv \
  --output test_output_filtering/connor_handedness/analysis_run1/ \
  --bin 4
```

This creates:
- `cdf_plots.png` with 6 subplots (3 metrics × 2 event types)
- Binned data files with score annotations

---

## Output Files

### Directory Structure

```
test_output_filtering/
├── splice_region_analysis/          # From splice_site_distribution.py
│   ├── donor/
│   │   ├── counts/                  # Count histograms
│   │   └── proportion/              # Proportion histograms
│   └── acceptor/
│       ├── counts/
│       └── proportion/
│
├── splice_grid_output/              # From splice_site_grid_analysis.py
│   └── [run_name]/
│       ├── splice_site_grid_merged.png
│       ├── splice_site_analysis.tsv
│       ├── splice_site_analysis_longer.tsv
│       └── splice_site_analysis_shorter.tsv
│
└── connor_handedness/               # From create_score_cdf.py
    └── [run_name]/
        ├── cdf_plots.png
        ├── binned_longer.tsv
        └── binned_shorter.tsv
```

---

## Dependencies

### Python Packages

Create conda environment from `workflow/envs/`:
```bash
# For splice_site_grid_analysis.py (sequence extraction)
conda activate maxEnt  # or relevant environment with pyfaidx
pip install pyfaidx

# For all scripts
pip install pandas numpy matplotlib seaborn
```

### Required Files

1. **Productive/Poison flanking files** - From upstream Snakemake pipeline
2. **Reference genome FASTA** - Indexed with samtools faidx
3. **Splice site score reference** - conor_handedness.tsv

---

## Troubleshooting

### Common Issues

**Issue:** "pyfaidx not found"
```bash
conda activate maxEnt
pip install pyfaidx
```

**Issue:** "Reference genome not indexed"
```bash
samtools faidx /path/to/reference.fa
```

**Issue:** "No events found after filtering"
- Check chromosome names match between files
- Verify junction_type is present in data (altDonor vs altAcceptor)
- Try removing filters (--chromosome All, no --IntronExon-range)

**Issue:** "Duplicate values in quantile binning"
```bash
# Use rank-based binning instead
python workflow/scripts/create_score_cdf.py ... --no-quantile
```

---

## Citation

If using this analysis pipeline, please cite:
- MaxEnt splice site scoring: Yeo & Burge (2004)
- LUC7 scoring: [Source reference]
- Handedness/Balance metrics: Conor et al. (pending)

---

**Last Updated:** 2026-02-09
**Pipeline Version:** 1.0
**Contact:** Dylan Stermer
