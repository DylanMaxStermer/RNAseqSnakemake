# Snakemake Workflow – RNAseq

**Author:** Ben Fair, Dylan Stermer  
**Affiliation:** University of Chicago, Genetics  
**Contact:** dylanstermer@uchicago.edu  
**Version:** 2.0
**Last Updated:** 2026-02-10

---

## Description

A comprehensive RNA-seq analysis Snakemake workflow covering alignment through splicing characterization. Handles multiple species simultaneously, as defined in `config/STAR_Genome_List.tsv` and `config/samples.tsv`.

**Core capabilities:**
- Genome indexing (STAR) and read alignment
- QC (fastp, MultiQC) and gene expression quantification (featureCounts)
- BigWig generation (unstranded and grouped)
- Splicing analysis via Leafcutter2: junction counting, PSI quantification, junction classification, differential splicing
- Poison event identification: poison exons, poison donor/acceptor junctions, poison skipped exons
- Productive alternative splicing identification: productive donor/acceptor junctions
- Exon characterization: MaxEntScan splice site scoring, phyloP/phastCons conservation, flanking CDS annotation, splice site distance distributions
- Reciprocal liftover: mapping exon coordinates from hg38 to orthologous positions in other species

---

## Repository Structure

```
project_root/
├── workflow/
│   ├── Snakefile                         # Main entry point; defines all targets
│   ├── rules/
│   │   ├── common.py                     # Shared helpers, wildcard setup, config parsing
│   │   ├── IndexGenome.smk               # STAR genome indexing
│   │   ├── PreprocessAndAlign.smk        # fastp trimming + STAR alignment
│   │   ├── QC.smk                        # MultiQC, read count summaries
│   │   ├── ExpressionAnalysis.smk        # featureCounts gene expression
│   │   ├── MakeBigwigs.smk               # BigWig generation (unstranded + grouped)
│   │   ├── SplicingLC2.smk               # Leafcutter2 junction counting, PSI,
│   │   │                                 #   junction classification, differential splicing,
│   │   │                                 #   poison/productive event identification
│   │   ├── ExonCharacteristic.smk        # MaxEntScan scoring, conservation analysis,
│   │   │                                 #   flanking CDS annotation, splice site
│   │   │                                 #   distance distributions, constitutive exons
│   │   └── LiftOver.smk                  # Reciprocal liftover (hg38 → other species)
│   ├── envs/                             # Conda environment YAML files
│   │   ├── leafcutter2.yml
│   │   ├── maxEnt.yml
│   │   ├── phastcons_analysis.yml
│   │   └── wiggletools.yaml
│   └── scripts/
│       ├── leafcutter2/                  # Leafcutter2 submodule
│       ├── ExonCharacteristics/          # Exon characterization scripts
│       │   ├── calculate_maxent_scores_fast.py
│       │   ├── extractSpliceSite.py
│       │   ├── find_cassette_exons_from_leafcutter_leafcutter2.py
│       │   ├── plot_phastC_phyloP_Exon_and_Intron.py
│       │   ├── standardize_bed_format.py
│       │   ├── standardize_constitutive_exons.py
│       │   ├── single_intron_cds_finder.py
│       │   └── annotate_constitutive_exon_phase.py
│       ├── FlankingCDS/
│       │   ├── flanking_CDS_poisonDonorAcceptor.py
│       │   └── flanking_CDS_productiveDonorAcceptor.py
│       ├── 04_reciprocal_liftover_workflow.py
│       ├── poisonExon_ID_leafcutter2.py
│       ├── productiveDonorAcceptor_ID_leafcutter2.py
│       ├── splice_site_distribution.py   # Poison vs. productive splice site distance plots
│       └── plot_maxent_score_cdfs.py     # Ad-hoc: CDF plots of MaxEntScan scores by exon type
├── config/
│   ├── config.yaml                       # Main config: paths, liftover targets, scratch dir
│   ├── samples.tsv                       # Sample metadata and FASTQ paths
│   ├── STAR_Genome_List.tsv              # Genome builds and reference paths
│   └── contrast_groupfiles/             # Sample group files for differential splicing
├── results/
│   ├── Alignments/                       # BAM files and indices
│   ├── QC/                              # MultiQC reports, read count summaries
│   ├── featureCounts/                   # Gene expression count matrices
│   ├── bigwigs/                         # Coverage tracks (unstranded + grouped)
│   ├── SplicingAnalysis/
│   │   ├── leafcutter/                  # Junction counts and PSI tables (BED format)
│   │   ├── leafcutter2/                 # Leafcutter2 cluster ratios and junction counts
│   │   ├── ClassifyJuncs/               # Leafcutter2 junction classifications
│   │   ├── ObservedJuncsAnnotations/    # Junction annotations with splice site scores
│   │   ├── differential_splicing_tidy/  # Differential splicing results per contrast
│   │   ├── PoisonEvent/
│   │   │   ├── Exon/PE/                 # Poison exon calls
│   │   │   ├── DonorAcceptor/           # Poison donor/acceptor junctions + flanking CDS
│   │   │   └── Skipped/                 # Poison skipped exon junctions
│   │   └── productive/
│   │       ├── Constitutive/            # Constitutive exons from single-intron clusters
│   │       └── ProductiveDonorAcceptor/ # Productive alt donor/acceptor + flanking CDS
│   ├── ExonCharacteristics/
│   │   ├── {GenomeName}/
│   │   │   ├── StandardFormat/          # Standardized exon TSVs per exon type
│   │   │   ├── SpliceSites/             # MaxEntScan-scored donor/acceptor BED files
│   │   │   └── phyloP_phastCons/        # Conservation score summaries and per-exon data
│   │   └── plots/                       # Splice site distance distribution plots
│   └── liftOver/                        # Reciprocal liftover results per exon type × target species
│       └── {GenomeName}/{exon_type}/{target}/
│           ├── one_way/                 # Standard liftover output
│           ├── reciprocal/              # Reciprocal-filtered liftover output
│           └── ortholog_mapping/        # 1-to-1 and 1-to-many ortholog tables
├── logs/                                # Per-rule log files
├── analysis/                            # Rmd notebooks and ad-hoc analysis outputs
├── snakemake_profiles/slurm/            # SLURM cluster submission profile
├── .gitignore
└── README.md
```

---

## Usage

### Step 1: Configure workflow

Edit `config/config.yaml` to set paths and parameters.

Create `config/samples.tsv` with your sample metadata 

Create `config/STAR_Genome_List.tsv` with genome information 

### Step 2: Install dependencies

Create and activate the conda environment:
```bash
conda env create -f rnaSeqAnalysis.yaml
conda activate rnaSeqAnalysis
```

### Step 3: Execute workflow

Test your configuration by performing a dry-run via:
```bash
snakemake -n
```

Execute the workflow locally via:
```bash
snakemake --cores $N
```

using `$N` cores or run it in a cluster environment via the included slurm snakemake profile:
```bash
snakemake --profile snakemake_profiles/slurm  --conda-frontend conda
```

---

## Configuration Files

### config/config.yaml

Key settings:

| Key | Description |
|-----|-------------|
| `samples` | Path to `samples.tsv` |
| `STAR_genomes` | Path to `STAR_Genome_List.tsv` |
| `GenomesPrefix` | Root directory for all reference genome files |
| `contrast_group_files_prefix` | Path prefix for differential splicing group files |
| `scratch` | Scratch directory for temporary large files |
| `liftover.source_assembly` | Assembly name of the input exon files (e.g. `hg38`) |
| `liftover.targets.{name}.forward_chain` | Chain file: source → target species |
| `liftover.targets.{name}.reverse_chain` | Chain file: target → source species (for reciprocal filtering) |

Available liftover targets (chain files must exist in `ChainFiles/`):
`mm39` (Mouse), `rheMac10` (Macaque), `galGal6` (Chicken), `monDom5` (Opossum), `oryCun2` (Rabbit), `rn7` (Rat)

---

## Input File Formats

### samples.tsv

Tab-separated file defining samples and sequencing data locations.

**Required columns:**
- `sample` - Unique sample identifier (used in output file names)
- `STARGenomeName` - Reference genome name (must match an entry in `STAR_Genome_List.tsv`)
- `Strandedness` - Library strandedness: `U` (unstranded), `FR` (forward/second-strand), or `RF` (reverse/first-strand)
- `Aligner` - Alignment tool: `STAR` or `minimap2`
- `R1` - Path to forward read FASTQ file (or leave blank if using `SRA_accession`)
- `R2` - Path to reverse read FASTQ file (or leave blank for single-end or SRA download)

**Optional columns:**
- `SRA_accession` - SRA accession number for automatic download (use instead of providing `R1`/`R2` paths)
- `StudyFirstAuthor` - First author of study
- `cell_type` - Cell type or tissue
- `Approach` - Experimental approach
- `Description` - Sample description
- `Platform` - Sequencing platform
- `R1_link` - URL to forward read file
- `R2_link` - URL to reverse read file

**Example:**
```tsv
sample	STARGenomeName	R1	R2	SRA_accession	Strandedness	Aligner	cell_type
HFF_DMSO_rep1	GRCh38	/data/reads/HFF_DMSO_1_R1.fq.gz	/data/reads/HFF_DMSO_1_R2.fq.gz		RF	STAR	fibroblast
HFF_SMG1i_rep1	GRCh38	/data/reads/HFF_SMG1i_1_R1.fq.gz	/data/reads/HFF_SMG1i_1_R2.fq.gz		RF	STAR	fibroblast
SRR_sample		SRR123456			U	STAR	
```

---

### Contrast Group Files

Tab-separated files defining sample groupings for differential analysis. Place these files in `config/contrast_groupfiles/`.

**Required columns:**
- `sample` - Sample name (must match `sample` column in `samples.tsv`)
- `contrast` - Group assignment for this sample

**File naming convention:** `<condition1>_vs_<condition2>.txt`

**Example:** `config/contrast_groupfiles/SMG1i_24h_vs_DMSO_0h.txt`
- 2 columns sample &	contrast
```tsv
SMG1i_HFF_24h_1	NMD_enriched
SMG1i_HFF_24h_2	NMD_enriched
SMG1i_HFF_24h_3	NMD_enriched
DMSO_HFF_0h_1	Control
DMSO_HFF_0h_2	Control
DMSO_HFF_0h_3	Control
```

**Notes:**
- Group names (e.g., `NMD_enriched`, `Control`) are arbitrary but must be consistent within each file
- Each contrast file defines one pairwise comparison
- Samples not listed in a contrast file are excluded from that comparison

### analysis

- This is where Rmd should be kept for each project.
- New scripts can also be made here

## Overview 

![Workflow DAG](config/rulegraph.svg)


### Setup Leafcutter Environment

#### Create and Configure Conda Environment
```bash
# Create environment with R and core dependencies
conda create -n leafcutter -c conda-forge \
  r-base=4.2.1 \
  r-rcpp=1.0.10 \
  r-rstan=2.21.7 \
  r-stanheaders=2.21.0-7 \
  r-devtools

# Activate environment
conda activate leafcutter
```

#### Install Required R Packages

Launch R and install packages not available through conda:
```r
# Launch R
R

# Install Bioconductor packages
install.packages("BiocManager")
BiocManager::install("Biobase")
BiocManager::install("DirichletMultinomial")

# Install specific versions of oompa packages
install.packages("remotes")

remotes::install_version(
  "oompaBase",
  version = "3.2.9",
  repos = "https://cloud.r-project.org",
  type = "source"
)

remotes::install_version(
  "oompaData",
  version = "3.1.4",
  repos = "https://cloud.r-project.org",
  type = "source"
)

remotes::install_version(
  "TailRank",
  version = "3.2.2",
  repos = "https://cloud.r-project.org",
  type = "source"
)

# Install leafcutter from GitHub
devtools::install_github("davidaknowles/leafcutter/leafcutter", ref = "psi_2019")
```

#### Configure Snakemake Rule for Environment Conflicts

If running on a cluster with system-wide R modules that conflict with your conda environment, modify the `leafcutter_ds_contrasts` rule in your Snakefile to explicitly set environment variables:

**Remove the `conda:` directive from the rule** and add these exports to the shell command:
```snakemake
rule leafcutter_ds_contrasts:
    # ... (input, output, params remain unchanged)
    shell:
        """
        # Force conda environment to take precedence over system modules
        export PATH=<your_conda_path>/envs/leafcutter/bin:$PATH
        export LD_LIBRARY_PATH=<your_conda_path>/envs/leafcutter/lib:$LD_LIBRARY_PATH
        export R_LIBS_USER=""
        export R_LIBS=""
        
        mkdir -p {output.outputdir}
        Rscript workflow/scripts/leafcutter/scripts/leafcutter_ds.R \
          -p {threads} \
          -o {output.outputdir}/leaf \
          {params.ExtraParams} \
          -i {params.MinGroupSize} \
          -g {params.MinGroupSize} \
          {input.numers} \
          {input.groupfile} &> {log}
        """
```

**Replace `<your_conda_path>` with your actual conda installation path** (e.g., `/project/yangili1/dylan_stermer/miniconda3`).

**Why these exports are necessary:**
- `PATH`: Ensures your conda R is found before system R installations
- `LD_LIBRARY_PATH`: Ensures correct C++ standard library version for compiled R packages
- `R_LIBS_USER` and `R_LIBS`: Prevents R from loading incompatible packages from system locations

**Note:** If you added a `.libPaths()` call to the leafcutter R script, remove it - the environment variables handle library paths correctly.
