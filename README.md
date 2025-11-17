# Snakemake Workflow – RNAseq

**Author:** Ben Fair, Dylan Stermer  
**Affiliation:** University of Chicago, Genetics  
**Contact:** dylanstermer@uchicago.edu  
**Version:** 1.0  
**Last Updated:** 2025-11-08

---

## Description

This workflow contains rules to download genome files, index genomes and align reads with STAR, perform basic QC, count splice junction reads (regtools) and gene reads (featureCounts). Can handle different samples from different species, as defined in `config/STAR_Genome_List.tsv` and `config/samples.tsv`. Because this is often just the start of an RNA-seq analysis, this workflow might be best used as a module in a Snakemake workflow that further extends this work.

---

## Repository Structure
```
project_root/
├── workflow/
│   ├── Snakefile
│   ├── rules/
│   │   ├── common.py
│   │   ├── IndexGenome.smk
│   │   ├── PreprocessAndAlign.smk
│   │   ├── QC.smk
│   │   ├── ExpressionAnalysis.smk
│   │   └── SplicingAnalysis.smk
│   ├── envs/
│   └── scripts/
│       ├── ExtractIntronsFromGtf.py
│       ├── CountsToExpressionMatrix.py
│       ├── leafcutter/
│       └── ...
├── config/
│   ├── config.yaml
│   ├── samples.tsv
│   ├── contrast_groupfiles
│   └── STAR_Genome_List.tsv
├── results/
│   ├── Alignments/
│   ├── QC/
│   ├── featureCounts/
│   ├── SplicingAnalysis/
│   └── ...
├── logs/
├── analysis/
│   └── plots
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
```tsv
sample	contrast
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