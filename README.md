# Snakemake Workflow – RNAseq

**Author:** Dylan Stermer  
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
│   │   ├── fastp.yml
│   │   ├── subread_featureCounts.yml
│   │   ├── qualimap.yml
│   │   └── ...
│   └── scripts/
│       ├── ExtractIntronsFromGtf.py
│       ├── CountsToExpressionMatrix.py
│       ├── leafcutter/
│       └── ...
├── config/
│   ├── config.yaml
│   ├── samples.tsv
│   └── STAR_Genome_List.tsv
├── results/
│   ├── Alignments/
│   ├── QC/
│   ├── featureCounts/
│   └── SplicingAnalysis/
├── logs/
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

### samples.tsv

Required columns: `sample`, `STARGenomeName`, `Strandedness`, `Aligner`, `R1`, `R2`

Optional columns: `StudyFirstAuthor`, `cell_type`, `Approach`, `Description`, `SRA_accession`, `Platform`, `R1_link`, `R2_link`

The `sample` column defines output file names. `STARGenomeName` specifies which genome to use from `STAR_Genome_List.tsv`. `R1` and `R2` provide paths to local FASTQ files, or leave blank and provide `SRA_accession` to download from SRA. `Strandedness` should be `U` (unstranded), `FR` (forward), or `RF` (reverse). `Aligner` should be `STAR` or `minimap2`.

Example:
```
sample	STARGenomeName	R1	R2	SRA_accession	Strandedness	Aligner
sample1	hg38	/path/to/R1.fq.gz	/path/to/R2.fq.gz		U	STAR
```

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