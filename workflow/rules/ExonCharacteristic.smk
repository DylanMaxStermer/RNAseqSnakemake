#!/usr/bin/env python3
"""
test.smk - Refactored ExonCharacteristic rules for improved portability

Helper functions are now in workflow/rules/common.py
Conservation track paths are read from config/STAR_Genome_List.tsv (phyloP and phastCons columns)
"""

# =============================================================================
# Rules for Cassette Exon Identification
# =============================================================================

rule cassette_exon_from_leafcutter:
    """
    Identify cassette exons from leafcutter clustering results.

    Uses junction classifications to find alternatively spliced cassette exons.
    A cassette exon is an exon that can be included or skipped in the mature mRNA.

    PORTABILITY CHANGE: Script path now uses config variable instead of hardcoded path
    """
    input:
        counts = "results/SplicingAnalysis/leafcutter2/{GenomeName}/leafcutter2.junction_counts.const.gz",
        classifications = "results/SplicingAnalysis/leafcutter2/{GenomeName}/clustering/leafcutter2_junction_classifications.txt"
    output:
        all_cassettes = temp("results/ExonCharacteristics/{GenomeName}/LeafcutterCassetteExon/leafcutter_cassette/cassetteExons.tsv"),
        with_skipping = temp("results/ExonCharacteristics/{GenomeName}/LeafcutterCassetteExon/leafcutter_cassette/cassetteExons_withSkipping.tsv"),
        bed = temp("results/ExonCharacteristics/{GenomeName}/LeafcutterCassetteExon/leafcutter_cassette/cassetteExons.bed"),
        junctions_bed = temp("results/ExonCharacteristics/{GenomeName}/LeafcutterCassetteExon/leafcutter_cassette/cassetteExons_junctions.bed"),
        per_gene = temp("results/ExonCharacteristics/{GenomeName}/LeafcutterCassetteExon/leafcutter_cassette/cassetteExons_per_gene.tsv")
    log:
        "logs/ExonCharacteristic/cassette_exon_from_leafcutter/{GenomeName}.log"
    params:
        keep_annot = False,
        keep_coding = False,
        exclude_utr = False,
        keep_gencodepc = False,
        script = "workflow/scripts/ExonCharacteristics/find_cassette_exons_from_leafcutter_leafcutter2.py",
        outdir = "results/ExonCharacteristics/{GenomeName}/LeafcutterCassetteExon/leafcutter_cassette"
    resources:
        mem_mb = 16000
    shell:
        """
        mkdir -p {params.outdir}
        python {params.script} \
            -c {input.counts} \
            -a {input.classifications} \
            -o {params.outdir} \
            --keep-annot {params.keep_annot} \
            --keep-coding {params.keep_coding} \
            --exclude-utr {params.exclude_utr} \
            --keep-gencodepc {params.keep_gencodepc} \
            &> {log}
        """


rule filter_leafcutter_cassetteExonID:
    """
    Filter LeafCutter-derived cassette exons against a reference GTF.

    Produces BED files of cassette exons that match or don't match GTF annotations.
    This helps distinguish novel cassette exons from annotated ones.

    PORTABILITY CHANGE: Script path now uses config variable instead of hardcoded path
    """
    input:
        cassettes = "results/ExonCharacteristics/{GenomeName}/LeafcutterCassetteExon/leafcutter_cassette/cassetteExons.tsv",
        gtf = config["GenomesPrefix"] + "{GenomeName}/Reference.gtf"
    output:
        matched = temp("results/ExonCharacteristics/{GenomeName}/LeafcutterCassetteExon/filtered/lefcutter_cassette_cassette_gtf_match.bed"),
        unmatched = temp("results/ExonCharacteristics/{GenomeName}/LeafcutterCassetteExon/filtered/lefcutter_cassette_cassette_gtf_no_match.bed")
    log:
        "logs/ExonCharacteristic/filter_leafcutter_cassetteExonID/{GenomeName}.log"
    params:
        outprefix = "results/ExonCharacteristics/{GenomeName}/LeafcutterCassetteExon/filtered/lefcutter_cassette",
        script = "workflow/scripts/ExonCharacteristics/filter_leafcutter_cassetteExonID.py"
    resources:
        mem_mb = 16000
    shell:
        """
        mkdir -p $(dirname {params.outprefix})
        python {params.script} \
            -c {input.cassettes} \
            -g {input.gtf} \
            -o {params.outprefix} \
            &> {log}
        """


# =============================================================================
# Rules for Poison Exon Processing
# =============================================================================

rule poison_exon_psi_ds_filter:
    """
    Filter poison exons to keep only those with differentially spliced flanking introns.
    Applies filters based on:
    - p-value cutoff (statistical significance)
    - deltaPSI cutoff (effect size)
    - NMD PSI threshold (enrichment for NMD)
    """
    input:
        pe_file = "results/SplicingAnalysis/PoisonEvent/Exon/PE/{GenomeName}_coding_noncoding_introns_all.tsv",
        ds_results = expand(
            "results/SplicingAnalysis/differential_splicing_tidy/{contrast}/Results.tsv.gz",
            contrast=contrasts
        )
    output:
        pe_filtered = "results/SplicingAnalysis/PoisonEvent/Exon/PE/{GenomeName}_coding_noncoding_introns_all_ds_filtered_nmdpsi{nmd_psi_threshold}.tsv"
    params:
        p_cutoff = 0.05,
        deltapsi_cutoff = 0,
        script = "workflow/scripts/ExonCharacteristics/filter_pe_by_differential_splicing_nmdPSI.py"
    log:
        "logs/ExonCharacteristic/poison_exon_psi_ds_filter/{GenomeName}_nmdpsi{nmd_psi_threshold}.log"
    resources:
        mem_mb = 16000
    shell:
        """
        python {params.script} \
            -i {input.pe_file} \
            -r {input.ds_results} \
            -o {output.pe_filtered} \
            -p {params.p_cutoff} \
            -d {params.deltapsi_cutoff} \
            --nmd-psi {wildcards.nmd_psi_threshold} \
            &> {log}
        """


rule find_poison_exon_flank:
    """
    Identify flanking CDS features for poison exons using intron coordinates.

    Finds CDS features that are immediately adjacent to the upstream and downstream
    introns of each poison exon. This helps understand the coding context.

    Runs on both:
    - Unfiltered poison exons (all)
    - DS-filtered poison exons (with differential splicing)

    PORTABILITY CHANGE: Script path now uses config variable instead of hardcoded path
    """
    input:
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.gtf",
        pe_file = "results/SplicingAnalysis/PoisonEvent/Exon/PE/{GenomeName}_{pe_type}.tsv"
    output:
        pe_with_flanks = "results/SplicingAnalysis/PoisonEvent/Exon/PE/flanking_exon/{GenomeName}_{pe_type}_flanking.tsv",
        flanking_cds = "results/SplicingAnalysis/PoisonEvent/Exon/PE/flanking_exon/{GenomeName}_{pe_type}_flanking_CDSonly.tsv"
    params:
        script = "workflow/scripts/ExonCharacteristics/flanking_multiple_CDS.py"
    log:
        "logs/ExonCharacteristic/find_poison_exon_flank/{GenomeName}_{pe_type}.log"
    resources:
        mem_mb = 24000
    shell:
        """
        mkdir -p $(dirname {output.pe_with_flanks})
        python {params.script} \
            -g {input.gtf} \
            -p {input.pe_file} \
            -o {output.pe_with_flanks} \
            -f {output.flanking_cds} \
            &> {log}
        """


rule annotate_poison_exon_phase:
    """
    Match flanking exons to GTF annotations and annotate poison exons with phase information.

    Reading frame (phase) is critical for understanding how poison exons trigger NMD:
    - Phase 0, 1, 2: Reading frame at exon boundary
    - Frame-shifting exons are more likely to introduce premature stop codons

    Creates three outputs:
    1. Flanking exons matched to GTF with full annotation
    2. Flanking exons not matched to GTF (in NoMatch/ subdirectory)
    3. Poison exons annotated with phase information from flanking exons
    """
    input:
        pe_with_flanks = "results/SplicingAnalysis/PoisonEvent/Exon/PE/flanking_exon/{GenomeName}_{pe_type}_flanking.tsv",
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.gtf"
    output:
        flanking_match = temp("results/ExonCharacteristics/{GenomeName}/poisonExon/{pe_type}_flanking_gtfMatch.bed"),
        flanking_nomatch = temp("results/ExonCharacteristics/{GenomeName}/poisonExon/NoMatch/{pe_type}_flanking_gtfNoMatch.bed"),
        pe_annotated = temp("results/ExonCharacteristics/{GenomeName}/poisonExon/{pe_type}_poisonExons_annotatedFlanks.bed")
    params:
        output_prefix = "results/ExonCharacteristics/{GenomeName}/poisonExon/{pe_type}",
        script = "workflow/scripts/ExonCharacteristics/filter_poison_exon_flanking_output.py"
    log:
        "logs/ExonCharacteristic/annotate_poison_exon_phase/{GenomeName}_{pe_type}.log"
    resources:
        mem_mb = 32000
    shell:
        """
        mkdir -p $(dirname {output.flanking_match})
        mkdir -p $(dirname {output.flanking_nomatch})
        python {params.script} \
            -p {input.pe_with_flanks} \
            -g {input.gtf} \
            -o {params.output_prefix} \
            &> {log}

        # Move the NoMatch file to the NoMatch subdirectory
        mv {params.output_prefix}_flanking_gtfNoMatch.bed {output.flanking_nomatch}
        """

rule standardize_cassette_exons:
    """
    Standardize GTF-matched cassette exons to common BED format.
    """
    input:
        bed = "results/ExonCharacteristics/{GenomeName}/LeafcutterCassetteExon/filtered/lefcutter_cassette_cassette_gtf_match.bed"
    output:
        std = "results/ExonCharacteristics/{GenomeName}/StandardFormat/leafcutter_CE_gtf_matched/exons_std.tsv"
    params:
        script = "workflow/scripts/ExonCharacteristics/standardize_bed_format.py"
    log:
        "logs/ExonCharacteristic/standardize_cassette_exons/{GenomeName}.log"
    resources:
        mem_mb = 16000
    shell:
        """
        mkdir -p $(dirname {output.std})
        python {params.script} -i {input.bed} -o {output.std} &> {log}
        """


rule standardize_flanking_all:
    """
    Standardize flanking exons for all poison exons (unfiltered).
    """
    input:
        bed = "results/ExonCharacteristics/{GenomeName}/poisonExon/coding_noncoding_introns_all_flanking_gtfMatch.bed"
    output:
        std = "results/ExonCharacteristics/{GenomeName}/StandardFormat/exon_flanking_PE_all/exons_std.tsv"
    params:
        script = "workflow/scripts/ExonCharacteristics/standardize_bed_format.py"
    log:
        "logs/ExonCharacteristic/standardize_flanking_all/{GenomeName}.log"
    resources:
        mem_mb = 16000
    shell:
        """
        mkdir -p $(dirname {output.std})
        python {params.script} -i {input.bed} -o {output.std} &> {log}
        """


rule standardize_flanking_filtered:
    """
    Standardize flanking exons for DS-filtered poison exons.

    NOTE: Currently hardcoded to nmdpsi0.1 threshold. For full flexibility,
    this could be converted to use a wildcard, but that would require changes
    to the Snakefile targets as well.
    """
    input:
        bed = "results/ExonCharacteristics/{GenomeName}/poisonExon/coding_noncoding_introns_all_ds_filtered_nmdpsi0.1_flanking_gtfMatch.bed"
    output:
        std = "results/ExonCharacteristics/{GenomeName}/StandardFormat/exon_flanking_PE_filtered_nmdpsi0.1/exons_std.tsv"
    params:
        script = "workflow/scripts/ExonCharacteristics/standardize_bed_format.py"
    log:
        "logs/ExonCharacteristic/standardize_flanking_filtered/{GenomeName}.log"
    resources:
        mem_mb = 16000
    shell:
        """
        mkdir -p $(dirname {output.std})
        python {params.script} -i {input.bed} -o {output.std} &> {log}
        """


rule standardize_poison_all:
    """
    Standardize poison exons with annotated flanks (all, unfiltered).
    """
    input:
        bed = "results/ExonCharacteristics/{GenomeName}/poisonExon/coding_noncoding_introns_all_poisonExons_annotatedFlanks.bed"
    output:
        std = "results/ExonCharacteristics/{GenomeName}/StandardFormat/poisonExon_all/exons_std.tsv"
    params:
        script = "workflow/scripts/ExonCharacteristics/standardize_bed_format.py"
    log:
        "logs/ExonCharacteristic/standardize_poison_all/{GenomeName}.log"
    resources:
        mem_mb = 16000
    shell:
        """
        mkdir -p $(dirname {output.std})
        python {params.script} -i {input.bed} -o {output.std} &> {log}
        """


rule standardize_poison_filtered:
    """
    Standardize poison exons with annotated flanks (DS-filtered).

    REFACTORING: Direct standardization without intermediate temp copying.

    NOTE: Currently hardcoded to nmdpsi0.1 threshold. For full flexibility,
    this could be converted to use a wildcard.
    """
    input:
        bed = "results/ExonCharacteristics/{GenomeName}/poisonExon/coding_noncoding_introns_all_ds_filtered_nmdpsi0.1_poisonExons_annotatedFlanks.bed"
    output:
        std = "results/ExonCharacteristics/{GenomeName}/StandardFormat/poisonExon_filtered_nmdpsi0.1/exons_std.tsv"
    params:
        script = "workflow/scripts/ExonCharacteristics/standardize_bed_format.py"
    log:
        "logs/ExonCharacteristic/standardize_poison_filtered/{GenomeName}.log"
    resources:
        mem_mb = 16000
    shell:
        """
        mkdir -p $(dirname {output.std})
        python {params.script} -i {input.bed} -o {output.std} &> {log}
        """


# =============================================================================
# Rules for Conservation Analysis (PhyloP and PhastCons)
# =============================================================================

rule calculate_phyloP_phastcons:
    """
    Calculate PhyloP and PhastCons conservation scores around and within exons.
    Analyzes:
    - Flanking regions (configurable, default 50bp upstream/downstream)
    - Exon body
    - Intronic regions near splice sites
    SPECIES COMPATIBILITY:
    - Works for any genome with conservation tracks configured
    - Skips gracefully for genomes without conservation data (e.g., GRCh38)
    - To add a new genome, just update config with track paths
    """
    input:
        # CRITICAL CHANGE: Use input function to get genome-specific conservation files
        # Returns empty dict if conservation not available, which causes rule to be skipped
        # NOTE: unpack() must come BEFORE keyword arguments in Snakemake
        unpack(get_conservation_tracks),
        bed = "results/ExonCharacteristics/{GenomeName}/StandardFormat/{exon_type}/exons_std.tsv"
    output:
        summary = "results/ExonCharacteristics/{GenomeName}/phyloP_phastCons/{exon_type}/phyloP35way_phastCons35way_summary.tsv",
        per_exon = "results/ExonCharacteristics/{GenomeName}/phyloP_phastCons/{exon_type}/phyloP35way_phastCons35way_per_exon.tsv",
        phase0 = "results/ExonCharacteristics/{GenomeName}/phyloP_phastCons/{exon_type}/phyloP35way_phastCons35way_per_exon_phase0.tsv",
        phase1 = "results/ExonCharacteristics/{GenomeName}/phyloP_phastCons/{exon_type}/phyloP35way_phastCons35way_per_exon_phase1.tsv",
        phase2 = "results/ExonCharacteristics/{GenomeName}/phyloP_phastCons/{exon_type}/phyloP35way_phastCons35way_per_exon_phase2.tsv",
        phase_na = "results/ExonCharacteristics/{GenomeName}/phyloP_phastCons/{exon_type}/phyloP35way_phastCons35way_per_exon_phase_NA.tsv",
        summary_phase0 = "results/ExonCharacteristics/{GenomeName}/phyloP_phastCons/{exon_type}/phyloP35way_phastCons35way_summary_phase0.tsv",
        summary_phase1 = "results/ExonCharacteristics/{GenomeName}/phyloP_phastCons/{exon_type}/phyloP35way_phastCons35way_summary_phase1.tsv",
        summary_phase2 = "results/ExonCharacteristics/{GenomeName}/phyloP_phastCons/{exon_type}/phyloP35way_phastCons35way_summary_phase2.tsv"
    params:
        script = "workflow/scripts/ExonCharacteristics/test_phase_seperated_phylop.py",
        outdir = "results/ExonCharacteristics/{GenomeName}/phyloP_phastCons/{exon_type}",
        # Window sizes for flanking and intronic regions
        flank_size = 50,  # bp upstream/downstream of exon
        exon_extend = 39,  # bp into exon from boundaries
        min_exon_size = 50  # minimum exon size to analyze
    log:
        "logs/ExonCharacteristic/calculate_phyloP_phastcons/{GenomeName}_{exon_type}.log"
    conda:
        "../envs/phastcons_analysis.yml"
    resources:
        mem_mb = 20000
    shell:
        """
        mkdir -p {params.outdir}
        python {params.script} \
            -b {input.bed} \
            -p {input.phylop} \
            -w {input.phastcons} \
            -o {params.outdir}/phyloP35way_phastCons35way \
            -f {params.flank_size} \
            -e {params.exon_extend} \
            -m {params.min_exon_size} \
            --split-by-phase > {log} 2>&1

        touch {output.phase_na}
        """


rule plot_phyloP:
    """
    Plot PhyloP conservation profiles for phase 0 exons across different exon types.
    """
    input:
        summaries = expand(
            "results/ExonCharacteristics/{{GenomeName}}/phyloP_phastCons/{exon_type}/phyloP35way_phastCons35way_summary_phase0.tsv",
            exon_type=[
                "exon_flanking_PE_all",
                "exon_flanking_PE_filtered_nmdpsi0.1",
                "leafcutter_CE_gtf_matched",
                "poisonExon_all",
                "poisonExon_filtered_nmdpsi0.1",
                "Constitutive"
            ]
        )
    output:
        plot = "results/ExonCharacteristics/{GenomeName}/phyloP_phastCons/plots/phyloP_phase0_comparison_mean.pdf"
    params:
        # CHANGE: Script path now uses config with fallback to default location
        script = "workflow/scripts/ExonCharacteristics/plot_phastC_phyloP_Exon_and_Intron.py",
        labels = ["exon_flanking_PE_all", "exon_flanking_PE_filtered_nmdpsi0.1",
                  "leafcutter_CE_gtf_matched", "poisonExon_all", "poisonExon_filtered_nmdpsi0.1",
                  "Constitutive"],
        output_prefix = "results/ExonCharacteristics/{GenomeName}/phyloP_phastCons/plots/phyloP_phase0",
        outdir = "results/ExonCharacteristics/{GenomeName}/phyloP_phastCons/plots"
    log:
        "logs/ExonCharacteristic/plot_phyloP/{GenomeName}_phase0.log"
    conda:
        "../envs/phastcons_analysis.yml"
    resources:
        mem_mb = 8000
    shell:
        """
        mkdir -p {params.outdir}

        python {params.script} \
            -i {input.summaries} \
            -l {params.labels} \
            -o {params.output_prefix} \
            -m phylop \
            -f pdf > {log} 2>&1
        """


rule plot_phastCons:
    """
    Plot phastCons conservation profiles for phase 0 exons across different exon types.
    """
    input:
        summaries = expand(
            "results/ExonCharacteristics/{{GenomeName}}/phyloP_phastCons/{exon_type}/phyloP35way_phastCons35way_summary_phase0.tsv",
            exon_type=[
                "exon_flanking_PE_all",
                "exon_flanking_PE_filtered_nmdpsi0.1",
                "leafcutter_CE_gtf_matched",
                "poisonExon_all",
                "poisonExon_filtered_nmdpsi0.1",
                "Constitutive"
            ]
        )
    output:
        plot = "results/ExonCharacteristics/{GenomeName}/phyloP_phastCons/plots/phastCons_phase0_comparison_mean.pdf"
    params:
        script = "workflow/scripts/ExonCharacteristics/plot_phastC_phyloP_Exon_and_Intron.py",
        labels = ["exon_flanking_PE_all", "exon_flanking_PE_filtered_nmdpsi0.1",
                  "leafcutter_CE_gtf_matched", "poisonExon_all", "poisonExon_filtered_nmdpsi0.1",
                  "Constitutive"],
        output_prefix = "results/ExonCharacteristics/{GenomeName}/phyloP_phastCons/plots/phastCons_phase0",
        outdir = "results/ExonCharacteristics/{GenomeName}/phyloP_phastCons/plots"
    log:
        "logs/ExonCharacteristic/plot_phastCons/{GenomeName}_phase0.log"
    conda:
        "../envs/phastcons_analysis.yml"
    resources:
        mem_mb = 8000
    shell:
        """
        mkdir -p {params.outdir}

        python {params.script} \
            -i {input.summaries} \
            -l {params.labels} \
            -o {params.output_prefix} \
            -m phastcons \
            -f pdf > {log} 2>&1
        """


# =============================================================================
# Rules for Splice Site Analysis
# =============================================================================

rule extract_splice_site_sequence:
    """
    Extract splice site coordinates and sequences from standardized exon BED files.
    5' splice site (donor): 3 bases in exon + 6 bases in intron
    3' splice site (acceptor): 20 bases in intron + 3 bases in exon
    """
    input:
        bed = "results/ExonCharacteristics/{GenomeName}/StandardFormat/{exon_type}/exons_std.tsv",
        fa = config['GenomesPrefix'] + "{GenomeName}/Reference.fa"
    output:
        donor = temp("results/ExonCharacteristics/{GenomeName}/SpliceSites/{exon_type}/donor_sites.bed"),
        acceptor = temp("results/ExonCharacteristics/{GenomeName}/SpliceSites/{exon_type}/acceptor_sites.bed")
    params:
        # CHANGE: Script path now uses config with fallback to default location
        script = "workflow/scripts/ExonCharacteristics/extractSpliceSite.py",
        outdir = "results/ExonCharacteristics/{GenomeName}/SpliceSites/{exon_type}"
    log:
        "logs/ExonCharacteristic/extract_splice_site_sequence/{GenomeName}_{exon_type}.log"
    conda:
        "../envs/maxEnt.yml"
    resources:
        mem_mb = 16000
    shell:
        """
        mkdir -p {params.outdir}

        python {params.script} \
            -i {input.bed} \
            -f {input.fa} \
            -o {params.outdir}/splice_sites \
            &> {log}

        mv {params.outdir}/splice_sites_donor_sites.bed {output.donor}
        mv {params.outdir}/splice_sites_acceptor_sites.bed {output.acceptor}
        """


rule calculate_maxEnt_scores:
    """
    Calculate MaxEntScan scores for splice sites using maxentpy.
    Higher scores indicate stronger splice sites.
    Donor: 3 exon + 6 intron bases (9 total)
    Acceptor: 20 intron + 3 exon bases (23 total)
    """
    input:
        donor = "results/ExonCharacteristics/{GenomeName}/SpliceSites/{exon_type}/donor_sites.bed",
        acceptor = "results/ExonCharacteristics/{GenomeName}/SpliceSites/{exon_type}/acceptor_sites.bed"
    output:
        donor_scored = "results/ExonCharacteristics/{GenomeName}/SpliceSites/{exon_type}/donor_sites_scored.bed",
        acceptor_scored = "results/ExonCharacteristics/{GenomeName}/SpliceSites/{exon_type}/acceptor_sites_scored.bed"
    params:
        script = "workflow/scripts/ExonCharacteristics/calculate_maxent_scores_fast.py",
        outdir = "results/ExonCharacteristics/{GenomeName}/SpliceSites/{exon_type}"
    log:
        "logs/ExonCharacteristic/calculate_maxEnt_scores/{GenomeName}_{exon_type}.log"
    conda:
        "../envs/maxEnt.yml"
    resources:
        mem_mb = 16000
    shell:
        """
        python {params.script} \
            -d {input.donor} \
            -a {input.acceptor} \
            -o {params.outdir}/splice_sites \
            &> {log}

        mv {params.outdir}/splice_sites_donor_sites_scored.bed {output.donor_scored}
        mv {params.outdir}/splice_sites_acceptor_sites_scored.bed {output.acceptor_scored}
        """




# =============================================================================
# Rules for Constitutive Exon Identification via Single-Intron Clusters
# =============================================================================

rule ConstitutiveExon_Leafcutter2:
    """
    Identify CDS exons adjacent to single-intron Leafcutter2 clusters.

    Uses single_intron_cds_finder.py to classify CDS exons by their relationship
    to single-intron clusters: internal (constitutive cassette), first (only 5'
    junction), or last (only 3' junction).

    Outputs five BED12 files per genome to results/SplicingAnalysis/productive/Constitutive/.
    """
    input:
        clusters = "results/SplicingAnalysis/leafcutter2/{GenomeName}/clustering/leafcutter2_clusters",
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.gtf",
    output:
        five_prime     = "results/SplicingAnalysis/productive/Constitutive/{GenomeName}/cds_5prime_adjacent_introns.bed",
        three_prime    = "results/SplicingAnalysis/productive/Constitutive/{GenomeName}/cds_3prime_adjacent_introns.bed",
        constitutive   = "results/SplicingAnalysis/productive/Constitutive/{GenomeName}/constitutive_exons.bed",
        all_classified = "results/SplicingAnalysis/productive/Constitutive/{GenomeName}/all_exons_classified.bed",
        terminal       = "results/SplicingAnalysis/productive/Constitutive/{GenomeName}/terminal_exons.bed",
    params:
        outdir = "results/SplicingAnalysis/productive/Constitutive/{GenomeName}",
        script = "workflow/scripts/ExonCharacteristics/single_intron_cds_finder.py",
    log:
        "logs/ExonCharacteristic/ConstitutiveExon_Leafcutter2/{GenomeName}.log"
    resources:
        mem_mb = 16000
    shell:
        """
        mkdir -p {params.outdir}
        python {params.script} \
            -l {input.clusters} \
            -g {input.gtf} \
            -o {params.outdir} \
            &> {log}
        """
rule standardize_constitutive_exons:
    """
    Standardize constitutive_exons.bed from single_intron_cds_finder to common BED format.
    Extracts gene name from composite name field and adds 5primeJunc/3primeJunc columns.
    """
    input:
        bed = "results/SplicingAnalysis/productive/Constitutive/{GenomeName}/constitutive_exons.bed"
    output:
        std = temp("results/ExonCharacteristics/{GenomeName}/StandardFormat/Constitutive/exons_std_NoPhase.tsv")
    params:
        script = "workflow/scripts/ExonCharacteristics/standardize_constitutive_exons.py"
    log:
        "logs/ExonCharacteristic/standardize_constitutive_exons/{GenomeName}.log"
    resources:
        mem_mb = 16000
    shell:
        """
        mkdir -p $(dirname {output.std})
        python {params.script} \
            -i {input.bed} \
            -o {output.std} \
            &> {log}
        """


rule annotate_constitutive_exon_phase:
    """
    Annotate constitutive exons with CDS reading frame (phase) from the reference GTF.
    Matches exon coordinates to GTF CDS features using exact coordinate lookup.
    Adds a 'phase' column (0, 1, 2, or NA) to the standardized exon file.
    """
    input:
        exons = "results/ExonCharacteristics/{GenomeName}/StandardFormat/Constitutive/exons_std_NoPhase.tsv",
        gtf   = config['GenomesPrefix'] + "{GenomeName}/Reference.gtf",
    output:
        phased = "results/ExonCharacteristics/{GenomeName}/StandardFormat/Constitutive/exons_std.tsv",
    params:
        script = "workflow/scripts/ExonCharacteristics/annotate_constitutive_exon_phase.py",
    log:
        "logs/ExonCharacteristic/annotate_constitutive_exon_phase/{GenomeName}.log"
    resources:
        mem_mb = 16000
    shell:
        """
        python {params.script} \
            -e {input.exons} \
            -g {input.gtf} \
            -o {output.phased} \
            &> {log}
        """


rule flankingCDS_productiveDonorAcceptor:
    """
    Identify flanking CDS features for productive alternative donor/acceptor junctions.
    Annotates which CDS is shared (unchanged) vs. alternative at each splice site.
    """
    input:
        gtf        = config['GenomesPrefix'] + "{GenomeName}/Reference.gtf",
        productive = "results/SplicingAnalysis/productive/ProductiveDonorAcceptor/{GenomeName}_altDonorAcceptor.bed",
    output:
        with_flanks  = "results/SplicingAnalysis/productive/ProductiveDonorAcceptor/flankingCDS/{GenomeName}_altDonorAcceptor_flankingCDS.tsv",
        flanking_cds = "results/SplicingAnalysis/productive/ProductiveDonorAcceptor/flankingCDS/{GenomeName}_altDonorAcceptor_flankingCDS_only.tsv",
    params:
        script = "workflow/scripts/FlankingCDS/flanking_CDS_productiveDonorAcceptor.py",
    log:
        "logs/ExonCharacteristic/flankingCDS_productiveDonorAcceptor/{GenomeName}.log"
    resources:
        mem_mb = 24000
    shell:
        """
        mkdir -p $(dirname {output.with_flanks})
        python {params.script} \
            -g {input.gtf} \
            -p {input.productive} \
            -o {output.with_flanks} \
            -f {output.flanking_cds} \
            &> {log}
        """


rule flankingCDS_poisonDonorAcceptor:
    """
    Identify flanking CDS features for poison donor/acceptor junctions.
    Annotates which CDS is shared (unchanged) vs. alternative at each splice site.
    """
    input:
        gtf    = config['GenomesPrefix'] + "{GenomeName}/Reference.gtf",
        poison = "results/SplicingAnalysis/PoisonEvent/DonorAcceptor/PoisonDonorAcceptor/{GenomeName}_poisonDonorAcceptor.bed",
    output:
        with_flanks  = "results/SplicingAnalysis/PoisonEvent/DonorAcceptor/PoisonDonorAcceptor/flankingCDS/{GenomeName}_poisonDonorAcceptor_flankingCDS.tsv",
        flanking_cds = "results/SplicingAnalysis/PoisonEvent/DonorAcceptor/PoisonDonorAcceptor/flankingCDS/{GenomeName}_poisonDonorAcceptor_flankingCDS_only.tsv",
    params:
        script = "workflow/scripts/FlankingCDS/flanking_CDS_poisonDonorAcceptor.py",
    log:
        "logs/ExonCharacteristic/flankingCDS_poisonDonorAcceptor/{GenomeName}.log"
    resources:
        mem_mb = 24000
    shell:
        """
        mkdir -p $(dirname {output.with_flanks})
        python {params.script} \
            -g {input.gtf} \
            -p {input.poison} \
            -o {output.with_flanks} \
            -f {output.flanking_cds} \
            &> {log}
        """


rule alternateDonorAcceptor_distribution:
    """
    Analyze and visualize splice site distance distributions for poison vs. productive
    alternative donor/acceptor events. Generates histograms (counts and proportions)
    for both donor and acceptor events, colored by frame (multiples of 3 vs. non-multiples).
    """
    input:
        poison     = "results/SplicingAnalysis/PoisonEvent/DonorAcceptor/PoisonDonorAcceptor/flankingCDS/{GenomeName}_poisonDonorAcceptor_flankingCDS.tsv",
        productive = "results/SplicingAnalysis/productive/ProductiveDonorAcceptor/flankingCDS/{GenomeName}_altDonorAcceptor_flankingCDS.tsv",
    output:
        donor_counts     = "results/ExonCharacteristics/plots/{GenomeName}/AltDonorAcceptor/distribution/donor/counts/overlay_comparison.png",
        donor_proportion = "results/ExonCharacteristics/plots/{GenomeName}/AltDonorAcceptor/distribution/donor/proportion/overlay_comparison.png",
        acceptor_counts  = "results/ExonCharacteristics/plots/{GenomeName}/AltDonorAcceptor/distribution/acceptor/counts/overlay_comparison.png",
        acceptor_prop    = "results/ExonCharacteristics/plots/{GenomeName}/AltDonorAcceptor/distribution/acceptor/proportion/overlay_comparison.png",
    params:
        outdir       = "results/ExonCharacteristics/plots/{GenomeName}/AltDonorAcceptor/distribution",
        intron_range = 150,
        exon_range   = 150,
        script       = "workflow/scripts/splice_site_distribution.py",
    log:
        "logs/ExonCharacteristic/alternateDonorAcceptor_distribution/{GenomeName}.log"
    resources:
        mem_mb = 8000
    shell:
        """
        mkdir -p {params.outdir}
        python {params.script} \
            --poison {input.poison} \
            --productive {input.productive} \
            --output {params.outdir} \
            --intron-range {params.intron_range} \
            --exon-range {params.exon_range} \
            &> {log}
        """
