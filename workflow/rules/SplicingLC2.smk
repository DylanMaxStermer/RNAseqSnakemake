wildcard_constraints:
    DonorsOrAcceptors = "Donors|Acceptors",

# Leafcutter2 Splicing Analysis 
# - Includes: Junction Classification

rule ExtractJuncs:
    input:
        bam = "results/Alignments/{sample}/Aligned.sortedByCoord.out.bam",
        index = "results/Alignments/{sample}/Aligned.sortedByCoord.out.bam.indexing_done",
    output:
        "results/SplicingAnalysis/juncfiles/{sample}.junc",
    params:
        strand = "0"
    conda:
        "../envs/regtools.yml"
    log:
        "logs/ExtractJuncs/{sample}.log"
    shell:
        """
        (regtools junctions extract -m 20 -s {params.strand} {input.bam} > {output}) &> {log}
        """

rule make_leafcutter_juncfile:
    input:
        ExpandAllSamplesInFormatStringFromGenomeNameWildcard("results/SplicingAnalysis/juncfiles/{sample}.junc"),
    output:
        "results/SplicingAnalysis/leafcutter/{GenomeName}/juncfilelist.txt"
    params:
        SamplesToRemove = ""
    run:
        import os
        if params.SamplesToRemove:
            SamplesToRemove = open(params.SamplesToRemove, 'r').read().split('\n')
        else:
            SamplesToRemove=[]
        with open(output[0], "w") as out:
            for filepath in input:
                samplename = os.path.basename(filepath).split(".junc")[0]
                if samplename not in  SamplesToRemove:
                    out.write(filepath + '\n')

rule leafcutter2:
    """
    TEST leafcutter2 for classification of juncts and for junction quantification
    Note: gene_type and transcript_type are the gencode names. gene_biotype and transcript_biotype are the other ....
    """
    input:
        juncs = ExpandAllSamplesInFormatStringFromGenomeNameWildcard("results/SplicingAnalysis/juncfiles/{sample}.junc"),
        juncfile_list = "results/SplicingAnalysis/leafcutter/{GenomeName}/juncfilelist.txt",

        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.gtf",
        fa = config['GenomesPrefix'] + "{GenomeName}/Reference.fa",
    output:
        outdir = directory("results/SplicingAnalysis/leafcutter2/{GenomeName}"),
        cluster_ratios = "results/SplicingAnalysis/leafcutter2/{GenomeName}/leafcutter2.cluster_ratios.const.gz",
        junction_counts = "results/SplicingAnalysis/leafcutter2/{GenomeName}/leafcutter2.junction_counts.const.gz",
        classifications = "results/SplicingAnalysis/leafcutter2/{GenomeName}/clustering/leafcutter2_junction_classifications.txt",
        clusters = "results/SplicingAnalysis/leafcutter2/{GenomeName}/clustering/leafcutter2_clusters",
    log:
        "logs/leafcutter2/{GenomeName}.log"
    conda:
        "../envs/leafcutter2.yml"
    resources:
        mem_mb = 20000
    shell:
        """
        mkdir -p {output.outdir}
        python workflow/scripts/leafcutter2/scripts/leafcutter2.py\
            -j {input.juncfile_list} \
            -r {output.outdir} \
            -o leafcutter2 \
            -A {input.gtf} --max_juncs 250 --log-level DEBUG \
            -G {input.fa} -D 0 -C &> {log}
        """

# Poison Exon, Donor, Acceptor & Skipped Identification Using Leafcutter2 Output

rule poisonExonID_leafcutter2:
    """
    This rule takes leafcutter annotation and count data and identifies potential poison exons.
    It also colors them in a BED file.

    The script creates two directories:
    - PE/ containing poison exon analysis files
    - ColoredBed/ containing colored BED files for visualization
    If I want to use output files of this for any other rule I will need to call those files directly.

    Will need to change the prefix part of the script
    """
    input:
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.gtf",
        cluster_ratios = "results/SplicingAnalysis/leafcutter2/{GenomeName}/leafcutter2.cluster_ratios.const.gz",
        classifications = "results/SplicingAnalysis/leafcutter2/{GenomeName}/clustering/leafcutter2_junction_classifications.txt",
    output:
        coding_noncoding_introns = "results/SplicingAnalysis/PoisonEvent/Exon/PE/{GenomeName}_coding_noncoding_introns_all.tsv",
    params:
        output_prefix = "results/SplicingAnalysis/PoisonEvent/Exon/{GenomeName}",
        max_distance = 500,  # Maximum distance between non-coding introns
    log:
        "logs/poisonExonID_lc2/{GenomeName}.log"
    resources:
        mem_mb = 28000
    shell:
        """
        python workflow/scripts/poisonExon_ID_leafcutter2.py \
            -g {input.gtf} \
            -c {input.cluster_ratios} \
            -l {input.classifications} \
            -o {params.output_prefix} \
            -d {params.max_distance} \
            &> {log}
        """



rule poisonDonorAcceptor_ID:
    """
    Identify poison donor and acceptor junctions from Leafcutter2 output.
    Detects alternative poison junctions where splice sites partially match CDS boundaries.
    Filters out junctions that are poison exon upstream/downstream introns.
    """
    input:
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.gtf",
        cluster_ratios = "results/SplicingAnalysis/leafcutter2/{GenomeName}/leafcutter2.cluster_ratios.const.gz",
        classifications = "results/SplicingAnalysis/leafcutter2/{GenomeName}/clustering/leafcutter2_junction_classifications.txt",
        poison_exon_introns = "results/SplicingAnalysis/PoisonEvent/Exon/PE/{GenomeName}_coding_noncoding_introns_all.tsv",
    output:
        donor_acceptor_bed = "results/SplicingAnalysis/PoisonEvent/DonorAcceptor/PoisonDonorAcceptor/{GenomeName}_poisonDonorAcceptor.bed",
    params:
        output_prefix = "results/SplicingAnalysis/PoisonEvent/DonorAcceptor/{GenomeName}",
    log:
        "logs/poisonDonorAcceptor_ID_lc2/{GenomeName}.log"
    resources:
        mem_mb = 28000
    shell:
        """
        python /project/yangili1/dylan_stermer/GitHubTesting/PoisonID/scripts/poisonDonorAcceptor_ID_leafcutter2.py \
            -g {input.gtf} \
            -l {input.classifications} \
            -c {input.cluster_ratios} \
            -p {input.poison_exon_introns} \
            -o {params.output_prefix} \
            &> {log}
        """

rule poisonSkipped_ID:
    """
    Identify poison skipped exon junctions from Leafcutter2 output.
    Detects non-coding junctions that span from one CDS boundary to another,
    representing skipping of normally productive exons.
    Filters out junctions that are poison exon upstream/downstream introns.
    """
    input:
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.gtf",
        cluster_ratios = "results/SplicingAnalysis/leafcutter2/{GenomeName}/leafcutter2.cluster_ratios.const.gz",
        classifications = "results/SplicingAnalysis/leafcutter2/{GenomeName}/clustering/leafcutter2_junction_classifications.txt",
        poison_exon_introns = "results/SplicingAnalysis/PoisonEvent/Exon/PE/{GenomeName}_coding_noncoding_introns_all.tsv",
    output:
        skipped_bed = "results/SplicingAnalysis/PoisonEvent/Skipped/PoisonSkipped/{GenomeName}_poisonSkipped.bed",
    params:
        output_prefix = "results/SplicingAnalysis/PoisonEvent/Skipped/{GenomeName}",
    log:
        "logs/poisonSkipped_ID_lc2/{GenomeName}.log"
    resources:
        mem_mb = 28000
    shell:
        """
        python /project/yangili1/dylan_stermer/GitHubTesting/PoisonID/scripts/poisonSkipped_ID_leafcutter2.py \
            -g {input.gtf} \
            -l {input.classifications} \
            -c {input.cluster_ratios} \
            -p {input.poison_exon_introns} \
            -o {params.output_prefix} \
            &> {log}
        """ 

rule productiveDonorAcceptor_ID:
    """
    Identify productive (coding) alternative donor and acceptor junctions
    from Leafcutter2 output.
    Detects alternative splice sites that preserve the reading frame relative
    to the canonical CDS boundaries.
    """
    input:
        gtf             = config['GenomesPrefix'] + "{GenomeName}/Reference.gtf",
        cluster_ratios  = "results/SplicingAnalysis/leafcutter2/{GenomeName}/leafcutter2.cluster_ratios.const.gz",
        classifications = "results/SplicingAnalysis/leafcutter2/{GenomeName}/clustering/leafcutter2_junction_classifications.txt",
    output:
        alt_donor_acceptor_bed = "results/SplicingAnalysis/productive/ProductiveDonorAcceptor/{GenomeName}_altDonorAcceptor.bed",
    params:
        output_prefix = "results/SplicingAnalysis/productive/{GenomeName}",
    log:
        "logs/ExonCharacteristic/productiveDonorAcceptor_ID/{GenomeName}.log"
    resources:
        mem_mb = 28000
    shell:
        """
        python workflow/scripts/productiveDonorAcceptor_ID_leafcutter2.py \
            -g {input.gtf} \
            -l {input.classifications} \
            -c {input.cluster_ratios} \
            -o {params.output_prefix} \
            &> {log}
        """


# Differential Splicing Analysis With Leafcutter2
# - The following 3 rules generate annotated junction files needed for tidy_leafcutter_differentialSplicing

rule ConcatJuncFilesAndKeepUniq:
    """
    Concatenate all junction files and keep unique junctions.
    Required for: tidy_leafcutter_differentialSplicing (generates intron_info input)
    """
    input:
        ExpandAllSamplesInFormatStringFromGenomeNameWildcard("results/SplicingAnalysis/juncfiles/{sample}.junc"),
    output:
        "results/SplicingAnalysis/ObservedJuncsAnnotations/{GenomeName}.uniq.junc.gz"
    log:
        "logs/ConcatJuncFilesAndKeepUniq/{GenomeName}.log"
    resources:
        mem_mb = GetMemForSuccessiveAttempts(24000, 48000)
    shell:
        """
        (awk -v OFS='\\t' '{{ split($11, blockSizes, ","); JuncStart=$2+blockSizes[1]; JuncEnd=$3-blockSizes[2]; print $0, JuncStart, JuncEnd }}' {input} | sort -k1,1 -k6,6 -k13,13n -k14,14n | python workflow/scripts/AggregateSortedCattedJuncBeds.py | bedtools sort -i - | bgzip -c /dev/stdin  > {output}) &> {log}
        """

rule AnnotateConcatedUniqJuncFile_basic:
    """
    Annotate unique junctions using regtools.
    Required for: tidy_leafcutter_differentialSplicing (generates intron_info input)
    Note: regtools end coordinate is off by one. The right coordinate should be adjusted down 1 for proper viewing of junc in IGV.
    This command breaks when gtf has trailing spaces, like in some NCBI gtfs: https://github.com/griffithlab/regtools/issues/92
    """
    input:
        junc = "results/SplicingAnalysis/ObservedJuncsAnnotations/{GenomeName}.uniq.junc.gz",
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.basic.gtf",
        fa = config['GenomesPrefix'] + "{GenomeName}/Reference.fa"
    output:
        "results/SplicingAnalysis/ObservedJuncsAnnotations/{GenomeName}.uniq.annotated.tsv.gz"
    log:
        "logs/AnnotateConcatedUniqJuncFile_hg38Basic.{GenomeName}.log"
    conda:
        "../envs/regtools.yml"
    shell:
        """
        ( zcat {input.junc} \
            | grep -v 'random' \
            | grep -v 'chrUn' \
            | regtools junctions annotate - {input.fa} {input.gtf} \
            | gzip -c > {output} ) &> {log}
        """

rule Add_splice_site_scores_to_regtools_annotate:
    """
    Add splice site scores to annotated junctions.
    Required for: tidy_leafcutter_differentialSplicing (generates intron_info input)
    """
    input:
        annotated_junctions="results/SplicingAnalysis/ObservedJuncsAnnotations/{GenomeName}.uniq.annotated.tsv.gz",
        reference_fasta=config['GenomesPrefix'] + "{GenomeName}/Reference.fa",
        fai = config['GenomesPrefix'] + "{GenomeName}/Reference.fa.fai",
        AnnotatedIntronsWithSS = config['GenomesPrefix'] + "{GenomeName}/Reference.Introns.bed.gz"
    output:
        "results/SplicingAnalysis/ObservedJuncsAnnotations/{GenomeName}.uniq.annotated.with_ss_scores.tsv.gz"
    conda:
        "../scripts/leafcutter2/scripts/Reformat_gtf.conda_env.yml"
    log:
        "logs/Add_splice_site_scores_to_regtools_annotate.{GenomeName}.log"
    resources:
        mem_mb = GetMemForSuccessiveAttempts(24000, 48000)
    shell:
        """
        python workflow/scripts/Add_SS_To_RegtoolsAnnotate.py \
            --input {input.annotated_junctions} \
            --reference {input.reference_fasta} \
            --introns {input.AnnotatedIntronsWithSS} \
            --output {output} &> {log}
        """

rule rename_lc2_rows:
    """
    This rule renames the rows of the leafcutter2 junction counts file to match the leafcutter1 format.
    This is needed because the differential splicing script expects the leafcutter1 format.
    """
    input:
        junction_counts = "results/SplicingAnalysis/leafcutter2/{GenomeName}/leafcutter2.junction_counts.const.gz",
    output:
        renamed_counts = "results/SplicingAnalysis/leafcutter2/{GenomeName}/leafcutter2_renamed.junction_counts.gz",
    conda:
        "../envs/leafcutter2.yml"
    log:
        "logs/rename_lc2_rows/{GenomeName}.log"
    shell:
        """
        zcat {input.junction_counts} | \
        awk -F'\\t' 'BEGIN {{OFS="\\t"}} NR==1 {{print; next}} {{sub(/:[^:]+$/, "", $1); print}}' | \
        gzip > {output.renamed_counts} 2> {log}
        """

rule leafcutter_ds_contrasts:
    """
    Differential splicing analysis using leafcutter_ds.R with leafcutter2 junction counts.
    Uses the renamed LC2 junction counts that match the LC1 format.
    """
    input:
        groupfile = lambda wildcards: os.path.abspath(config['contrast_group_files_prefix'] + "{contrast}.txt"),
        numers = lambda wildcards: get_filled_path_from_contrast(wildcards, "results/SplicingAnalysis/leafcutter2/{GenomeName}/leafcutter2_renamed.junction_counts.gz"),
    output:
        outputdir = directory("results/SplicingAnalysis/differential_splicing/{contrast}/"),
        effect_sizes = "results/SplicingAnalysis/differential_splicing/{contrast}/leaf_effect_sizes.txt",
        pvalues = "results/SplicingAnalysis/differential_splicing/{contrast}/leaf_cluster_significance.txt",
    threads: 4
    wildcard_constraints:
        treatment = "|".join(contrasts)
    resources:
        ntasks = 5,
        mem_mb = 24000
    params:
        ExtraParams = "",
        MinGroupSize = min_samples_per_group_for_contrast
    log:
        "logs/leafcutter_ds/{contrast}.log"
    shell:
        """
        export PATH=/project/yangili1/dylan_stermer/miniconda3/conda-envs/fkoompa/bin:$PATH
        export LD_LIBRARY_PATH=/project/yangili1/dylan_stermer/miniconda3/conda-envs/fkoompa/lib
        export R_LIBS_USER=""
        export R_LIBS=""

        mkdir -p {output.outputdir}
        Rscript workflow/scripts/leafcutter/scripts/leafcutter_ds.R -p {threads} -o {output.outputdir}/leaf {params.ExtraParams} -i {params.MinGroupSize} -g {params.MinGroupSize} {input.numers} {input.groupfile} &> {log}
        """

rule tidy_leafcutter_differentialSplicing:
    """
    Join effect sizes, p values, and intron info into a single table.
    Uses leafcutter2 classifications for junction annotation.
    """
    input:
        intron_info = lambda wildcards: get_filled_path_from_contrast(wildcards, "results/SplicingAnalysis/ObservedJuncsAnnotations/{GenomeName}.uniq.annotated.with_ss_scores.tsv.gz"),
        Classifications = lambda wildcards: get_filled_path_from_contrast(wildcards, "results/SplicingAnalysis/leafcutter2/{GenomeName}/clustering/leafcutter2_junction_classifications.txt"),
        effect_sizes = "results/SplicingAnalysis/differential_splicing/{contrast}/leaf_effect_sizes.txt",
        pvalues = "results/SplicingAnalysis/differential_splicing/{contrast}/leaf_cluster_significance.txt",
    output:
        "results/SplicingAnalysis/differential_splicing_tidy/{contrast}/Results.tsv.gz",
    log:
        "logs/tidy_leafcutter_differentialSplicing/{contrast}.log"
    resources:
        mem_mb = 48000
    conda:
        "../envs/r_2.yml"
    shell:
        """
        Rscript workflow/scripts/leafcutter_tidyDifferentialSplicingResults.R {input.intron_info} {input.Classifications} {input.effect_sizes} {input.pvalues} {output} &> {log}
        """