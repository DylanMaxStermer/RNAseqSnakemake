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
        cluster_ratios = "results/SplicingAnalysis/leafcutter2/{GenomeName}/leafcutter2.cluster_ratios.gz",
        junction_counts = "results/SplicingAnalysis/leafcutter2/{GenomeName}/leafcutter2.junction_counts.gz",
        classifications = "results/SplicingAnalysis/leafcutter2/{GenomeName}/clustering/leafcutter2_junction_classifications.txt",
    log:
        "logs/leafcutter2/{GenomeName}.log"
    conda:
        "/project/yangili1/dylan_stermer/GitHubTesting/RNAseq_Snakemake/workflow/scripts/leafcutter2_carlos2/scripts/add_on_scripts/lc2_reformat.yml"
    resources:
        mem_mb = 20000
    shell:
        """
        mkdir -p {output.outdir}
        python /project/yangili1/dylan_stermer/GitHubTesting/RNAseq_Snakemake/workflow/scripts/leafcutter2/scripts/leafcutter2.py\
            -j {input.juncfile_list} \
            -r {output.outdir} \
            -o leafcutter2 \
            -A {input.gtf} --max_juncs 250 --log-level DEBUG \
            -G {input.fa} &> {log}
        """

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
        cluster_ratios = "results/SplicingAnalysis/leafcutter2/{GenomeName}/leafcutter2.cluster_ratios.gz",
        classifications = "results/SplicingAnalysis/leafcutter2/{GenomeName}/clustering/leafcutter2_junction_classifications.txt",
    output:
        pe_dir = directory("results/SplicingAnalysis/poisonExon/{GenomeName}/lc2/PE"),
        colored_dir = directory("results/SplicingAnalysis/poisonExon/{GenomeName}/lc2/ColoredBed"),
    params:
        output_prefix = "results/SplicingAnalysis/poisonExon/{GenomeName}/lc2/{GenomeName}",
        max_distance = 500,  # Maximum distance between non-coding introns
    log:
        "logs/poisonExonID_lc2_test/{GenomeName}.log"
    resources:
        mem_mb = 48000
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
