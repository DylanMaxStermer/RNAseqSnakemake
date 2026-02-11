rule MakeBigwigs_NormalizedToGenomewideCoverage:
    """
    Scale bigwig to base coverage per billion chromosomal reads
    """
    input:
        fai = lambda wildcards: config['GenomesPrefix'] + samples.loc[samples['sample']==wildcards.sample]['STARGenomeName'].tolist()[0] + "/Reference.fa.fai",
        bam = "results/Alignments/{sample}/Aligned.sortedByCoord.out.bam",
        bai = "results/Alignments/{sample}/Aligned.sortedByCoord.out.bam.indexing_done",
        NormFactorsFile = "results/QC/ReadCountsPerSamples.tsv"
    params:
        GenomeCovArgs="-split",
        bw_minus = "bw_minus=",
        MKTEMP_ARGS = "-p " + config['scratch'],
        SORT_ARGS="-T " + config['scratch'],
        Region = "",
    shadow: "shallow"
    output:
        bw = "results/bigwigs/unstranded/{sample}.bw",
        bw_minus = []
    log:
        "logs/MakeBigwigs_unstranded/{sample}.log"
    resources:
        mem_mb = GetMemForSuccessiveAttempts(42000, 52000)
    shell:
        """
        ScaleFactor=$(bc <<< "scale=3;1000000000/$(grep '{wildcards.sample}' {input.NormFactorsFile} | awk 'NR==1 {{print $2}}')")
        workflow/scripts/BamToBigwig.sh {input.fai} {input.bam} {output.bw}  GENOMECOV_ARGS="{params.GenomeCovArgs} -scale ${{ScaleFactor}}" REGION='{params.Region}' MKTEMP_ARGS="{params.MKTEMP_ARGS}" SORT_ARGS="{params.SORT_ARGS}" {params.bw_minus}"{output.bw_minus}" &> {log}
        """

rule MakeBigwigs_Grouped_Mean:
    """
    Create grouped bigwigs by averaging normalized individual sample bigwigs.
    Groups defined by bigWigGroup column in samples.tsv.
    """
    input:
        bigwigs = get_samples_for_bigwig_group,
        fai = lambda wildcards: config['GenomesPrefix'] + get_genome_for_bigwig_group(wildcards) + "/Reference.fa.fai",
    output:
        bw = "results/bigwigs/grouped/{group}.bw",
    log:
        "logs/MakeBigwigs_grouped/{group}.log"
    conda:
        "../envs/wiggletools.yaml"
    shadow: "shallow"
    params:
        MKTEMP_ARGS = "-p " + config['scratch'],
        SORT_ARGS = "-T " + config['scratch'],
    resources:
        mem_mb = GetMemForSuccessiveAttempts(16000, 32000)
    shell:
        """
        wiggletools write_bg - mean {input.bigwigs} | \
        LC_COLLATE=C sort -k1,1 -k2,2n {params.SORT_ARGS} > temp_mean.bedGraph

        bedGraphToBigWig temp_mean.bedGraph {input.fai} {output.bw} &> {log}
        """