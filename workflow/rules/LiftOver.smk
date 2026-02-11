# ============================================================================
# Reciprocal LiftOver
# ============================================================================
# Lifts exons_std.tsv coordinates from hg38 to each target assembly using
# reciprocal liftover (forward lift then map-back filter for 1:1 orthologs).
#
# Config keys (config.yaml):
#   liftover.source_assembly               — source assembly name (e.g. "hg38")
#   liftover.targets.{name}.forward_chain  — chain file hg38 -> target
#   liftover.targets.{name}.reverse_chain  — chain file target -> hg38
#
# Wildcards:
#   {GenomeName}  — genome annotation build (from STAR_Genome_List)
#   {exon_type}   — exon category subdirectory under StandardFormat/
#   {target}      — target assembly name (key in liftover.targets config)
# ============================================================================

rule reciprocal_liftover:
    """
    Run reciprocal liftover for a single exon_type x target assembly combination.
    Outputs one_way/, reciprocal/, and ortholog_mapping/ subdirectories.
    """
    input:
        exons = "results/ExonCharacteristics/{GenomeName}/StandardFormat/{exon_type}/exons_std.tsv",
    output:
        reciprocal_lifted   = "results/liftOver/{GenomeName}/{exon_type}/{target}/reciprocal/hg38_to_{target}_reciprocal_lifted.tsv",
        reciprocal_unmapped = "results/liftOver/{GenomeName}/{exon_type}/{target}/reciprocal/hg38_to_{target}_reciprocal_unmapped.tsv",
        one_way_lifted      = "results/liftOver/{GenomeName}/{exon_type}/{target}/one_way/hg38_to_{target}_lifted.tsv",
        one_way_unmapped    = "results/liftOver/{GenomeName}/{exon_type}/{target}/one_way/hg38_to_{target}_unmapped.tsv",
        orthologs_1to1      = "results/liftOver/{GenomeName}/{exon_type}/{target}/ortholog_mapping/hg38_to_{target}_1to1_orthologs.tsv",
        orthologs_1tomany   = "results/liftOver/{GenomeName}/{exon_type}/{target}/ortholog_mapping/hg38_to_{target}_1tomany_orthologs.tsv",
        quality_metrics     = "results/liftOver/{GenomeName}/{exon_type}/{target}/ortholog_mapping/hg38_to_{target}_quality_metrics.tsv",
    params:
        script        = "workflow/scripts/04_reciprocal_liftover_workflow.py",
        outdir        = "results/liftOver/{GenomeName}/{exon_type}/{target}",
        source        = config['liftover']['source_assembly'],
        forward_chain = lambda wc: config['liftover']['targets'][wc.target]['forward_chain'],
        reverse_chain = lambda wc: config['liftover']['targets'][wc.target]['reverse_chain'],
    log:
        "logs/LiftOver/{GenomeName}/{exon_type}/{target}.log",
    resources:
        mem_mb = 8000,
    shell:
        """
        mkdir -p {params.outdir}
        python {params.script} \
            -i {input.exons} \
            -f {params.forward_chain} \
            -r {params.reverse_chain} \
            -o {params.outdir} \
            -s {params.source} \
            -t {wildcards.target} \
            &> {log}
        """
