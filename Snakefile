workdir:
    "scripts"
include:
    "scripts/allelic_methylation.snakefile"
include:
    "scripts/fvb_resolution.snakefile"
include:
    "scripts/nanopolish.snakefile"
include:
    "scripts/notebooks.snakefile"
include:
    "scripts/misc.snakefile"
include:
    "scripts/rna_seq.snakefile"
include:
    "scripts/rrbs.snakefile"
conda:
    "environment.yaml"

"""Inputs:
        "../genome_data/imprinted_genes.tsv",
        "../genome_data/primary_ICRs_mm10.tsv",
        "../genome_data/ICR_plot_regions.tsv",
"""

"""TODO:
Maybe do: quality analysis
"""

rule all:
    input:
        "../notebooks/b6_haplotype_analysis.html",
        "../notebooks/cast_haplotype_analysis.html",
        "../notebooks/castxb6_haplotype_analysis.html",
        "../notebooks/b6xcast_haplotype_analysis.html",
        "../notebooks/single_strain_auc.html",
        "../notebooks/differential_methylation_heatmaps.html",
        "../notebooks/visualise_dmr.html",
        "../notebooks/compare_detected_dmrs.html",
        "../notebooks/genome_level_methylation_summary.html",
        "../notebooks/nanopolish_methylation_validation.html",
        "../notebooks/rnaseq_analysis.html",

rule intermediate_download:
    input:
        "../nanopore/b6xcast.minion.intermediate_download",
        "../nanopore/b6.minion.intermediate_download",
        "../nanopore/cast.minion.intermediate_download",
        "../nanopore/castxb6.promethion.intermediate_download",

rule haplotype_analysis:
    input:
        "../notebooks/b6_haplotype_analysis.html",
        "../notebooks/cast_haplotype_analysis.html",
        "../notebooks/castxb6_haplotype_analysis.html",
        "../notebooks/b6xcast_haplotype_analysis.html",
        "../notebooks/single_strain_auc.html",

rule dmr_analysis:
    input:
        "../notebooks/differential_methylation_heatmaps.html",
        "../notebooks/visualise_dmr.html",
        "../notebooks/compare_detected_dmrs.html"

rule methylation_analysis:
    input:
        "../notebooks/genome_level_methylation_summary.html",
        "../notebooks/nanopolish_methylation_validation.html",

rule rnaseq_analysis:
    input:
        "../notebooks/rnaseq_analysis.html"
