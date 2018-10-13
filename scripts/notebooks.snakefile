rule b6_haplotype_analysis:
    input:
        "../RData/albacore_1.2.2.b6/haplotype_df.RData",
        "../RData/albacore_1.2.2.b6/summary_df.RData",
        "../RData/albacore_1.2.2.b6/fit_reads_df.RData",
        "../plots/"
    output:
        "../notebooks/b6_haplotype_analysis.html",
    shell:
        "Rscript render_notebook.R ../notebooks/b6_haplotype_analysis.Rmd"

rule cast_haplotype_analysis:
    input:
        "../RData/albacore_1.2.2.cast/haplotype_df.RData",
        "../RData/albacore_1.2.2.cast/summary_df.RData",
        "../RData/albacore_1.2.2.cast/fit_reads_df.RData",
        "../plots/"
    output:
        "../notebooks/cast_haplotype_analysis.html",
    shell:
        "Rscript render_notebook.R ../notebooks/cast_haplotype_analysis.Rmd"

rule forward_cross_haplotype_analysis:
    input:
        "../RData/albacore_1.2.2.b6xcast/haplotype_df.RData",
        "../RData/albacore_1.2.2.b6xcast/summary_df.RData",
        "../RData/albacore_1.2.2.b6xcast/fit_reads_df.RData",
        "../plots/"
    output:
        "../notebooks/b6xcast_haplotype_analysis.html",
    shell:
        "Rscript render_notebook.R ../notebooks/b6xcast_haplotype_analysis.Rmd"

rule reverse_cross_haplotype_analysis:
    input:
        "../RData/albacore_2.2.7.castxb6.promethion/haplotype_df.RData",
        "../RData/albacore_2.2.7.castxb6.promethion/summary_df.RData",
        "../RData/albacore_2.2.7.castxb6.promethion/fit_reads_df.RData",
        "../plots/"
    output:
        "../notebooks/castxb6_haplotype_analysis.html",
    shell:
        "Rscript render_notebook.R ../notebooks/castxb6_haplotype_analysis.Rmd"

rule haplotype_auc:
    input:
        "../RData/albacore_1.2.2.b6/haplotype_df.RData",
        "../RData/albacore_1.2.2.cast/haplotype_df.RData"
    output:
        "../notebooks/single_strain_auc.html",
    shell:
        "Rscript render_notebook.R ../notebooks/single_strain_auc.Rmd"

rule dmr_comparison:
    input:
        "../genome_data/ensembl_GRCm38.98.chr.genes.tsv",
        "../tables/dss_dmrlist.csv",
        "../rna_seq/parent_biased_genes_10pctFDR.tsv",
        "../rna_seq/strain_biased_genes_5pctFDR.tsv",
        "../plots/",
        "../tables/",
    output:
        "../notebooks/compare_detected_dmrs.html"
    shell:
        "Rscript render_notebook.R ../notebooks/compare_detected_dmrs.Rmd"

rule visualise_dmrs:
    input:
        "../RData/albacore_2.2.7.castxb6.promethion/fit_reads.RData",
        "../RData/albacore_2.2.7.castxb6.promethion/fit_reads_df.RData",
        "../tables/dss_dmrlist.csv",
        "../RData/albacore_1.2.2.b6xcast/fit_reads.RData",
        "../RData/albacore_1.2.2.b6xcast/fit_reads_df.RData",
        "../genome_data/CpG_coordinates_mm10.RData",
        '../bisulfite/CpG_context.combined_replicates.genome1.summary.tsv',
        '../bisulfite/CpG_context.combined_replicates.genome2.summary.tsv',
        "../RData/knownGene.RData",
        "../RData/rna_seq_with_reverse2.RData",
        "../genome_data/CGI_coordinates_mm10.masked.HMM.tsv",
        "../genome_data/primary_ICRs_mm10.tsv",
        "../genome_data/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf",
        "../plots/",
    output:
        "../notebooks/visualise_dmr.html"
    shell:
        "Rscript render_notebook.R ../notebooks/visualise_dmr.Rmd"

rule differential_methylation_heatmaps:
    input:
        "../nanopore/albacore_1.2.2.b6xcast.compare_haplotype_methylation.tsv",
        "../nanopore/albacore_2.2.7.castxb6.promethion.compare_haplotype_methylation.tsv",
        '../genome_data/ICR_plot_regions.tsv',
        "../plots/",
    output:
        "../notebooks/differential_methylation_heatmaps.html"
    shell:
        "Rscript render_notebook.R ../notebooks/differential_methylation_heatmaps.Rmd"

rule genome_level_methylation_summary:
    input:
        "../genome_data/ensembl_GRCm38.98.chr.gtf",
        "../genome_data/CGI_coordinates_mm10.masked.HMM.tsv",
        "../RData/nanopolish_df.RData",
        "../RData/bisulfite_df.RData",
        "../plots/",
    output:
        "../notebooks/genome_level_methylation_summary.html"
    shell:
        "Rscript render_notebook.R ../notebooks/genome_level_methylation_summary.Rmd"

rule nanopolish_methylation_validation:
    input:
        "../RData/nanopolish_df.RData",
        "../RData/bisulfite_df.RData",
        "../plots/",
    output:
        "../notebooks/nanopolish_methylation_validation.html"
    shell:
        "Rscript render_notebook.R ../notebooks/nanopolish_methylation_validation.Rmd"

rule rnaseq_analysis_notebook:
    input:
        "../rna_seq/repeat_analysis_Sep18_imprint.RData",
        "../rna_seq/imprinted_genes_Supp_file_1.tsv",
        "../rna_seq/repeat_analysis_Sep18_strain.RData",
        "../plots/",
    output:
        "../notebooks/rnaseq_analysis.html"
    shell:
        "Rscript render_notebook.R ../notebooks/rnaseq_analysis.Rmd"

rule notebook_setup:
    input:
        "../bisulfite/CpG_context_TB1_PlacentaE14.5_WT_R1_trimmed.fq_bismark_bt2_pe.summary.tsv",
        "../nanopore/albacore_1.2.2.b6xcast.methylation.summary.tsv",
        "../genome_data/ensembl_GRCm38.98.chr.gtf",
        "../genome_data/ensembl_GRCm38.98.chr.genes.tsv",
        "../RData/paired_DSS.RData",
    output:
        "../tables/dss_dmrlist.csv",
        "../RData/knownGene.RData",
        "../RData/nanopolish_df.RData",
        "../RData/bisulfite_df.RData",
        "../plots/",
    script:
        "notebook_setup.R"
