rule b6_haplotype_analysis:
    input:
        "../RData/albacore_1.2.2.b6/haplotype_df.RData",
        "../RData/albacore_1.2.2.b6/summary_df.RData",
        "../RData/albacore_1.2.2.b6/fit_reads_df.RData",
    output:
        "../notebooks/b6_haplotype_analysis.html",
        "../plots/b6.haplotype_score_combination.png",
    shell:
        "Rscript render_notebook.R ../notebooks/b6_haplotype_analysis.Rmd"

rule cast_haplotype_analysis:
    input:
        "../RData/albacore_1.2.2.cast/haplotype_df.RData",
        "../RData/albacore_1.2.2.cast/summary_df.RData",
        "../RData/albacore_1.2.2.cast/fit_reads_df.RData",
    output:
        "../notebooks/cast_haplotype_analysis.html",
        "../plots/cast.haplotype_score_combination.png",
    shell:
        "Rscript render_notebook.R ../notebooks/cast_haplotype_analysis.Rmd"

rule b6xcast_haplotype_analysis:
    input:
        "../RData/albacore_1.2.2.b6xcast/haplotype_df.RData",
        "../RData/albacore_1.2.2.b6xcast/summary_df.RData",
        "../RData/albacore_1.2.2.b6xcast/fit_reads_df.RData",
    output:
        "../notebooks/b6xcast_haplotype_analysis.html",
        "../plots/b6xcast.haplotype_score_combination.pdf",
    shell:
        "Rscript render_notebook.R ../notebooks/b6xcast_haplotype_analysis.Rmd"

rule castxb6_haplotype_analysis:
    input:
        "../RData/albacore_2.2.7.castxb6.promethion/haplotype_df.RData",
        "../RData/albacore_2.2.7.castxb6.promethion/summary_df.RData",
        "../RData/albacore_2.2.7.castxb6.promethion/fit_reads_df.RData",
    output:
        "../notebooks/castxb6_haplotype_analysis.html",
        "../plots/castxb6.haplotype_score_combination.png",
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
        "../genome_data/Mus_musculus.GRCm38_90.chr.genes.tsv",
        "../tables/dss_dmrlist.csv",
        "../rna_seq/parent_biased_genes_10pctFDR.tsv",
        "../rna_seq/strain_biased_genes_5pctFDR.tsv",
    output:
        "../notebooks/compare_detected_dmrs.html",
        "../plots/dss_distance/strain_linear.png"
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
        '../bisulfite/B6CastF1.combined_replicates.genome1.summary.tsv',
        '../bisulfite/B6CastF1.combined_replicates.genome2.summary.tsv',
        "../RData/knownGene.RData",
        "../RData/rna_seq_with_reverse.RData",
        "../genome_data/CGI_coordinates_mm10.masked.HMM.tsv",
        "../genome_data/primary_ICRs_mm10.tsv",
        "../genome_data/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf",
    output:
        "../notebooks/visualise_dmr.html",
        "../plots/dmr_examples/4933421O10Rik.pdf"
    shell:
        "Rscript render_notebook.R ../notebooks/visualise_dmr.Rmd"

rule differential_methylation_heatmaps:
    input:
        "../nanopore/albacore_1.2.2.b6xcast.compare_haplotype_methylation.tsv",
        "../nanopore/albacore_2.2.7.castxb6.promethion.compare_haplotype_methylation.tsv",
        '../genome_data/ICR_plot_regions.tsv',
    output:
        "../notebooks/differential_methylation_heatmaps.html",
        "../plots/primary_ICRs_heatmap.png",
    shell:
        "Rscript render_notebook.R ../notebooks/differential_methylation_heatmaps.Rmd"

rule genome_level_methylation_summary:
    input:
        "../genome_data/Mus_musculus.GRCm38_90.chr.gtf",
        "../genome_data/CGI_coordinates_mm10.masked.HMM.tsv",
        "../RData/nanopolish_df.RData",
        "../RData/bisulfite_df.RData",
    output:
        "../notebooks/genome_level_methylation_summary.html",
        "../plots/gene_body_metaplot.pdf",
    shell:
        "Rscript render_notebook.R ../notebooks/genome_level_methylation_summary.Rmd"

rule nanopolish_methylation_validation:
    input:
        "../RData/nanopolish_df.RData",
        "../RData/bisulfite_df.RData",
    output:
        "../notebooks/nanopolish_methylation_validation.html",
        "../plots/binned_site_concordance.pdf"
    shell:
        "Rscript render_notebook.R ../notebooks/nanopolish_methylation_validation.Rmd"

rule rnaseq_analysis_notebook:
    input:
        "../genome_data/Mus_musculus.GRCm38_90.chr.gtf",
        '../genome_data/known_imprinted_genes.RData',
        bam = ("../rna_seq/B6CastF1_2.hisat2.genome1.bam",
               "../rna_seq/B6CastF1_2.hisat2.genome2.bam",
               "../rna_seq/B6CastF1_3.hisat2.genome1.bam",
               "../rna_seq/B6CastF1_3.hisat2.genome2.bam",
               "../rna_seq/B6CastF1_4.hisat2.genome1.bam",
               "../rna_seq/B6CastF1_4.hisat2.genome2.bam",
               "../rna_seq/B6CastF1_5.hisat2.genome1.bam",
               "../rna_seq/B6CastF1_5.hisat2.genome2.bam",
               "../rna_seq/CastB6F1_2.hisat2.genome1.bam",
               "../rna_seq/CastB6F1_2.hisat2.genome2.bam",
               "../rna_seq/CastB6F1_3.hisat2.genome1.bam",
               "../rna_seq/CastB6F1_3.hisat2.genome2.bam",
               "../rna_seq/CastB6F1_4.hisat2.genome1.bam",
               "../rna_seq/CastB6F1_4.hisat2.genome2.bam",
               "../rna_seq/CastB6F1_5.hisat2.genome1.bam",
               "../rna_seq/CastB6F1_5.hisat2.genome2.bam")
    output:
        '../rna_seq/strain_biased_genes_5pctFDR.tsv',
        '../rna_seq/parent_biased_genes_10pctFDR.tsv',
        '../tables/imprinted_genes_Supp_file_1.tsv',
        "../notebooks/rnaseq_analysis.html",
        '../plots/allelic_bias.pdf',
    shell:
        "Rscript render_notebook.R ../notebooks/rnaseq_analysis.Rmd"

rule prepare_visualisation:
    input:
        "../genome_data/Mus_musculus.GRCm38_90.chr.gtf",
        "../rna_seq/all_runs_with_reverse_coverage.tsv"
    output:
        "../RData/knownGene.RData",
        "../RData/rna_seq_with_reverse.RData"
    script:
        "prepare_visualisation.R"
