rule b6_haplotype_analysis:
    input:
        "../RData/b6.minion/haplotype_df.RData",
        "../RData/b6.minion/summary_df.RData",
        "../RData/b6.minion/fit_reads_df.RData",
        notebook = "../notebooks/b6_haplotype_analysis.Rmd",
    output:
        "../notebooks/b6_haplotype_analysis.html",
        "../plots/b6.haplotype_score_combination.png",
    shell:
        "Rscript render_notebook.R {input.notebook}"

rule cast_haplotype_analysis:
    input:
        "../RData/cast.minion/haplotype_df.RData",
        "../RData/cast.minion/summary_df.RData",
        "../RData/cast.minion/fit_reads_df.RData",
        notebook = "../notebooks/cast_haplotype_analysis.Rmd",
    output:
        "../notebooks/cast_haplotype_analysis.html",
        "../plots/cast.haplotype_score_combination.png",
    shell:
        "Rscript render_notebook.R {input.notebook}"

rule b6xcast_haplotype_analysis:
    input:
        "../RData/b6xcast.minion/haplotype_df.RData",
        "../RData/b6xcast.minion/summary_df.RData",
        "../RData/b6xcast.minion/fit_reads_df.RData",
        notebook = "../notebooks/b6xcast_haplotype_analysis.Rmd",
    output:
        "../notebooks/b6xcast_haplotype_analysis.html",
        "../plots/b6xcast.haplotype_score_combination.pdf",
    shell:
        "Rscript render_notebook.R {input.notebook}"

rule castxb6_haplotype_analysis:
    input:
        "../RData/castxb6.promethion/haplotype_df.RData",
        "../RData/castxb6.promethion/summary_df.RData",
        "../RData/castxb6.promethion/fit_reads_df.RData",
        notebook = "../notebooks/castxb6_haplotype_analysis.Rmd",
    output:
        "../notebooks/castxb6_haplotype_analysis.html",
        "../plots/castxb6.haplotype_score_combination.png",
    shell:
        "Rscript render_notebook.R {input.notebook}"

rule haplotype_auc:
    input:
        "../RData/b6.minion/haplotype_df.RData",
        "../RData/cast.minion/haplotype_df.RData",
        notebook = "../notebooks/single_strain_auc.Rmd",
    output:
        "../notebooks/single_strain_auc.html",
    shell:
        "Rscript render_notebook.R {input.notebook}"

rule dmr_comparison:
    input:
        "../genome_data/Mus_musculus.GRCm38_90.chr.genes.tsv",
        "../tables/dss_dmrlist.csv",
        "../rna_seq/parent_biased_genes_10pctFDR.tsv",
        "../rna_seq/strain_biased_genes_5pctFDR.tsv",
        notebook = "../notebooks/compare_detected_dmrs.Rmd",
    output:
        "../notebooks/compare_detected_dmrs.html",
        "../plots/dss_distance/strain_linear.png"
    shell:
        "Rscript render_notebook.R {input.notebook}"

rule visualise_dmrs:
    input:
        "../RData/castxb6.promethion/fit_reads.RData",
        "../RData/castxb6.promethion/fit_reads_df.RData",
        "../tables/dss_dmrlist.csv",
        "../RData/b6xcast.minion/fit_reads.RData",
        "../RData/b6xcast.minion/fit_reads_df.RData",
        "../genome_data/GRCm38_90.cpg_coordinates.tsv",
        '../bisulfite/B6CastF1.combined_replicates.genome1.summary.tsv',
        '../bisulfite/B6CastF1.combined_replicates.genome2.summary.tsv',
        "../RData/knownGene.RData",
        "../RData/rna_seq_with_reverse.RData",
        "../genome_data/CGI_coordinates_mm10.masked.HMM.tsv",
        "../genome_data/primary_ICRs_mm10.tsv",
        "../genome_data/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf",
        notebook = "../notebooks/visualise_dmr.Rmd",
    output:
        "../notebooks/visualise_dmr.html",
        "../plots/dmr_examples/4933421O10Rik.pdf"
    shell:
        "Rscript render_notebook.R {input.notebook}"

rule differential_methylation_heatmaps:
    input:
        "../nanopore/b6xcast.minion.compare_haplotype_methylation.tsv",
        "../nanopore/castxb6.promethion.compare_haplotype_methylation.tsv",
        '../genome_data/ICR_plot_regions.tsv',
        notebook = "../notebooks/differential_methylation_heatmaps.Rmd",
    output:
        "../notebooks/differential_methylation_heatmaps.html",
        "../plots/primary_ICRs_heatmap.png",
    shell:
        "Rscript render_notebook.R {input.notebook}"

rule genome_level_methylation_summary:
    input:
        "../genome_data/Mus_musculus.GRCm38_90.chr.gtf",
        "../genome_data/CGI_coordinates_mm10.masked.HMM.tsv",
        "../RData/nanopolish_df.RData",
        "../RData/bisulfite_df.RData",
        notebook = "../notebooks/genome_level_methylation_summary.Rmd",
    output:
        "../notebooks/genome_level_methylation_summary.html",
        "../plots/gene_body_metaplot.pdf",
    shell:
        "Rscript render_notebook.R {input.notebook}"

rule nanopolish_methylation_validation:
    input:
        "../RData/nanopolish_df.RData",
        "../RData/bisulfite_df.RData",
        notebook = "../notebooks/nanopolish_methylation_validation.Rmd",
    output:
        "../notebooks/nanopolish_methylation_validation.html",
        "../plots/binned_site_concordance.pdf"
    shell:
        "Rscript render_notebook.R {input.notebook}"

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
               "../rna_seq/CastB6F1_5.hisat2.genome2.bam"),
        notebook = "../notebooks/rnaseq_analysis.Rmd",
    output:
        '../rna_seq/strain_biased_genes_5pctFDR.tsv',
        '../rna_seq/parent_biased_genes_10pctFDR.tsv',
        '../tables/imprinted_genes_Supp_file_1.tsv',
        "../notebooks/rnaseq_analysis.html",
        '../plots/allelic_bias.pdf',
    shell:
        "Rscript render_notebook.R {input.notebook}"

rule prepare_visualisation:
    input:
        "../genome_data/Mus_musculus.GRCm38_90.chr.gtf",
        "../rna_seq/all_runs_with_reverse_coverage.tsv",
    output:
        "../RData/knownGene.RData",
        "../RData/rna_seq_with_reverse.RData"
    script:
        "prepare_visualisation.R"
