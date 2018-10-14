rule merge_bisulfite_genome1:
    input:
        "../bisulfite/CpG_context_BC6.all.R1_val_1.fq_bismark_bt2_pe.sorted.genome1.txt.gz",
        "../bisulfite/CpG_context_BC7.all.R1_val_1.fq_bismark_bt2_pe.sorted.genome1.txt.gz",
        "../bisulfite/CpG_context_BC8.all.R1_val_1.fq_bismark_bt2_pe.sorted.genome1.txt.gz",
        "../bisulfite/CpG_context_BC9.all.R1_val_1.fq_bismark_bt2_pe.sorted.genome1.txt.gz"
    output:
        "../bisulfite/CpG_context.combined_replicates.genome1.summary.tsv"
    shell:
        "python summarize_bisulfite_methylation.py {output} {input}"

rule merge_bisulfite_genome2:
    input:
        "../bisulfite/CpG_context_BC6.all.R1_val_1.fq_bismark_bt2_pe.sorted.genome2.txt.gz",
        "../bisulfite/CpG_context_BC7.all.R1_val_1.fq_bismark_bt2_pe.sorted.genome2.txt.gz",
        "../bisulfite/CpG_context_BC8.all.R1_val_1.fq_bismark_bt2_pe.sorted.genome2.txt.gz",
        "../bisulfite/CpG_context_BC9.all.R1_val_1.fq_bismark_bt2_pe.sorted.genome2.txt.gz"
    output:
        "../bisulfite/CpG_context.combined_replicates.genome2.summary.tsv"
    shell:
        "python summarize_bisulfite_methylation.py {output} {input}"

rule merge_bisulfite_BC7:
    input:
        "../bisulfite/CpG_context_BC7.all.R1_val_1.fq_bismark_bt2_pe.txt.gz",
    output:
        "../bisulfite/CpG_context_BC7.all.R1_val_1.fq_bismark_bt2_pe.summary.tsv"
    shell:
        "python summarize_bisulfite_methylation.py {output} {input}"

rule ensembl_gtf_to_tsv:
    input:
        "../genome_data/ensembl_GRCm38.98.chr.gtf",
    output:
        "../genome_data/ensembl_GRCm38.98.chr.genes.tsv",
    script:
        "ensembl_gtf_to_tsv.R"

rule icr_plot_region_string:
    input:
        "../genome_data/ICR_plot_regions.tsv",
    output:
        "../genome_data/ICR_plot_regions_string.txt"
    shell:
        "python tsv_to_region.py {input} > {output}"
