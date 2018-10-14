rule suppdb_bam:
    input:
        bam = "../nanopore/{sample}.sorted.bam",
        bai = "../nanopore/{sample}.sorted.bam.bai",
    output:
        suppdb = "../nanopore/{sample}.sorted.bam.suppdb",
        summary = "../nanopore/{sample}.summary.tsv"
    shell:
        "python build_supplementary_db.py {input.bam} --summary {output.summary}"

rule split_methylation_by_alignment:
    input:
        bam = "../nanopore/{sample}.sorted.bam",
        meth = "../nanopore/{sample}.methylation.tsv",
        suppdb = "../nanopore/{sample}.sorted.bam.suppdb",
    output:
        "../nanopore/{sample}.methylation.split_supplementary.tsv"
    shell:
        "python split_methylation_by_alignment.py {input.bam} {input.meth} > {output}"

rule calculate_methylation_frequency:
    input:
        "../nanopore/{sample}.methylation.tsv"
    output:
        "../nanopore/{sample}.methylation.summary.tsv"
    shell:
        "python calculate_methylation_frequency.py -i {input} -p > {output}"

rule call_variant_proportion:
    input:
        bam = "../nanopore/{sample}.sorted.bam",
        bam_index = "../nanopore/{sample}.sorted.bam.bai",
        suppdb = "../nanopore/{sample}.sorted.bam.suppdb",
        phased_bam = "../nanopore/{sample}.phased_sorted.bam",
        phased_bam_index = "../nanopore/{sample}.phased_sorted.bam.bai",
        vcf = "../genome_data/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf",
    output:
        "../nanopore/{sample}.phased.tsv"
    shell:
        "python call_variant_proportion.py -b {input.bam} -p {input.phased_bam} -v {input.vcf} -o {output}"

rule sort_by_read:
    input:
        meth_split = "../nanopore/{sample}.methylation.split_supplementary.tsv",
        meth = "../nanopore/{sample}.methylation.tsv"
    output:
        "../nanopore/{sample}.methylation.sorted.by_read.tsv"
    shell:
        "mkdir -p {output}.tmp && tail -n +2 {input.meth_split} | "
        "sort -k4,4 -k1,1 -k2n,2n -T {output}.tmp | "
        "cat <(head -n 1 {input.meth}) - > {output} && "
        "rm -rf {output}.tmp"

rule sort_by_site:
    input:
        meth_split = "../nanopore/{sample}.methylation.split_supplementary.tsv",
        meth = "../nanopore/{sample}.methylation.tsv"
    output:
        "../nanopore/{sample}.methylation.sorted.by_site.tsv"
    shell:
        "mkdir -p {output}.tmp && tail -n +2 {input.meth_split} | "
        "sort -k1,1 -k2n,2n -k4,4 -T {output}.tmp | "
        "cat <(head -n 1 {input.meth}) - > {output} && "
        "rm -rf {output}.tmp"

rule fit_reads:
    input:
        "../nanopore/{sample}.summary.tsv",
        "../nanopore/{sample}.phased.tsv",
        "../nanopore/{sample}.methylation.sorted.by_read.tsv"
    output:
        "../RData/{sample}/summary_df.RData",
        "../RData/{sample}/haplotype_df.RData",
        "../RData/{sample}/fit_reads.RData",
        "../RData/{sample}/fit_reads_df.RData",
    params:
        sample = lambda wildcards, output: wildcards.sample
    shell:
        "Rscript fit_reads.R ../nanopore/{params.sample} ../RData/{params.sample}"

rule split_methylation_by_haplotype:
    input:
        meth = "../nanopore/{sample}.methylation.sorted.by_site.tsv",
        phase = "../nanopore/{sample}.phased.tsv"
    output:
        "../nanopore/{sample}.methylation.sorted.by_site.tsv.ref.tsv",
        "../nanopore/{sample}.methylation.sorted.by_site.tsv.alt.tsv"
    shell:
        "python split_methylation_by_haplotype.py -m {input.meth} -p {input.phase}"

rule calculate_allele_meth_freq:
    input:
        "../nanopore/{sample}.methylation.sorted.by_site.tsv.{allele}.tsv",
    output:
        "../nanopore/{sample}.methylation.{allele}_summary.tsv"
    shell:
        "python calculate_methylation_frequency.py -i {input} -p > {output}"

rule paired_dss:
    input:
        "../nanopore/albacore_1.2.2.b6xcast.methylation.ref_summary.tsv",
        "../nanopore/albacore_1.2.2.b6xcast.methylation.alt_summary.tsv",
        "../nanopore/albacore_2.2.7.castxb6.promethion.methylation.ref_summary.tsv",
        "../nanopore/albacore_2.2.7.castxb6.promethion.methylation.alt_summary.tsv",
    output:
        "../RData/paired_DSS.RData"
    script:
        "dss_paired.R"

rule compare_haplotype_methylation:
    input:
        meth = "../nanopore/{sample}.methylation.sorted.by_site.tsv",
        phase = "../nanopore/{sample}.phased.tsv",
        region = "../genome_data/ICR_plot_regions_string.txt"
    output:
        "../nanopore/{sample}.compare_haplotype_methylation.tsv"
    shell:
        "python compare_haplotype_methylation.py -m {input.meth} -p {input.phase} -b 11 -o 5 -r $(cat {input.region}) > {output}"

rule build_dmrlist:
    input:
        "../RData/paired_DSS.RData",
        "../genome_data/Mus_musculus.GRCm38_90.chr.genes.tsv",
    output:
        "../tables/dss_dmrlist.csv",
    script:
        "build_dss_dmrlist.R"

rule build_methylation_df:
    input:
        "../bisulfite/CpG_context_BC7.all.R1_val_1.fq_bismark_bt2_pe.summary.tsv",
        "../nanopore/albacore_1.2.2.b6xcast.methylation.summary.tsv",
    output:
        "../RData/nanopolish_df.RData",
        "../RData/bisulfite_df.RData",
    script:
        "build_methylation_df.R"
