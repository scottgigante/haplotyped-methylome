rule samtools_index_phased:
    input:
        "../nanopore/{sample}.phased.sorted.bam"
    output:
        "../nanopore/{sample}.phased.sorted.bam.bai"
    shell:
        "samtools index {input}"

rule suppdb_bam:
    input:
        bam = "../nanopore/{sample}.sorted.bam",
        bai = "../nanopore/{sample}.sorted.bam.bai",
    output:
        suppdb = "../nanopore/{sample}.sorted.bam.suppdb",
        summary = "../nanopore/{sample}.summary.tsv"
    shell:
        "python build_supplementary_db.py {input.bam} --summary {output.summary}"

rule suppdb_phased:
    input:
        bam = "../nanopore/{sample}.phased.sorted.bam",
        bai = "../nanopore/{sample}.phased.sorted.bam.bai",
    output:
        "../nanopore/{sample}.phased.sorted.bam.suppdb",
    shell:
        "python build_supplementary_db.py {input.bam}"

rule split_methylation_by_alignment:
    input:
        bam = "../nanopore/{sample}.sorted.bam",
        meth = "../nanopore/{sample}.methylation.tsv"
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


# def get_vcf(wildcards):
#     if "b6xcast" in wildcards.sample:
#         return ("../genome_data/FVB_NJ.mgp.v5.snps.dbSNP142.vcf",
#                 "../genome_data/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf")
#     else:
#         return "../genome_data/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf"


rule call_variant_proportion:
    input:
        bam = "../nanopore/{sample}.sorted.bam",
        phased_bam = "../nanopore/{sample}.phased.sorted.bam",
        vcf = "../genome_data/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf"
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
        "tail -n +2 {input.meth_split} | sort -k4,4 -k1,1 -k2n,2n | cat <(head -n 1 {input.meth}) - > {output}"

rule sort_by_site:
    input:
        meth_split = "../nanopore/{sample}.methylation.split_supplementary.tsv",
        meth = "../nanopore/{sample}.methylation.tsv"
    output:
        "../nanopore/{sample}.methylation.sorted.by_site.tsv"
    shell:
        "tail -n +2 {input.meth_split} | sort -k1,1 -k2n,2n -k4,4 | cat <(head -n 1 {input.meth}) - > {output}"

rule fit_reads:
    input:
        "../nanopore/{sample}.summary.tsv",
        "../nanopore/{sample}.phased.tsv",
        "../nanopore/{sample}.methylation.sorted.by_read.tsv"
    output:
        "../RData/{sample}/summary_df.RData",
        "../RData/{sample}/haplotype_df.RData",
        "../RData/{sample}/fit_reads.RData",
        "../RData/{sample}/fit_reads_df.RData"

rule split_methylation_by_haplotype:
    input:
        meth = "../nanopore/{sample}.methylation.sorted.by_site.tsv",
        phase = "../nanopore/{sample}.phased.tsv"
    output:
        "../nanopore/{sample}.methylation.sorted.by_site.tsv.ref.tsv",
        "../nanopore/{sample}.methylation.sorted.by_site.tsv.alt.tsv"
    shell:
        "python split_methylation_by_haplotype.py -m {input.meth} -p {input.phase}"

rule calculate_meth_ref_freq:
    input:
        "../nanopore/{sample}.methylation.sorted.by_site.tsv.ref.tsv",
    output:
        "../nanopore/{sample}.methylation.ref.summary.tsv"
    shell:
        "python calculate_methylation_frequency.py -i {input} -p > {output}"

rule calculate_meth_alt_freq:
    input:
        "../nanopore/{sample}.methylation.sorted.by_site.tsv.alt.tsv",
    output:
        "../nanopore/{sample}.methylation.alt.summary.tsv"
    shell:
        "python calculate_methylation_frequency.py -i {input} -p > {output}"

rule paired_dss:
    input:
        "../nanopore/albacore_1.2.2.b6xcast.methylation.ref.summary.tsv",
        "../nanopore/albacore_1.2.2.b6xcast.methylation.alt.summary.tsv",
        "../nanopore/albacore_2.2.7.castxb6.promethion.methylation.ref.summary.tsv",
        "../nanopore/albacore_2.2.7.castxb6.promethion.methylation.alt.summary.tsv",
    output:
        "../RData/paired_DSS.RData"
    script:
        "dss_paired.R"

rule compare_haplotype_methylation:
    input:
        meth = "../nanopore/{sample}.methylation.sorted.by_read.tsv",
        phase = "../nanopore/{sample}.phased.tsv",
        region = "../genome_data/ICR_plot_regions_string.txt"
    output:
        "../nanopore/{sample}.compare_haplotype_methylation.tsv"
    shell:
        "python compare_haplotype_methylation.py -m {input.meth} -p {input.phase} -b 11 -o 5 -r $(cat {input.region}) > {output}"
