"""TODO:

fastq from ENA
"""
rule download_rrbs:
    output:
        "../bisulfite/{sample}.fastq.gz"

rule bismark_create_path:
    input:
        "../genome_data/GRCm38.Cast_N-masked.fa",
    output:
        "../bismark_genome/GRCm38.Cast_N-masked.fa"
    shell:
        "cp {input} {output}"

rule bismark_prepare_genome:
    input:
        "../bismark_genome/GRCm38.Cast_N-masked.fa"
    output:
        "../bismark_genome/Bisulfite_Genome/"
    shell:
        "cd $(dirname {input}) && "
        "bismark_genome_preparation --bowtie2 ./ && "
        "cd ../scripts/"

rule bismark_align:
    input:
        genome = "../bismark_genome/Bisulfite_Genome/",
        r1 = "../bisulfite/{sample}_R1_val_1.fq.gz",
        r2 = "../bisulfite/{sample}_R2_val_2.fq.gz",
    output:
        log = "../bisulfite/{sample}.bismark.log",
        bam = "../bisulfite/{sample}.bismark_bt2_pe.bam",
    params:
        basename = lambda wildcards, output: "../bisulfite/{}".format(
            wildcards.sample)
    threads:
        4
    shell:
        "bismark --gzip --bam --bowtie2 -p {threads} -B {params.basename} "
        "$(dirname {input.genome}) -1 {input.r1} -2 {input.r2} &> {output.log}"

rule snpsplit_bismark:
    input:
        snp = "../genome_data/all_SNPs_CAST_EiJ_GRCm38.txt.gz",
        bam = "../bisulfite/{sample}.bismark_bt2_pe.bam",
    output:
        "../bisulfite/{sample}.bismark_bt2_pe.genome1.bam",
        "../bisulfite/{sample}.bismark_bt2_pe.genome2.bam",
        log = "../bisulfite/{sample}.snpsplit.log",
    shell:
        "SNPsplit --paired --bisulfite --snp_file {input.snp} {input.bam} &> {output.log}"

rule bismark_extract:
    input:
        "../bisulfite/{sample}.bam",
    output:
        "../bisulfite/{sample}.txt.gz",
        log = "../bisulfite/{sample}.bismark.log",
    threads:
        2
    shell:
        "bismark_methylation_extractor --ignore 13 --paired-end --multicore {threads} "
        "--comprehensive --merge_non_CpG --report --output methylation_extractor --gzip "
        "--bedGraph {input} &> {output.log}"
