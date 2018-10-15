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
        "../bismark_genome/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
        "../bismark_genome/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa"
    shell:
        "cd $(dirname {input}) && "
        "bismark_genome_preparation --bowtie2 ./ && "
        "cd ../scripts/"

rule bismark_align:
    input:
        "../bismark_genome/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
        "../bismark_genome/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa",
        r1 = "../bisulfite/{sample}_R1_val_1.fq.gz",
        r2 = "../bisulfite/{sample}_R2_val_2.fq.gz",
    output:
        log = "../bisulfite/{sample}.bismark.log",
        bam = "../bisulfite/{sample}.bismark_bt2_pe.bam",
    params:
        basename = lambda wildcards, output: "../bisulfite/{}".format(
            wildcards.sample),
        genome = "../bismark_genome/",
    threads:
        4
    shell:
        "bismark --gzip --bam --bowtie2 -p {threads} -B {params.basename} "
        "{params.genome} -1 {input.r1} -2 {input.r2} &> {output.log}"

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

rule merge_bisulfite_genome1:
    input:
        "../bisulfite/B6CastF1_1.bismark_bt2_pe.genome1.txt.gz",
        "../bisulfite/B6CastF1_2.bismark_bt2_pe.genome1.txt.gz",
        "../bisulfite/B6CastF1_5.bismark_bt2_pe.genome1.txt.gz",
        "../bisulfite/B6CastF1_6.bismark_bt2_pe.genome1.txt.gz"
    output:
        "../bisulfite/B6CastF1.combined_replicates.genome1.summary.tsv"
    shell:
        "python summarize_bisulfite_methylation.py {output} {input}"

rule merge_bisulfite_genome2:
    input:
        "../bisulfite/B6CastF1_1.bismark_bt2_pe.genome2.txt.gz",
        "../bisulfite/B6CastF1_2.bismark_bt2_pe.genome2.txt.gz",
        "../bisulfite/B6CastF1_5.bismark_bt2_pe.genome2.txt.gz",
        "../bisulfite/B6CastF1_6.bismark_bt2_pe.genome2.txt.gz"
    output:
        "../bisulfite/B6CastF1.combined_replicates.genome2.summary.tsv"
    shell:
        "python summarize_bisulfite_methylation.py {output} {input}"

rule merge_matched_bisulfite:
    input:
        "../bisulfite/B6CastF1_1.bismark_bt2_pe.txt.gz",
    output:
        "../bisulfite/B6CastF1_1.bismark_bt2_pe.summary.tsv"
    shell:
        "python summarize_bisulfite_methylation.py {output} {input}"
