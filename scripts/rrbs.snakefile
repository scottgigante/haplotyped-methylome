"""TODO:

fastq from ENA
"""
rule download_bisulfite_B6Cast1:
    output:
        "../bisulfite/B6CastF1_1_R{file}.fastq.gz"
    params:
        url = lambda wildcards, output: "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR263/002/ERR2639372/ERR2639372_{}.fastq.gz".format(
            wildcards.file),
    shell:
        "wget -q -O {output} {params.url}"

rule download_bisulfite_B6Cast2:
    output:
        "../bisulfite/B6CastF1_2_R{file}.fastq.gz"
    params:
        url = lambda wildcards, output: "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR263/003/ERR2639373/ERR2639373_{}.fastq.gz".format(
            wildcards.file),
    shell:
        "wget -q -O {output} {params.url}"

rule download_bisulfite_B6Cast5:
    output:
        "../bisulfite/B6CastF1_5_R{file}.fastq.gz"
    params:
        url = lambda wildcards, output: "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR263/007/ERR2639377/ERR2639377_{}.fastq.gz".format(
            wildcards.file),
    shell:
        "wget -q -O {output} {params.url}"

rule download_bisulfite_B6Cast6:
    output:
        "../bisulfite/B6CastF1_6_R{file}.fastq.gz"
    params:
        url = lambda wildcards, output: "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR263/009/ERR2639379/ERR2639379_{}.fastq.gz".format(
            wildcards.file),
    shell:
        "wget -q -O {output} {params.url}"

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
        "../bismark_genome/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa",
    params:
        log = "../bismark_genome/Bisulfite_Genome.log"
    shell:
        "cd $(dirname {input}) && "
        "bismark_genome_preparation --bowtie2 ./ &> {params.log} && "
        "cd ../scripts/"

rule bismark_align:
    input:
        "../bismark_genome/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
        "../bismark_genome/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa",
        r1 = "../bisulfite/{sample}_R1_val_1.fq.gz",
        r2 = "../bisulfite/{sample}_R2_val_2.fq.gz",
    output:
        bam = "../bisulfite/{sample}_pe.bam",
    params:
        log = lambda wildcards, output: "../bisulfite/{}.bismark.log".format(
            wildcards.sample),
        basename = lambda wildcards, output: "{}".format(
            wildcards.sample),
        genome = "../bismark_genome/",
    threads:
        8
    shell:
        "bismark --gzip --bam --bowtie2 -p {threads} -B {params.basename} "
        "-o ../bisulfite/ {params.genome} -1 {input.r1} -2 {input.r2} &> {params.log}"

rule sort_bisulfite:
    input:
        "../bisulfite/{sample}_pe.bam",
    output:
        "../bisulfite/{sample}_pe.sorted.bam",
    threads:
        8
    shell:
        "samtools sort -n -@ {threads} -T {input}.samtools.tmp -o {output} {input}"

rule snpsplit_bismark:
    input:
        snp = "../genome_data/all_SNPs_CAST_EiJ_GRCm38.txt.gz",
        bam = "../bisulfite/{sample}_pe.sorted.bam",
    output:
        genome1 = "../bisulfite/{sample}_pe.sorted.genome1.bam",
        genome2 = "../bisulfite/{sample}_pe.sorted.genome2.bam",
    params:
        log = lambda wildcards, output: "../bisulfite/{}.snpsplit.log".format(
            wildcards.sample),
    shell:
        "SNPsplit --paired --bisulfite --snp_file {input.snp} --no_sort {input.bam} -o ../bisulfite/ &> {params.log}"

rule bismark_extract:
    input:
        "../bisulfite/{sample}.bam",
    output:
        "../bisulfite/{sample}.bedGraph.gz",
        "../bisulfite/{sample}.bismark.cov.gz",
        txt = "../bisulfite/CpG_context_{sample}.txt.gz",
    params:
        log = lambda wildcards, output: "../bisulfite/{}.bismark.log".format(
            wildcards.sample),
    threads:
        4
    shell:
        "bismark_methylation_extractor --ignore 13 --paired-end --multicore {threads} "
        "--comprehensive --merge_non_CpG --report --output $(dirname {output.txt}) --gzip "
        "--bedGraph {input} &> {params.log}"

rule merge_bisulfite_genome1:
    input:
        "../bisulfite/CpG_context_B6CastF1_1_pe.sorted.genome1.txt.gz",
        "../bisulfite/CpG_context_B6CastF1_2_pe.sorted.genome1.txt.gz",
        "../bisulfite/CpG_context_B6CastF1_5_pe.sorted.genome1.txt.gz",
        "../bisulfite/CpG_context_B6CastF1_6_pe.sorted.genome1.txt.gz"
    output:
        "../bisulfite/B6CastF1.combined_replicates.genome1.summary.tsv"
    shell:
        "python summarize_bisulfite_methylation.py {output} {input}"

rule merge_bisulfite_genome2:
    input:
        "../bisulfite/CpG_context_B6CastF1_1_pe.sorted.genome2.txt.gz",
        "../bisulfite/CpG_context_B6CastF1_2_pe.sorted.genome2.txt.gz",
        "../bisulfite/CpG_context_B6CastF1_5_pe.sorted.genome2.txt.gz",
        "../bisulfite/CpG_context_B6CastF1_6_pe.sorted.genome2.txt.gz"
    output:
        "../bisulfite/B6CastF1.combined_replicates.genome2.summary.tsv"
    shell:
        "python summarize_bisulfite_methylation.py {output} {input}"

rule merge_matched_bisulfite:
    input:
        "../bisulfite/CpG_context_B6CastF1_1_pe.txt.gz",
    output:
        "../bisulfite/B6CastF1_1_pe.summary.tsv"
    shell:
        "python summarize_bisulfite_methylation.py {output} {input}"
