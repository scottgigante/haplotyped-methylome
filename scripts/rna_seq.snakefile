"""TODO:

fastq from ENA
"""
rule download_rnaseq:
    output:
        "../rna_seq/{sample}.fastq.gz"

rule trim_galore:
    input:
        r1 = "../{outdir}/{sample}_R1.fastq.gz",
        r2 = "../{outdir}/{sample}_R2.fastq.gz",
    output:
        r1 = "../{outdir}/{sample}_R1_val_1.fq.gz",
        r2 = "../{outdir}/{sample}_R2_val_2.fq.gz",
    params:
        outdir = lambda wildcards, output: "../{}/".format(wildcards.outdir)
    shell:
        "trim_galore --phred33 --fastqc --gzip --paired -o {params.outdir} {input.r1} {input.r2}"

rule snpsplit_create_path:
    input:
        "../genome_data/GRCm38_90.fa",
    output:
        "../snpsplit_prepare_genome/GRCm38_90.fa"
    shell:
        "cp {input} {output}"

rule snpsplit_prepare_genome:
    input:
        genome = "../snpsplit_prepare_genome/GRCm38_90.fa",
        vcf = "../genome_data/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf"
    output:
        "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr1.N-masked.fa",
        "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr2.N-masked.fa",
        "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr3.N-masked.fa",
        "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr4.N-masked.fa",
        "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr5.N-masked.fa",
        "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr6.N-masked.fa",
        "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr7.N-masked.fa",
        "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr8.N-masked.fa",
        "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr9.N-masked.fa",
        "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr10.N-masked.fa",
        "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr11.N-masked.fa",
        "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr12.N-masked.fa",
        "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr13.N-masked.fa",
        "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr14.N-masked.fa",
        "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr15.N-masked.fa",
        "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr16.N-masked.fa",
        "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr17.N-masked.fa",
        "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr18.N-masked.fa",
        "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr19.N-masked.fa",
        "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chrX.N-masked.fa",
        "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chrY.N-masked.fa",
        "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chrMT.N-masked.fa",
        "../snpsplit_prepare_genome/all_SNPs_CAST_EiJ_GRCm38.txt.gz",
    shell:
        "cd $(dirname {input.genome}) && "
        "SNPsplit_genome_preparation --vcf_file {input.vcf} --strain CAST_EiJ --reference_genome ./ && "
        "cd ../scripts/"

rule clean_up_snpsplit_genome:
    input:
        fasta = ("../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr1.N-masked.fa",
                 "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr2.N-masked.fa",
                 "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr3.N-masked.fa",
                 "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr4.N-masked.fa",
                 "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr5.N-masked.fa",
                 "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr6.N-masked.fa",
                 "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr7.N-masked.fa",
                 "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr8.N-masked.fa",
                 "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr9.N-masked.fa",
                 "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr10.N-masked.fa",
                 "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr11.N-masked.fa",
                 "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr12.N-masked.fa",
                 "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr13.N-masked.fa",
                 "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr14.N-masked.fa",
                 "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr15.N-masked.fa",
                 "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr16.N-masked.fa",
                 "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr17.N-masked.fa",
                 "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr18.N-masked.fa",
                 "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chr19.N-masked.fa",
                 "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chrX.N-masked.fa",
                 "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chrY.N-masked.fa",
                 "../snpsplit_prepare_genome/CAST_EiJ_N-masked/chrMT.N-masked.fa"),
        snp = "../snpsplit_prepare_genome/all_SNPs_CAST_EiJ_GRCm38.txt.gz"
    params:
        tempdir = "../snpsplit_prepare_genome/"
    output:
        fasta = "../genome_data/GRCm38.Cast_N-masked.fa",
        snp = "../genome_data/all_SNPs_CAST_EiJ_GRCm38.txt.gz"
    shell:
        "cat {input.fasta} > {output.fasta} && "
        "mv {input.snp} {output.snp} && "
        "rm -r {params.tempdir}"

rule build_hisat2:
    input:
        "../genome_data/{genome}.fa"
    params:
        base = lambda wildcards, output: "../genome_data/{}".format(
            wildcards.genome)
    output:
        "../genome_data/{genome}.1.ht2"
    shell:
        "hisat2-build {input} {params.base}"

rule map_hisat2:
    input:
        genome = "../genome_data/GRCm38.Cast_N-masked.fa",
        index = "../genome_data/GRCm38.Cast_N-masked.1.ht2",
        r1 = "../rna_seq/{sample}_R1_val_1.fq.gz",
        r2 = "../rna_seq/{sample}_R2_val_2.fq.gz",
    output:
        "../rna_seq/{sample}.hisat2.bam"
    threads:
        16
    shell:
        "hisat2 -p {threads} --no-softclip -x {input.genome} -1 {input.r1} -2 {input.r2} | "
        "samtools sort -@ {threads} -T {output}.samtools.tmp -o {output}"

rule snp_split_hisat2:
    input:
        snp = "../genome_data/all_SNPs_CAST_EiJ_GRCm38.txt.gz",
        bam = "../rna_seq/{sample}.hisat2.bam"
    output:
        "../rna_seq/{sample}.hisat2.genome1.bam",
        "../rna_seq/{sample}.hisat2.genome2.bam",
    shell:
        "SNPsplit --snp_file {input.snp} --paired {input.bam}"

strains = ['B6Cast', 'CastB6']
samples = [2, 3, 4, 5]
rule samtools_depth:
    input:
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
        index = ("../rna_seq/B6CastF1_2.hisat2.genome1.bam.bai",
                 "../rna_seq/B6CastF1_2.hisat2.genome2.bam.bai",
                 "../rna_seq/B6CastF1_3.hisat2.genome1.bam.bai",
                 "../rna_seq/B6CastF1_3.hisat2.genome2.bam.bai",
                 "../rna_seq/B6CastF1_4.hisat2.genome1.bam.bai",
                 "../rna_seq/B6CastF1_4.hisat2.genome2.bam.bai",
                 "../rna_seq/B6CastF1_5.hisat2.genome1.bam.bai",
                 "../rna_seq/B6CastF1_5.hisat2.genome2.bam.bai",
                 "../rna_seq/CastB6F1_2.hisat2.genome1.bam.bai",
                 "../rna_seq/CastB6F1_2.hisat2.genome2.bam.bai",
                 "../rna_seq/CastB6F1_3.hisat2.genome1.bam.bai",
                 "../rna_seq/CastB6F1_3.hisat2.genome2.bam.bai",
                 "../rna_seq/CastB6F1_4.hisat2.genome1.bam.bai",
                 "../rna_seq/CastB6F1_4.hisat2.genome2.bam.bai",
                 "../rna_seq/CastB6F1_5.hisat2.genome1.bam.bai",
                 "../rna_seq/CastB6F1_5.hisat2.genome2.bam.bai")
    output:
        "../rna_seq/all_runs_with_reverse_coverage.tsv"
    shell:
        "samtools depth {input.bam} > {output}"
