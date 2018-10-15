rule download_nanopore:
    output:
        "../nanopore/{sample}.fa.gz"

rule download_fast5:
    output:
        "../nanopore/{sample}.fast5/"

rule bwa_mem_align:
    input:
        genome = "../genome_data/GRCm38_90.CAST_masked.fa",
        index = "../genome_data/GRCm38_90.CAST_masked.fa.bwt",
        reads = "../nanopore/{sample}.fa.gz"
    output:
        "../nanopore/{sample}.sorted.bam",
    threads:
        16
    shell:
        "bwa mem -t {threads} {input.genome} {input.reads} | "
        "samtools sort -T {output}.samtools.tmp -@ {threads} -o {output}"

rule nanopolish_index:
    input:
        reads = "../nanopore/{sample}.fa.gz",
        fast5 = "../nanopore/{sample}.fast5/"
    output:
        "../nanopore/{sample}.fa.gz.readdb"
    shell:
        "nanopolish index -d {input.fast5} {input.reads}"

rule nanopolish_methylation:
    input:
        genome = "../genome_data/GRCm38_90.CAST_masked.fa",
        bam = "../nanopore/{sample}.sorted.bam",
        index = "../nanopore/{sample}.sorted.bam.bai",
        reads = "../nanopore/{sample}.fa.gz",
        readdb = "../nanopore/{sample}.fa.gz.readdb"
    output:
        "../nanopore/{sample}.methylation.tsv"
    threads:
        16
    shell:
        "nanopolish call-methylation -t {threads} -r {input.reads} -b {input.bam} -g {input.genome} > {output}"

rule nanopolish_phase:
    input:
        genome = "../genome_data/GRCm38_90.fa",
        bam = "../nanopore/{sample}.sorted.bam",
        index = "../nanopore/{sample}.sorted.bam.bai",
        reads = "../nanopore/{sample}.fa.gz",
        vcf = "../genome_data/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf",
        readdb = "../nanopore/{sample}.fa.gz.readdb"
    output:
        "../nanopore/{sample}.phased_sorted.bam"
    threads:
        16
    shell:
        "nanopolish phase-reads -t {threads} -r {input.reads} -b {input.bam} -g {input.genome} {input.vcf} | "
        "samtools sort -T {output}.samtools.tmp -@ {threads} -o {output}"
