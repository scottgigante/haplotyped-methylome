rule bwa_mem_align:
    input:
        genome = "../genome_data/GRCm38_68.CAST_masked.fa",
        index = "../genome_data/GRCm38_68.CAST_masked.fa.bwt",
        reads = "../nanopore/{sample}.fa.gz"
    output:
        "../nanopore/{sample}.sorted.bam",
    threads:
        16
    shell:
        "bwa mem -t {threads} {input.genome} {input.reads} | samtools sort -T {sample}.samtools.tmp -@ {threads} -o {output}"

rule samtools_index_bam:
    input:
        "../nanopore/{sample}.sorted.bam"
    output:
        "../nanopore/{sample}.sorted.bam.bai"
    shell:
        "samtools index {input}"

rule nanopolish_index:
    input:
        reads = "../nanopore/{sample}.fa.gz",
        fast5 = "../nanopore/{sample}.fast5/"
    output:
        "../nanopore/{sample}.fa.gz.readdb"
    shell:
        "nanopolish index -d {fast5} {reads}"

rule nanopolish_methylation:
    input:
        genome = "../genome_data/GRCm38_68.CAST_masked.fa",
        bam = "../nanopore/{sample}.sorted.bam",
        index = "../nanopore/{sample}.sorted.bam.bai",
        reads = "../nanopore/{sample}.fa.gz",
        readdb = "../nanopore/{sample}.fa.gz.readdb"
    output:
        "../nanopore/{sample}.sorted.bam"
    threads:
        16
    shell:
        "nanopolish call-methylation -t {threads} -r {reads} -b {bam} -g {genome} > {output}"

rule nanopolish_phase:
    input:
        genome = "../genome_data/GRCm38_68.fa",
        bam = "../nanopore/{sample}.sorted.bam",
        index = "../nanopore/{sample}.sorted.bam.bai",
        reads = "../nanopore/{sample}.fa.gz",
        vcf = "../genome_data/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf",
        readdb = "../nanopore/{sample}.fa.gz.readdb"
    output:
        "../nanopore/{sample}.phased.sorted.bam"
    threads:
        16
    shell:
        "nanopolish phase-reads -t {threads} -r {reads} -b {bam} -g {genome} {vcf} | samtools sort -T {sample}.samtools.tmp -@ {threads} -o {output}"
