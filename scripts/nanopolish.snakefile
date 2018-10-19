rule download_fast5_b6xcast_0503:
    output:
        temp("../nanopore/2017_05_03.b6xcast.minion.fast5.tar.gz")
    params:
        url = "ftp://ftp.sra.ebi.ac.uk/vol1/ERA153/ERA1533189/oxfordnanopore_native/2017_05_03_MOUSE_WGS_ONT.fast5.tar.gz",
    shell:
        "wget -q -O {output} {params.url}"

rule download_fast5_b6xcast_0512:
    output:
        temp("../nanopore/2017_05_12.b6xcast.minion.fast5.tar.gz")
    params:
        url = "ftp://ftp.sra.ebi.ac.uk/vol1/ERA153/ERA1533189/oxfordnanopore_native/2017_05_12_MOUSE_WGS_ONT.fast5.tar.gz",
    shell:
        "wget -q -O {output} {params.url}"

rule download_fast5_b6xcast_0525:
    output:
        temp("../nanopore/2017_05_25.b6xcast.minion.fast5.tar.gz")
    params:
        url = "ftp://ftp.sra.ebi.ac.uk/vol1/ERA153/ERA1533189/oxfordnanopore_native/2017_05_25_MOUSE_WGS_ONT.fast5.tar.gz",
    shell:
        "wget -q -O {output} {params.url}"

rule download_fast5_castxb6:
    output:
        temp("../nanopore/castxb6.promethion.fast5.tar.gz")
    params:
        url = "ftp://ftp.sra.ebi.ac.uk/vol1/ERA153/ERA1533189/oxfordnanopore_native/20180410_0355_BLEWITT_CASTB6_LSK109.fast5.tar.gz",
    shell:
        "wget -q -O {output} {params.url}"

rule download_fast5_b6:
    output:
        temp("../nanopore/b6.minion.fast5.tar.gz")
    params:
        url = "ftp://ftp.sra.ebi.ac.uk/vol1/ERA153/ERA1533189/oxfordnanopore_native/Black6_WGS_ONT.fast5.tar.gz",
    shell:
        "wget -q -O {output} {params.url}"

rule download_fast5_cast:
    output:
        temp("../nanopore/cast.minion.fast5.tar.gz")
    params:
        url = "ftp://ftp.sra.ebi.ac.uk/vol1/ERA153/ERA1533189/oxfordnanopore_native/Cast_WGS_ONT.fast5.tar.gz",
    shell:
        "wget -q -O {output} {params.url}"

rule fast5_md5:
    input:
        "../nanopore/{archive}.tar.gz",
        md5 = ancient("../md5/nanopore/{archive}.tar.gz.md5")
    output:
        temp("../nanopore/{archive}.tar.gz.md5_ok")
    shell:
        "md5sum -c {input.md5} && touch {output}"

rule untar:
    input:
        tar = "../nanopore/{archive}.tar.gz",
        md5 = "../nanopore/{archive}.tar.gz.md5_ok",
    output:
        directory("../nanopore/{archive}/")
    shell:
        "mkdir -p {output} && "
        "tar xzf {input.tar} -C {output} && "
        "find {output} -mindepth 2 -type d -links 2  -exec mv -t {output} {{}} \\+ && "
        "for i in $(find {output} -name '*.tar.gz'); do "
        "  tar -xzf $i -C {output} && rm $i; "
        "done"

rule merge_b6xcast:
    input:
        "../nanopore/2017_05_03.b6xcast.minion.fast5/",
        "../nanopore/2017_05_12.b6xcast.minion.fast5/",
        "../nanopore/2017_05_25.b6xcast.minion.fast5/",
    output:
        directory("../nanopore/b6xcast.minion.fast5/"),
    shell:
        "mkdir -p {output} && "
        "for i in {input}; do mv $i {output}; done"

rule albacore_minion:
    input:
        "../nanopore/{sample}.minion.fast5/",
    output:
        directory("../nanopore/{sample}.minion.albacore/workspace"),
        "../nanopore/{sample}.minion.albacore/sequencing_summary.txt",
    params:
        outdir = lambda wildcards, output: "../nanopore/{}.minion.albacore".format(
            wildcards.sample)
    threads:
        16
    shell:
        "read_fast5_basecaller.py -c r94_450bps_linear.cfg -o fastq -i {input} -r -s {params.outdir} -t {threads} -q 0"

rule albacore_promethion:
    input:
        "../nanopore/{sample}.promethion.fast5/",
    output:
        directory("../nanopore/{sample}.promethion.albacore/workspace"),
        "../nanopore/{sample}.minion.albacore/sequencing_summary.txt",
    params:
        outdir = lambda wildcards, output: "../nanopore/{}.promethion.albacore".format(
            wildcards.sample)
    threads:
        16
    shell:
        "read_fast5_basecaller.py -c r941_450bps_linear_prom.cfg -o fastq -i {input} -r -s {params.outdir} -t {threads} -q 0"

rule albacore_merge_output:
    input:
        "../nanopore/{sample}.albacore/workspace",
    output:
        "../nanopore/{sample}.fastq",
    shell:
        "find {input} -name '*.fastq' -exec cat {{}} \\+ > {output}"

rule bwa_mem_align:
    input:
        genome = "../genome_data/GRCm38_90.CAST_masked.fa",
        index = "../genome_data/GRCm38_90.CAST_masked.fa.bwt",
        reads = "../nanopore/{sample}.fastq"
    output:
        "../nanopore/{sample}.sorted.bam",
    threads:
        16
    log:
        "../nanopore/{sample}.sorted.bam.log"
    shell:
        "bwa mem -t {threads} {input.genome} {input.reads} | "
        "samtools sort -T {output}.samtools.tmp -@ {threads} -o {output} &> {log}"

rule nanopolish_index:
    input:
        reads = "../nanopore/{sample}.fastq",
        fast5 = "../nanopore/{sample}.fast5/"
    output:
        "../nanopore/{sample}.fastq.index.readdb"
    params:
        log = lambda wildcards, output: "{}.log".format(output)
    shell:
        "nanopolish index -d {input.fast5} {input.reads} &> {params.log}"

rule nanopolish_methylation:
    input:
        genome = "../genome_data/GRCm38_90.CAST_masked.fa",
        bam = "../nanopore/{sample}.sorted.bam",
        index = "../nanopore/{sample}.sorted.bam.bai",
        reads = "../nanopore/{sample}.fastq",
        readdb = "../nanopore/{sample}.fastq.index.readdb"
    output:
        "../nanopore/{sample}.methylation.tsv"
    threads:
        1
    shell:
        "nanopolish call-methylation -t 16 -r {input.reads} -b {input.bam} -g {input.genome} > {output}"

rule nanopolish_phase:
    input:
        genome = ancient("../genome_data/GRCm38_90.fa"),
        bam = "../nanopore/{sample}.sorted.bam",
        index = "../nanopore/{sample}.sorted.bam.bai",
        reads = "../nanopore/{sample}.fastq",
        vcf = ancient("../genome_data/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf"),
        readdb = "../nanopore/{sample}.fastq.index.readdb"
    output:
        "../nanopore/{sample}.phased_sorted.bam"
    threads:
        1
    shell:
        "nanopolish phase-reads -t 16 -r {input.reads} -b {input.bam} -g {input.genome} {input.vcf} | "
        "samtools sort -T {output}.samtools.tmp -@ {threads} -o {output}"

rule intermediate_download_bam:
    input:
        ancient("../genome_data/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf"),
        ancient("../genome_data/GRCm38_90.fa"),
        "../genome_data/GRCm38_90.CAST_masked.fa",
    output:
        temp("../nanopore/{sample}.intermediate_bam"),
    params:
        sample = lambda wildcards, output: wildcards.sample
    shell:
        "touch ../nanopore/{params.sample}.intermediate_bam && "
        "touch ../nanopore/{params.sample}.sorted.bam && "  # actually download
        "touch ../nanopore/{params.sample}.sorted.bam.bai"  # actually index


rule intermediate_download_nanopolish:
    input:
        "../nanopore/{sample}.intermediate_bam",
    params:
        sample = lambda wildcards, output: wildcards.sample
    output:
        temp("../nanopore/{sample}.intermediate_download"),
    shell:
        # actually download
        "touch ../nanopore/{params.sample}.phased_sorted.bam && "
        # actually download
        "touch ../nanopore/{params.sample}.methylation.tsv && "
        "touch ../nanopore/{params.sample}.intermediate_download"
