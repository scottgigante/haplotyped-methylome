rule download_B6_genome:
    output:
        "../genome_data/GRCm38_90.fa.gz",
    params:
        url = "ftp://ftp.ensembl.org/pub/release-90/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"
    shell:
        "wget -q {params.url} -O {output}"

rule download_gtf:
    output:
        "../genome_data/Mus_musculus.GRCm38_90.chr.gtf.gz"
    params:
        url = "ftp://ftp.ensembl.org/pub/release-90/gtf/mus_musculus/Mus_musculus.GRCm38.90.chr.gtf.gz"
    shell:
        "wget -q {params.url} -O {output}"

rule gunzip_fa:
    input:
        "{file}.fa.gz"
    output:
        "{file}.fa"
    shell:
        "gunzip {input}"

rule gunzip_vcf:
    input:
        "{file}.vcf.gz"
    output:
        "{file}.vcf"
    shell:
        "gunzip {input}"

rule gunzip_gtf:
    input:
        "{file}.gtf.gz"
    output:
        "{file}.gtf"
    shell:
        "gunzip {input}"

rule bwa_index:
    input:
        "../genome_data/{genome}.fa"
    output:
        "../genome_data/{genome}.fa.bwt",
        "../genome_data/{genome}.fa.ann",
        "../genome_data/{genome}.fa.amb",
        "../genome_data/{genome}.fa.pac",
        "../genome_data/{genome}.fa.sa",
    shell:
        "bwa index {input}"

rule samtools_index:
    input:
        "{file}.bam"
    output:
        "{file}.bam.bai"
    shell:
        "samtools index {input}"

rule download_vcf:
    output:
        "../genome_data/{strain}.mgp.v5.snps.dbSNP142.vcf.gz"
    params:
        url = lambda wildcards, output: "ftp://ftp-mouse.sanger.ac.uk/current_snps/strain_specific_vcfs/{}.mgp.v5.snps.dbSNP142.vcf.gz".format(
            wildcards.strain)
    shell:
        "wget -q {params.url} -O {output}"

rule mask_with_cast_genome:
    input:
        fasta = "../genome_data/{genome}.fa",
        vcf = "../genome_data/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf"
    output:
        "../genome_data/{genome}.CAST_masked.fa"
    shell:
        "python mask_genome_variants.py -i {input.fasta} -o {output} -v {input.vcf}"

rule ensembl_gtf_to_tsv:
    input:
        "../genome_data/Mus_musculus.GRCm38_90.chr.gtf",
    output:
        "../genome_data/Mus_musculus.GRCm38_90.chr.genes.tsv",
    shell:
        "Rscript ensembl_gtf_to_tsv.R {input}"

rule icr_plot_region_string:
    input:
        "../genome_data/ICR_plot_regions.tsv",
    output:
        "../genome_data/ICR_plot_regions_string.txt"
    shell:
        "python tsv_to_region.py {input} 10000 > {output}"

rule genome_cpg_coordinates:
    input:
        "../genome_data/GRCm38_90.fa"
    output:
        "../genome_data/GRCm38_90.cpg_coordinates.tsv"
    shell:
        "python count_cpg.py {input} > {output}"
