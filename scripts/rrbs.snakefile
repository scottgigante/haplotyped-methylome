"""TODO:

fastq from ENA
"""

rule bismark_prepare_genome:
    input:
        "../genome_data/GRCm38.Cast_N-masked.fa"
    params:
        tempdir = "../genome_data/GRCm38.Cast_N-masked"
    output:
