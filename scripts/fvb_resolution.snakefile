rule call_variant_proportion_threeway:
    input:
        bam = "../nanopore/albacore_1.2.2.b6xcast.sorted.bam",
        bam_index = "../nanopore/albacore_1.2.2.b6xcast.sorted.bam.bai",
        suppdb = "../nanopore/albacore_1.2.2.b6xcast.sorted.bam.suppdb",
        phased_bam = "../nanopore/albacore_1.2.2.b6xcast.phased_sorted.bam",
        phased_bam_index = "../nanopore/albacore_1.2.2.b6xcast.phased_sorted.bam.bai",
        cast_vcf = "../genome_data/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf",
        fvb_vcf = "../genome_data/FVB_NJ.mgp.v5.snps.dbSNP142.vcf",
    output:
        "../nanopore/albacore_1.2.2.b6xcast.threeway_phased.tsv"
    shell:
        "python call_variant_proportion.py -b {input.bam} -p {input.phased_bam} -v {input.cast_vcf} {input.fvb_vcf} -o {output}"

rule fvb_resolution:
    input:
        "../nanopore/albacore_1.2.2.b6xcast.summary.tsv",
        "../nanopore/albacore_1.2.2.b6xcast.threeway_phased.tsv",
        "../plots/"
    output:
        "../genome_data/fvb_regions.txt",
    shell:
        "Rscript render_notebook.R ../nanopolish_threeway_haplotype_analysis.Rmd"

rule merge_threeway_variants:
    input:
        phase = "../nanopore/albacore_1.2.2.b6xcast.threeway_phased.tsv",
        summary = "../nanopore/albacore_1.2.2.b6xcast.summary.tsv",
        regions = "../genome_data/fvb_regions.txt",
    output:
        "../nanopore/albacore_1.2.2.b6xcast.phased.tsv"
    shell:
        "python merge_recombined_haplotypes.py -i {input.phase} -s {input.summary} -a alt2 -r $(cat {regions}) -o {output}"
