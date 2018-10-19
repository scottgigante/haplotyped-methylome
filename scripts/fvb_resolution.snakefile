rule call_variant_proportion_threeway:
    input:
        bam = "../nanopore/b6xcast.minion.sorted.bam",
        bam_index = "../nanopore/b6xcast.minion.sorted.bam.bai",
        suppdb = "../nanopore/b6xcast.minion.sorted.bam.suppdb",
        phased_bam = "../nanopore/b6xcast.minion.phased_sorted.bam",
        phased_bam_index = "../nanopore/b6xcast.minion.phased_sorted.bam.bai",
        cast_vcf = ancient("../genome_data/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf"),
        fvb_vcf = ancient("../genome_data/FVB_NJ.mgp.v5.snps.dbSNP142.vcf"),
    output:
        "../nanopore/b6xcast.minion.threeway_phased.tsv"
    shell:
        "python call_variant_proportion.py -b {input.bam} -p {input.phased_bam} -v {input.cast_vcf} {input.fvb_vcf} -o {output}"

rule fvb_resolution:
    input:
        "../nanopore/b6xcast.minion.sorted.bam.summary.tsv",
        "../nanopore/b6xcast.minion.threeway_phased.tsv",
    output:
        "../genome_data/fvb_regions.txt",
        "../plots/rpart_fvb_resolution.png",
        "../notebooks/nanopolish_fvb_resolution.html",
    params:
        log = "../notebooks/nanopolish_fvb_resolution.log",
    shell:
        "Rscript render_notebook.R ../notebooks/nanopolish_fvb_resolution.Rmd &> {params.log}"

rule merge_threeway_variants:
    input:
        phase = "../nanopore/b6xcast.minion.threeway_phased.tsv",
        summary = "../nanopore/b6xcast.minion.sorted.bam.summary.tsv",
        regions = "../genome_data/fvb_regions.txt",
    output:
        "../nanopore/b6xcast.minion.phased.tsv"
    shell:
        "python merge_recombined_haplotypes.py -i {input.phase} -s {input.summary} -a alt2 -r $(cat {input.regions}) -o {output}"
