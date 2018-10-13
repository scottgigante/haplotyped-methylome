library(GenomicRanges)
library(tidyverse)
library(data.table)


if (!file.exists("../RData/bisulfite_df.RData")) {
  bisulfite_df <- read_tsv("../bisulfite/CpG_context_TB1_PlacentaE14.5_WT_R1_trimmed.fq_bismark_bt2_pe.summary.tsv", col_names=c("chr","pos","percentMeth", "meth", "coverage"), col_types='ciddd') %>%
    arrange(chr, pos) %>%
    dplyr::rename(start=pos) %>%
    mutate(end=start,
           chr=sub("chr","",chr))
  save(bisulfite_df, file="../RData/bisulfite_df.RData")
}

if (!file.exists("../RData/nanopolish_df.RData")) {
  nanopolish_df <- read.table("../nanopore/albacore_1.2.2.b6xcast.methylation.summary.tsv",
                              header=TRUE, sep="\t", stringsAsFactors = FALSE) %>%
    mutate(start=start+1,
           end=end+2) %>%
    dplyr::rename(chr=chromosome,
                  percentMeth=methylated_frequency) %>%
    arrange(chr, start)
  save(nanopolish_df, file="../RData/nanopolish_df.RData")
}

if (!file.exists("../RData/knownGene.RData")) {
  gtf <- rtracklayer::import("../genome_data/ensembl_GRCm38.98.chr.gtf")
  knownGene <- gtf %>%
    as_data_frame() %>%
    rename(chr=seqnames) %>%
    select(chr, start, end, strand, type, gene_id, gene_name, exon_number) %>%
    filter(type %in% c("exon", "five_prime_utr", "three_prime_utr")) %>%
    filter(end > start) %>%
    mutate(type=fct_recode(type, utr="three_prime_utr", utr="five_prime_utr")) %>%
    unique() %>%
    group_by(gene_name, start) %>%
    filter(end == min(end)) %>%
    ungroup() %>%
    arrange(gene_name, start, end) %>%
    makeGRangesFromDataFrame(keep.extra.columns=TRUE)
  save(knownGene, file="../RData/knownGene.RData")
}

if (!file.exists("../tables/dss_dmrlist.csv")) {
  load("../RData/paired_DSS.RData")
  combined_dmr <- dmr_parent %>%
    mutate(type="imprinted") %>%
    bind_rows(dmr_strain %>%
                mutate(type="strain")) %>%
    mutate(chr = as.character(chr)) %>%
    arrange(-abs(areaStat)) %>%
    mutate(dmr_start=start,
           dmr_end=end,
           start = dmr_start-5000,
           end = dmr_end + 5000,
           id=row_number()) %>%
    select(chr, start, end, everything()) %>%
    fuzzyjoin::genome_left_join(read_tsv("../genome_data/ensembl_GRCm38.98.chr.genes.tsv", col_types="ccccdd") %>%
                                  select(chr, start, end, gene_name)) %>%
    select(-starts_with("start"), -starts_with("end"), -chr.y) %>%
    select(id, chr=chr.x, start=dmr_start, end=dmr_end, everything()) %>%
    group_by(id) %>%
    mutate(overlapping_genes=paste0(gene_name, collapse=",")) %>%
    distinct(overlapping_genes, .keep_all = TRUE) %>%
    ungroup() %>%
    select(-gene_name)
  combined_dmr %>%
    write_csv("../tables/dss_dmrlist.csv")
}