library(GenomicRanges)
library(tidyverse)
library(data.table)

if (!dir.exists("../plots")) dir.create("../plots")
if (!dir.exists("../plots/all_dmrs")) dir.create("../plots/all_dmrs")
if (!dir.exists("../plots/icrs")) dir.create("../plots/icrs")

if (!dir.exists("../tables")) dir.create("../tables")

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