library(magrittr)
library(dplyr)
library(readr)

gtf <- rtracklayer::import("../genome_data/ensembl_GRCm38.98.chr.gtf")
gene_list <- gtf %>%
  as_data_frame() %>%
  rename(chr=seqnames) %>%
  select(chr, start, end, strand, gene_id, gene_name) %>%
  group_by(chr, gene_id, gene_name, strand) %>%
  summarise(start=min(start),
            end=max(end))
gene_list %>%
  write_tsv("../genome_data/ensembl_GRCm38.98.chr.genes.tsv")
