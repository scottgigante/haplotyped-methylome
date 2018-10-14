library(magrittr)
library(dplyr)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
gtf <- rtracklayer::import(args[1])
gene_list <- gtf %>%
  as_data_frame() %>%
  rename(chr=seqnames) %>%
  select(chr, start, end, strand, gene_id, gene_name) %>%
  group_by(chr, gene_id, gene_name, strand) %>%
  summarise(start=min(start),
            end=max(end))
gene_list %>%
  write_tsv(sub("gtf$", "genes.tsv", args[1]))
