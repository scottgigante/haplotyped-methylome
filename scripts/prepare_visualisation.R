suppressMessages(suppressPackageStartupMessages(library(GenomicRanges)))
suppressMessages(suppressPackageStartupMessages(library(tidyverse)))
suppressMessages(suppressPackageStartupMessages(library(data.table)))

gtf <- rtracklayer::import("../genome_data/Mus_musculus.GRCm38_90.chr.gtf")
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

rna_seq <- read_tsv("../rna_seq/all_runs_with_reverse_coverage.tsv", 
                    col_names=c("chr", "pos", "fwd1_mat", "fwd1_pat", 
                                "fwd2_mat", "fwd2_pat", "fwd3_mat", "fwd3_pat", 
                                "fwd4_mat", "fwd4_pat", "rev1_pat", "rev1_mat", 
                                "rev2_pat", "rev2_mat", "rev3_pat", "rev3_mat", 
                                "rev4_pat", "rev4_mat"),
                    col_types = 'ciiiiiiiiiiiiiiiii') %>%
  mutate(fwd_mat=(fwd1_mat+fwd2_mat+fwd3_mat+fwd4_mat)/4,
         fwd_pat=(fwd1_pat+fwd2_pat+fwd3_pat+fwd4_pat)/4,
         rev_pat=(rev1_pat+rev2_pat+rev3_pat+rev4_pat)/4,
         rev_mat=(rev1_mat+rev2_mat+rev3_mat+rev4_mat)/4) %>% 
  select(chr, pos, fwd_mat, fwd_pat, rev_mat, rev_pat)
save(rna_seq, file="../RData/rna_seq_with_reverse.RData")