library(tidyverse)

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
  fuzzyjoin::genome_left_join(read_tsv("../genome_data/Mus_musculus.GRCm38_90.chr.genes.tsv", col_types="ccccdd") %>%
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