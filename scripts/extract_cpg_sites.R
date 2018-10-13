library(tidyverse)
dir <- "../txk/cpg_data/GRCm38/"
cg_df <- map_df(list.files(dir, pattern = "CG_.*tsv"), ~ 
                  read_tsv(paste(dir, ., sep="/"), col_names = c("chr", "start"), col_types="ci_") %>%
                  mutate(chr = str_remove(chr, "chr")))
save(cg_df, file = "../genome_data/CpG_coordinates_mm10.RData")
