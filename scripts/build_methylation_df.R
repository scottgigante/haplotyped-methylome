library(tidyverse)
library(data.table)

bisulfite_df <- read_tsv("../bisulfite/CpG_context_BC7.all.R1_val_1.fq_bismark_bt2_pe.summary.tsv", col_names=c("chr","pos","percentMeth", "meth", "coverage"), col_types='ciddd') %>%
  arrange(chr, pos) %>%
  dplyr::rename(start=pos) %>%
  mutate(end=start,
         chr=sub("chr","",chr))
save(bisulfite_df, file="../RData/bisulfite_df.RData")

nanopolish_df <- read.table("../nanopore/albacore_1.2.2.b6xcast.methylation.summary.tsv",
                            header=TRUE, sep="\t", stringsAsFactors = FALSE) %>%
  mutate(start=start+1,
         end=end+2) %>%
  dplyr::rename(chr=chromosome,
                percentMeth=methylated_frequency) %>%
  arrange(chr, start)
save(nanopolish_df, file="../RData/nanopolish_df.RData")
