suppressMessages(suppressPackageStartupMessages(library(tidyverse)))
suppressMessages(suppressPackageStartupMessages(library(data.table)))

bisulfite_df <- read_tsv("../bisulfite/B6CastF1_1_pe.summary.tsv", col_names=c("chr","pos","percentMeth", "meth", "coverage"), col_types='ciddd') %>%
  arrange(chr, pos) %>%
  dplyr::rename(start=pos) %>%
  mutate(end=start,
         chr=sub("chr","",chr))
save(bisulfite_df, file="../RData/bisulfite_df.RData")

nanopolish_df <- read.table("../nanopore/b6xcast.minion.methylation.summary.tsv",
                            header=TRUE, sep="\t", stringsAsFactors = FALSE) %>%
  mutate(start=start+1,
         end=end+2) %>%
  dplyr::rename(chr=chromosome,
                percentMeth=methylated_frequency) %>%
  arrange(chr, start)
save(nanopolish_df, file="../RData/nanopolish_df.RData")
