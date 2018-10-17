# Fit longitudinal methylation
# Usage: Rscript fit_reads.R ../nanopore/file_prefix.phased.tsv ../RData/outdir/haplotype_df.RData

suppressMessages(library(tidyverse))

min_coverage <- 5
args = commandArgs(trailingOnly=TRUE)
infile <- args[1]
outfile <- args[2]
dir.create(outdir, showWarnings = FALSE)

haplotype_df <- read_tsv(infile, col_types='ccddiddi') %>%
  mutate(signal_coverage = signal_ref + signal_alt,
         signal_ratio = signal_ref / signal_coverage,
         base_coverage = base_ref + base_alt,
         base_ratio=base_ref / base_coverage) %>% 
  mutate(info=case_when(.$signal_coverage < min_coverage & .$base_coverage < min_coverage ~ "low_coverage",
                        .$signal_ratio < 0.5 & .$base_ratio < 0.5 ~ "pass",
                        .$signal_ratio > 0.5 & .$base_ratio > 0.5 ~ "pass",
                        .$base_coverage > 3*.$signal_coverage ~ "signal_low_coverage",
                        .$signal_coverage > 3*.$base_coverage ~ "base_low_coverage",
                        abs(0.5-.$signal_ratio)*3 < abs(0.5-.$base_ratio) ~ "signal_uncertain",
                        abs(0.5-.$base_ratio)*3 < abs(0.5-.$signal_ratio) ~ "base_uncertain",
                        TRUE ~ "fail"),
         base_genotype=ifelse(base_ratio>0.5, "ref", "alt"),
         signal_genotype=ifelse(signal_ratio>0.5, "ref", "alt"),
         genotype="fail",
         genotype=ifelse(startsWith(info, "signal") | info == "pass", base_genotype, genotype),
         genotype=ifelse(startsWith(info, "base"), signal_genotype, genotype)) %>%
  rename(read_name=read)
save(haplotype_df, file=outfile)
