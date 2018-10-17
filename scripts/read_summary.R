# Fit longitudinal methylation
# Usage: Rscript fit_reads.R ../nanopore/file_prefix.sorted.bam.summary.tsv ../RData/outdir/summary_df.RData
# Assumes existence of the following:
# - file_prefix.sorted.bam.summary.tsv

suppressMessages(library(tidyverse))

args = commandArgs(trailingOnly=TRUE)
infile <- args[1]
outfile <- args[2]

summary_df <- read_tsv(infile,
                       col_names=c("read_name", "chr", "start", "end", "qual"), 
                       col_types='cciid',
                       skip=1) %>% 
  unique()
save(summary_df, file=outfile)
