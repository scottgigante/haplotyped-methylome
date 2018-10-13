# Detect DMRs in nanopore data
# Usage: Rscript 

suppressMessages(library(tidyverse))
suppressMessages(library(DSS))

args <- commandArgs(trailingOnly=TRUE)
meth_fn <- args[1]
outdir <- args[2]

ref_df <- read.table(paste0(meth_fn, ".ref.summary.tsv"),
                              header=TRUE, sep="\t", stringsAsFactors = FALSE) %>%
    mutate(start=start+1,
           end=end+2) %>%
    dplyr::rename(chr=chromosome,
                  percentMeth=methylated_frequency) %>%
    arrange(chr, start)
alt_df <- read.table(paste0(meth_fn, ".alt.summary.tsv"),
                              header=TRUE, sep="\t", stringsAsFactors = FALSE) %>%
    mutate(start=start+1,
           end=end+2) %>%
    dplyr::rename(chr=chromosome,
                  percentMeth=methylated_frequency) %>%
    arrange(chr, start)

bs <- makeBSseqData(dat = list(ref_df %>% select(chr, pos=start, N=called_sites, X=called_sites_methylated) %>%
                                 mutate(X=round(X)),
                               alt_df %>% select(chr, pos=start, N=called_sites, X=called_sites_methylated) %>%
                                 mutate(X=round(X))),
                    sampleNames=c("ref", "alt"))
dml <- DMLtest(BSobj = bs, group1="ref", group2="alt", equal.disp = TRUE, smoothing=TRUE)

dmr <- callDMR(dml, p.threshold=0.001)
save(dmr, file=file.path(outdir, "DSS_dmr.RData"))

