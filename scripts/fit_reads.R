# Fit longitudinal methylation
# Usage: Rscript fit_reads.R ../nanopore/file_prefix ../RData/outdir
# Assumes existence of the following:
# - file_prefix.sorted.bam.summary.tsv
# - file_prefix.phased.tsv
# - file_prefix.methylation.sorted.by_read.tsv

suppressMessages(library(tidyverse))
suppressMessages(library(ggExtra))
suppressMessages(library(broom))
suppressMessages(library(data.table))

n_min <- 7
args = commandArgs(trailingOnly=TRUE)
infile <- args[1]
outdir <- args[2]
dir.create(outdir, showWarnings = FALSE)

load(file.path(outdir, "summary_df.RData"))
load(file.path(outdir, "haplotype_df.RData"))

fit_loess <- function(chr, read_name, start, percentMeth, ...) {
  x = data.frame(start=start, percentMeth=percentMeth)
  fit <- loess(percentMeth ~ start, x, control=loess.control(surface = "interpolate",
                                                             statistics = "none",
                                                             trace.hat = "approximate",
                                                             cell=0.5), ...)
  fit <- list(kd=fit$kd,
              start=min(start),
              end=max(start),
              chr=chr,
              read_name=read_name)
  fit
}

process_file <- function(filepath) {
  con = file(filepath, "r")
  reads=list()
  header = readLines(con, n = 1)
  
  # process first line
  line = readLines(con, n = 1)
  line <- strsplit(line, "\\t")[[1]]
  chr=line[1]
  read_name = line[4]
  pos=as.integer(line[2])
  meth = 1/(1+exp(-as.numeric(line[5])))
  i=1
  start=c(pos)
  percentMeth=c(meth)
  
  # process the rest
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    line <- strsplit(line, "\\t")[[1]]
    pos=as.integer(line[2])
    name = line[4]
    meth = 1/(1+exp(-as.numeric(line[5])))
    if (name != read_name) {
      # new line
      if (length(start) > n_min){
        len <- max(start) - min(start)
        span <- 0.1 + 8e-11*((max(-len + 1e5, 0)))^2
        try({
          reads[[i]] = fit_loess(chr, read_name, start, percentMeth, span=span)
          i <- i + 1
        })
      }
      chr=line[1]
      read_name = name
      start=c(pos)
      percentMeth=c(meth)
    } else {
      start=c(start, pos)
      percentMeth=c(percentMeth, meth)
    }
  }
  close(con)
  reads
}

fit_reads <- process_file(paste0(infile, ".methylation.sorted.by_read.tsv"))
save(fit_reads, file=file.path(outdir, "fit_reads.RData"))
# load(file.path(outdir, "fit_reads.RData"))

fit_reads_df <- data_frame(read_name=sapply(fit_reads, function(x) { x$read_name })) %>% 
  mutate(id=row_number()) %>%
  left_join(summary_df, by="read_name") %>%
  arrange(chr, start, end) %>%
  left_join(haplotype_df %>% select(read_name, genotype), by="read_name") %>%
  dplyr::select(chr, start, end, id, read_name, genotype, qual)
save(fit_reads_df, file=file.path(outdir, "fit_reads_df.RData"))
