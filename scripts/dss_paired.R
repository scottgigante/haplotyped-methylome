suppressMessages(suppressPackageStartupMessages(library(tidyverse)))
suppressMessages(suppressPackageStartupMessages(library(DSS)))
suppressMessages(suppressPackageStartupMessages(library(bsseq)))

meth_fn <- "../nanopore/b6xcast.minion.methylation"

forward_ref_df <- read.table(paste0(meth_fn, ".ref_summary.tsv"),
                             header=TRUE, sep="\t", stringsAsFactors = FALSE) %>%
  mutate(start=start+1,
         end=end+2,
         sampleName="b6xcast.mat") %>%
  dplyr::rename(chr=chromosome,
                percentMeth=methylated_frequency) %>%
  arrange(chr, start)
forward_alt_df <- read.table(paste0(meth_fn, ".alt_summary.tsv"),
                             header=TRUE, sep="\t", stringsAsFactors = FALSE) %>%
  mutate(start=start+1,
         end=end+2,
         sampleName="b6xcast.pat") %>%
  dplyr::rename(chr=chromosome,
                percentMeth=methylated_frequency) %>%
  arrange(chr, start)

meth_fn <- "../nanopore/castxb6.promethion.methylation"

reverse_ref_df <- read.table(paste0(meth_fn, ".ref_summary.tsv"),
                              header=TRUE, sep="\t", stringsAsFactors = FALSE) %>%
    mutate(start=start+1,
           end=end+2,
         sampleName="castxb6.pat") %>%
    dplyr::rename(chr=chromosome,
                  percentMeth=methylated_frequency) %>%
    arrange(chr, start)
reverse_alt_df <- read.table(paste0(meth_fn, ".alt_summary.tsv"),
                              header=TRUE, sep="\t", stringsAsFactors = FALSE) %>%
    mutate(start=start+1,
           end=end+2,
         sampleName="castxb6.mat") %>%
    dplyr::rename(chr=chromosome,
                  percentMeth=methylated_frequency) %>%
    arrange(chr, start)

tidy_df <- function(df, sampleName) {
  Cov <- rlang::sym(paste0("Cov.",sampleName))
  M <- rlang::sym(paste0("M.",sampleName))
  df %>%
    select(chr, start, end, 
           !!Cov:=called_sites, 
           !!M:=called_sites_methylated)
}

make_bs <- function(df) {
  M <- df %>% 
    select(starts_with("M.")) %>%
    rename_all(function(x) gsub("M.", "", x)) %>%
    mutate_all(function(x) ifelse(is.na(x), 0, x)) %>%
    as.matrix()
  Cov <- df %>% 
    select(starts_with("Cov.")) %>%
    rename_all(function(x) gsub("Cov.", "", x)) %>%
    mutate_all(function(x) ifelse(is.na(x), 0, x)) %>%
    as.matrix()
  BSseq(chr = df$chr, 
        pos = df$start,
        Cov = Cov,
        M = M,
        sampleNames = colnames(M))
}

bs <- tidy_df(forward_ref_df, "b6xcast.mat") %>%
  full_join(tidy_df(forward_alt_df, "b6xcast.pat"),
            by=c("chr", "start", "end")) %>%
  full_join(tidy_df(reverse_ref_df, "castxb6.pat"),
            by=c("chr", "start", "end")) %>%
  full_join(tidy_df(reverse_alt_df, "castxb6.mat"),
            by=c("chr", "start", "end")) %>%
  group_by(chr, start) %>%
  summarise_if(is.numeric, sum) %>%
  ungroup() %>%
  make_bs()
pData(bs)$strain <- c("b6xcast", "b6xcast", "castxb6", "castxb6")
pData(bs)$parent <- c("mat", "pat", "pat", "mat")

loci.idx <- which(rowSums(getCoverage(bs, type="Cov")==0) == 0)
dml_parent <- DMLtest(BSobj = bs[loci.idx,], 
               group1=c("b6xcast.mat", "castxb6.mat"), 
               group2=c("b6xcast.pat", "castxb6.pat"), 
               equal.disp = TRUE, smoothing=TRUE)
dml_strain <- DMLtest(BSobj = bs[loci.idx,], 
               group1=c("b6xcast.mat", "castxb6.pat"), 
               group2=c("castxb6.mat", "b6xcast.pat"), 
               equal.disp = TRUE, smoothing=TRUE)

dmr_parent <- callDMR(dml_parent, p.threshold=1e-5, dis.merge=1500)
dmr_strain <- callDMR(dml_strain, p.threshold=1e-5, dis.merge=1500)

save(dml_parent, dml_strain, dmr_parent, dmr_strain, file="../RData/paired_DSS.RData")

