#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  cat("Usage:\n\t ./arch_finder.R <African population code> <DAF cutoff> <output file>\n")
  quit()
}

afr_pop <- args[1]
afr_cutoff <- as.numeric(args[2])
output_file <- args[3]

library(tidyverse)
library(glue)
library(data.table)

sample_info <- read_csv("input/sample_info.csv")

outgroups <- c("bonobo", "chimp", "gorilla")
archaics <- c("AltaiNeandertal", "Vindija33.19", "Chagyrskaya")
pops <- unique(filter(sample_info, dataset == "1000genomes")$pop)
superpops <- unique(filter(sample_info, dataset %in% c("1000genomes", "SGDP"))$superpop) %>% c("Africa2", "AFR2")

if (! afr_pop %in% pops)
  stop("Specified African population missing in the frequency table")

for (chr in 1:22) {
  vcf_path <- glue("/mnt/sequencedb/gendivdata/4_processed/ALT_freqs_from_merged_VCF/af_high_sgdp_g1000_apes/freq_var_chr{chr}.tab.gz")

  all_sites <-
    fread(vcf_path, na.strings = ".") %>%
    select(
      chrom=`#CHROM`, pos=POS,
      one_of(outgroups), one_of(archaics),
      UstIshim = Ust_Ishim, Loschbour, Stuttgart = LBK,
      one_of(pops), one_of(superpops)
    )

  info_sites <- all_sites[
    bonobo == 0 & chimp == 0 & gorilla == 0 &
    AltaiNeandertal == 1 & Vindija33.19 == 1 & Chagyrskaya == 1 &
    all_sites[[afr_pop]] <= afr_cutoff
  ]

  fwrite(info_sites, output_file, append = TRUE, sep = "\t")
}
