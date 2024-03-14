#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Enrichment on TFBS for differentially methylated positions/loci
# @software version: R=4.2.2

library(dplyr)
library(tidyr)
library(data.table)

#Set path 
setwd(system("pwd", intern = T)) #If in linux
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
# setwd("..")

print("Reading CpG annotation")
Sys.time()
first_dir <- "/gpfs/"
# first_dir <- "~/Documents/mn4/"
data_path <- paste0(first_dir, "scratch/bsc83/bsc83535/GTEx/v9/Oliva/")
annotation <- read.csv(paste0(data_path, "GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv"))
Sys.time()

print("Reading chipatlas annotation")
Sys.time()
# chip <- read.delim(paste0(data_path, "remap2022_nr_macs2_hg19_v1_0.bed.gz"), header = F) #It takes almost 20 minutes
chip <- fread(paste0(data_path, "Oth.Lng.05.AllAg.AllCell.bed"), header = F)
Sys.time()

dim(chip)

#preprocess remap data
colnames(chip) <- c('chrom','start', 'end', 'to_split', 'score', 'strand', 'start_summit', 'end_summit', 'coord')
chip$TF <- gsub('.*=','', gsub('%.*', '', chip$to_split))
# chip$chrom <- sapply(chip$chrom, function(x) gsub("\\d+:\\s*", "", x))

#preprocess CpG annotation
ann_bed <- annotation[!is.na(annotation$MAPINFO) & !is.na(annotation$CHR),] %>%
  dplyr::select(chrom=CHR, start=MAPINFO, end=MAPINFO, name=IlmnID) %>% distinct()
ann_bed$chrom <- paste0('chr',ann_bed$chrom)
# head(ann_bed)

Sys.time()
print("About to merge")

### overlap TFBS and cpgs tested #### 
library(valr)
chip_1 <- chip[,c(1:3, 10)]
chip_cpgs <- bed_intersect(ann_bed, chip_1, suffix = c("_ann", "_tfbs"))

dim(chip_cpgs)
write.csv(chip_cpgs, paste0(data_path, "chip_seq_processed.csv"))
