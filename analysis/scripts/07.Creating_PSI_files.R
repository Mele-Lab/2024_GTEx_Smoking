#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to merge PSI files of each event type into one file per tissue 
# @software version: R=4.2.2

#Set path 
setwd(system("pwd", intern = T)) #If in linux
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
# setwd("..")

#Load libraries
library(optparse)
options(warn=-1)

# Parsing
parser <- OptionParser()
parser <- add_option(parser, opt_str=c("-t", "--tissue"), type="character",
                     dest="tissue",
                     help="Tissue")
options=parse_args(parser)
tissue=options$tissue
print(tissue)

# 1. Gene annotation
gene_annotation <- read.csv("data/public/gene_annotation.csv") #PC and lincRNA genes (no PAR)

# 2. Selected tissue samples
metadata <- readRDS(paste0("tissues/",tissue,"/metadata.rds"))

# 3. PSI values
splicing_events <- c("SE","MX","AF","AL","A5","A3","RI")

# Read psi value per sample and splicing event
psi <- lapply(splicing_events, function(e)
read.delim(paste0("SUPPA/PSI_values/",tissue,"/",tissue,".",e,".psi")))
names(psi) <- splicing_events

# Concatenate
PSI <- do.call(rbind.data.frame, psi)
# There are some 1.000002 and 0.99999999 values in the .tpm file that are read and seen as 1
# but internally are not changed. We should round those values
PSI <- round(PSI,2)

# Fix rownames
rownames(PSI) <- sapply(rownames(PSI), function(event_id)
  paste(unlist(strsplit(event_id,split = "\\."))[2:3],collapse = "."))
PSI$ensembl.id <- sapply(rownames(PSI), function(event_id)
  unlist(strsplit(event_id,split = ";"))[[1]])

# Keep splicing events that affect protein-coding and lincRNA genes
PSI <- PSI[PSI$ensembl.id %in% gene_annotation$ensembl.id,]

# Keep tissue samples with no missing metadata
colnames(PSI) <- gsub("\\.","-", colnames(PSI))
PSI <- PSI[,metadata$Sample] #I already did this in a previous code, but I do it again to double check and to remove the variabel gene_annotation

print("Starting with tpm")
# 4. Splicing events TPM values
# Read TPM value per sample and splicing event
tpm <- lapply(splicing_events, function(e)
  read.delim(paste0("SUPPA/PSI_values/",tissue,"/",tissue,".",e,".tpm")))
names(tpm) <- splicing_events

# Concatenate --
TPM <- do.call(rbind.data.frame, tpm)

# Fix rownames
rownames(TPM) <- sapply(rownames(TPM), function(event_id)
  paste(unlist(strsplit(event_id,split = "\\."))[2:3],collapse = "."))
TPM$ensembl.id <- sapply(rownames(TPM), function(event_id)
  unlist(strsplit(event_id,split = ";"))[[1]])

# Keep splicing events that affect protein-coding and lincRNA genes --
TPM <- TPM[TPM$ensembl.id %in% gene_annotation$ensembl.id,]

# Keep tissue samples (AFR and EUR sampes with no missing covariate metadata) --
colnames(TPM) <- gsub("\\.","-", colnames(TPM))
TPM <- TPM[,metadata$Sample]

# 5. Save data ####
saveRDS(PSI, paste0("tissues/",tissue,"/psi_splicing_events.rds"))
saveRDS(TPM, paste0("tissues/",tissue,"/tpm_splicing_events.rds"))
