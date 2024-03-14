#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to split the DNA methylation data per tissue
# @software version: R=4.2.2

#Set path 
setwd(system("pwd", intern = T)) #If in linux
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
# setwd("..")

data_path <- "/gpfs/scratch/bsc83/bsc83535/GTEx/" #Path where I have downloaded the raw beta values from https://www.nature.com/articles/s41588-022-01248-z 
# data_path <- "~/Documents/mn4/scratch/bsc83/bsc83535/GTEx/" #Path where I have downloaded the raw beta values from https://www.nature.com/articles/s41588-022-01248-z 
print("Reading whole data")
Sys.time()
data <- read.table(paste0(data_path, "v9/Oliva/GSE213478_methylation_DNAm_noob_final_BMIQ_all_tissues_987.txt.gz")) #It takes one hour
Sys.time()

#split data into tissues based on metadata
metadata <- read.delim(paste0(data_path, "v9/Oliva/eGTExDNA_Pierce_Jan18.09-11-2021.tsv")) #Data downloaded from Oliva et al.
names(metadata)[1] <- "ID"
names(metadata)[15] <- "Tissue"

samples <- strsplit(data[1,2], ",")[[1]]
samples <- samples[2:length(samples)]

samples_2 <- sapply(samples, function(x) substr(x, 2, nchar(x)-1))
names(samples_2) <- NULL

tissues <- sapply(samples_2, function(x) metadata$Tissue[metadata$ID==x])

print("Splitting whole data")

n <- nrow(data)
library(parallel)
Sys.time()
splitted_data <- mclapply(2:n, function(x) strsplit(data[x,2], ",")[[1]], mc.cores=8) #We use mclapply to speed up the computation
Sys.time()
names(splitted_data) <- data[2:n, 1]

splitted_data_frame <- as.data.frame(do.call(rbind, splitted_data))
splitted_data_frame <- splitted_data_frame[,2:ncol(splitted_data_frame)] #to remove the first column, which is always "", maybe we could have read the data frame in another way to solve that
colnames(splitted_data_frame) <- samples_2

print("Adding smoking metadata")
library("readxl")
smoking <- read_excel("data/protected/exSmoker_annotation.xlsx") #814
smoking <- smoking[!smoking$MHSMKTP %in% c("Cigar", "Pipe", "Other"),] #792    8 cigar smokers or ex-smokers or unknown, and 6 pipe ex-smokers/unkown and 8 others (ex-smokers/unkown)
smoking <- smoking[!smoking$SmokerStatus=="unknown",] #718
smoking$Smoking <- NA
smoking$Smoking[smoking$SmokerStatus=="non-smoker"] <- 0 #never-smokers
smoking$Smoking[smoking$SmokerStatus=="smoker"] <- 2
smoking$Smoking[smoking$SmokerStatus=="ex-smoker"] <- 1
smoking <- smoking[,c("SUBJID", "Smoking")]

translations <- read.csv("data/public/tissue_abreviation.txt")
for(tissue in unique(tissues)){
  print(tissue)
  if(tissue=="Breast - Mammary Tissue"){
    tissue_name <- "BreastMammaryTissue"
  } else{
    tissue_name <- translations$tissue[translations$SMTSD==tissue]
  }
  trues <- tissues == tissue
  data_to_save <- splitted_data_frame[,trues]

  #Reading metadata provided by Oliva et al. with the PEER values
  oliva_metadata <- read.delim(paste0(data_path, "v9/Oliva/", tissue_name, ".covariates.txt"))
  rownames(oliva_metadata) <- oliva_metadata[,1]
  oliva_metadata <- oliva_metadata[,2:ncol(oliva_metadata)]
  oliva_metadata <- as.data.frame(t(oliva_metadata))
  oliva_metadata$SUBJID <- rownames(oliva_metadata)
  oliva_metadata$SUBJID <- sapply(oliva_metadata$SUBJID, function(x) gsub("\\.", "-", x))

  #Adding metadata on inferred ancestry:
  ancestry <- read.delim("data/protected/inferred_ancestry_838donors.txt")
  names(ancestry) <- c("SUBJID", "Ancestry")
  metadata_to_save <- merge(oliva_metadata, ancestry, by="SUBJID")
  
  #Adding the rest of donor's metadata
  donor_metadata <- read.delim("data/protected/GTEx_Subject_Phenotypes.GRU.txt.gz")
  donor_metadata <- donor_metadata[,c("SUBJID", "SEX", "AGE", "BMI", "TRISCHD", "DTHHRDY")]
  metadata_to_save <- merge(metadata_to_save, donor_metadata, by="SUBJID")
  saveRDS(metadata_to_save, file = paste0("tissues/", tissue_name, "/methylation_metadata_no_smoking.rds"))
  
  #Adding smoking metadata
  metadata_to_save_s <- merge(metadata_to_save, smoking, by="SUBJID") #In their paper for mQTLs they use 5 PEERs for tissues with less than 50 samples, and 20 for tissues with more than 50. But they only show correlations with top 3 PEERs
  
  saveRDS(metadata_to_save_s, file = paste0("tissues/", tissue_name, "/methylation_metadata.rds"))
  
  #Data to save only for all donors in the given tissue:
  colnames(data_to_save) <- sapply(colnames(data_to_save), function(x) paste0(strsplit(x, "-")[[1]][1:2], collapse = "-")) #Now the column name is the subject id
  data_to_save <- data_to_save[,colnames(data_to_save) %in% metadata_to_save$SUBJID]
  saveRDS(data_to_save, file = paste0("tissues/", tissue_name, "/methylation_data_no_smoking.rds"))
  
  #Data to save only for donors with smoking annotation:
  data_to_save <- data_to_save[,colnames(data_to_save) %in% metadata_to_save_s$SUBJID]
  saveRDS(data_to_save, file = paste0("tissues/", tissue_name, "/methylation_data.rds"))
}
