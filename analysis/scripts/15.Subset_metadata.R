#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to subset the metadata in the lung for both expression and methylation
# @software version: R=4.2.2

#Set path
setwd(system("pwd", intern = T)) #If in linux
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")

for (tissue in c("BreastMammaryTissue", "ColonTransverse", "Lung", "MuscleSkeletal", "Ovary", "Prostate", "Testis", "WholeBlood")){
  print(tissue)
  expression_metadata <- readRDS(paste0("tissues/", tissue,"/metadata.rds"))
  methylation_metadata <- readRDS(paste0("tissues/", tissue,"/methylation_metadata.rds"))
  
  expression_metadata <- expression_metadata[expression_metadata$Donor %in% methylation_metadata$SUBJID,]
  methylation_metadata <- methylation_metadata[methylation_metadata$SUBJID %in% expression_metadata$Donor,]
  
  #We will model 118 samples in the lung
  table(expression_metadata$Smoking)
  
  saveRDS(expression_metadata, paste0("tissues/", tissue,"/metadata_subset.rds"))
  saveRDS(methylation_metadata, paste0("tissues/", tissue,"/methylation_metadata_subset.rds"))
}

