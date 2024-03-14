#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to downsample
# @software version: R=4.2.2


#Set path 
# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")

tissues <- c("BreastMammaryTissue", "ColonTransverse", "Lung", "MuscleSkeletal", "Ovary", "Prostate", "Testis", "WholeBlood")
metadata <- list()
for(tissue in tissues){
  print(tissue)
  metadata_s <- readRDS(paste0("tissues/", tissue, "/methylation_metadata.rds"))
  metadata[[tissue]] <- metadata_s
}

samples <- matrix(1,8,3)
for(i in 1:length(metadata)){
  samples[i,] <- table(metadata[[i]]$Smoking)
}

apply(samples,2,min) #I will downsample to 11-11, 25-25, 35-35, 50-50 

set.seed(12345)
for(tissue in tissues){
  print(tissue)
  dir.create(paste0("Downsampling/11/", tissue, "/"), recursive = T)
  for(iter in 1:50){
    subset_ex <- metadata[[tissue]][metadata[[tissue]]$Smoking==1,"SUBJID"]
    subset_never <- sample(metadata[[tissue]][metadata[[tissue]]$Smoking==0,"SUBJID"], 11, replace=F)
    subset_smoker <- sample(metadata[[tissue]][metadata[[tissue]]$Smoking==2,"SUBJID"], 11, replace=F)
    subset <- metadata[[tissue]][metadata[[tissue]]$SUBJID %in% c(subset_never, subset_ex, subset_smoker),]
    saveRDS(subset, paste0("Downsampling/11/", tissue, "/metadata_downsampling_methylation_", iter, ".rds"))
    if(sum(metadata[[tissue]]$Smoking==2)>=50 & sum(metadata[[tissue]]$Smoking==0)>=50){
      dir.create(paste0("Downsampling/50/", tissue, "/"), recursive = T)
      subset_never <- sample(metadata[[tissue]][metadata[[tissue]]$Smoking==0,"SUBJID"], 50, replace=F)
      subset_smoker <- sample(metadata[[tissue]][metadata[[tissue]]$Smoking==2,"SUBJID"], 50, replace=F)
      subset <- metadata[[tissue]][metadata[[tissue]]$SUBJID %in% c(subset_never, subset_ex, subset_smoker),]
      saveRDS(subset, paste0("Downsampling/50/", tissue, "/metadata_downsampling_methylation_", iter, ".rds"))
    }
    if(sum(metadata[[tissue]]$Smoking==2)>=35 & sum(metadata[[tissue]]$Smoking==0)>=35){
      dir.create(paste0("Downsampling/35/", tissue, "/"), recursive = T)
      subset_never <- sample(metadata[[tissue]][metadata[[tissue]]$Smoking==0,"SUBJID"], 35, replace=F)
      subset_smoker <- sample(metadata[[tissue]][metadata[[tissue]]$Smoking==2,"SUBJID"], 35, replace=F)
      subset <- metadata[[tissue]][metadata[[tissue]]$SUBJID %in% c(subset_never, subset_ex, subset_smoker),]
      saveRDS(subset, paste0("Downsampling/35/", tissue, "/metadata_downsampling_methylation_", iter, ".rds"))
    }
    if(sum(metadata[[tissue]]$Smoking==2)>=25 & sum(metadata[[tissue]]$Smoking==0)>=25){
      dir.create(paste0("Downsampling/25/", tissue, "/"), recursive = T)
      subset_never <- sample(metadata[[tissue]][metadata[[tissue]]$Smoking==0,"SUBJID"], 25, replace=F)
      subset_smoker <- sample(metadata[[tissue]][metadata[[tissue]]$Smoking==2,"SUBJID"], 25, replace=F)
      subset <- metadata[[tissue]][metadata[[tissue]]$SUBJID %in% c(subset_never, subset_ex, subset_smoker),]
      saveRDS(subset, paste0("Downsampling/25/", tissue, "/metadata_downsampling_methylation_", iter, ".rds"))
    }
  }
}
