#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to downsample
# @software version: R=4.2.2


#Set path 
# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")

metadata <- readRDS("output/metadata.rds")
# 
# samples <- matrix(1,46,3)
# for(i in 1:length(metadata)){
#   samples[i,] <- table(metadata[[i]]$Smoking)
# }
# 
# apply(samples, 2, min)

set.seed(12345)
for(tissue in names(metadata)){
  print(tissue)
  dir.create(paste0("Downsampling/35/", tissue, "/"), recursive = T)
  dir.create(paste0("Downsampling/11/", tissue, "/"), recursive = T)
  dir.create(paste0("Downsampling/25/", tissue, "/"), recursive = T)
  for(iter in 1:50){
    subset_ex <- metadata[[tissue]][metadata[[tissue]]$Smoking==1,"Donor"]
    subset_never <- sample(metadata[[tissue]][metadata[[tissue]]$Smoking==0,"Donor"], 35, replace=F)
    subset_smoker <- sample(metadata[[tissue]][metadata[[tissue]]$Smoking==2,"Donor"], 35, replace=F)
    subset <- metadata[[tissue]][metadata[[tissue]]$Donor %in% c(subset_never, subset_ex, subset_smoker),]
    saveRDS(subset, paste0("Downsampling/35/", tissue, "/metadata_downsampling_", iter, ".rds"))
    
    subset_never <- sample(metadata[[tissue]][metadata[[tissue]]$Smoking==0,"Donor"], 25, replace=F)
    subset_smoker <- sample(metadata[[tissue]][metadata[[tissue]]$Smoking==2,"Donor"], 25, replace=F)
    subset <- metadata[[tissue]][metadata[[tissue]]$Donor %in% c(subset_never, subset_ex, subset_smoker),]
    saveRDS(subset, paste0("Downsampling/25/", tissue, "/metadata_downsampling_", iter, ".rds"))
    
    subset_never <- sample(metadata[[tissue]][metadata[[tissue]]$Smoking==0,"Donor"], 11, replace=F)
    subset_smoker <- sample(metadata[[tissue]][metadata[[tissue]]$Smoking==2,"Donor"], 11, replace=F)
    subset <- metadata[[tissue]][metadata[[tissue]]$Donor %in% c(subset_never, subset_ex, subset_smoker),]
    saveRDS(subset, paste0("Downsampling/11/", tissue, "/metadata_downsampling_", iter, ".rds"))
    
    # if(sum(metadata[[tissue]]$Smoking==2)>=50 & sum(metadata[[tissue]]$Smoking==0)>=50){
    #   dir.create(paste0("Downsampling/50/", tissue, "/"), recursive = T)
    #   subset_never <- sample(metadata[[tissue]][metadata[[tissue]]$Smoking==0,"Donor"], 50, replace=F)
    #   subset_smoker <- sample(metadata[[tissue]][metadata[[tissue]]$Smoking==2,"Donor"], 50, replace=F)
    #   subset <- metadata[[tissue]][metadata[[tissue]]$Donor %in% c(subset_never, subset_ex, subset_smoker),]
    #   saveRDS(subset, paste0("Downsampling/50/", tissue, "/metadata_downsampling_", iter, ".rds"))
    # }
    # if(sum(metadata[[tissue]]$Smoking==2)>=80 & sum(metadata[[tissue]]$Smoking==0)>=80){
    #   dir.create(paste0("Downsampling/80/", tissue, "/"), recursive = T)
    #   subset_never <- sample(metadata[[tissue]][metadata[[tissue]]$Smoking==0,"Donor"], 80, replace=F)
    #   subset_smoker <- sample(metadata[[tissue]][metadata[[tissue]]$Smoking==2,"Donor"], 80, replace=F)
    #   subset <- metadata[[tissue]][metadata[[tissue]]$Donor %in% c(subset_never, subset_ex, subset_smoker),]
    #   saveRDS(subset, paste0("Downsampling/80/", tissue, "/metadata_downsampling_", iter, ".rds"))
    # }
    # if(sum(metadata[[tissue]]$Smoking==2)>=100 & sum(metadata[[tissue]]$Smoking==0)>=100){
    #   dir.create(paste0("Downsampling/100/", tissue, "/"), recursive = T)
    #   subset_never <- sample(metadata[[tissue]][metadata[[tissue]]$Smoking==0,"Donor"], 100, replace=F)
    #   subset_smoker <- sample(metadata[[tissue]][metadata[[tissue]]$Smoking==2,"Donor"], 100, replace=F)
    #   subset <- metadata[[tissue]][metadata[[tissue]]$Donor %in% c(subset_never, subset_ex, subset_smoker),]
    #   saveRDS(subset, paste0("Downsampling/100/", tissue, "/metadata_downsampling_", iter, ".rds"))
    # }
  }
}
