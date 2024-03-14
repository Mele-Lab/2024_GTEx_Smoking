#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez and Raquel Garcia-Perez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to run residuals on gene expression
# @software version: R=4.2.2
rm(list=ls())

#Set path 
# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")


metadata <- readRDS("output/metadata.rds")

table <- as.data.frame(cbind(names(metadata), "NA"))
names(table) <- c("Tissue", "Diseases")

for(tissue in names(metadata)){
  diseases <- names(metadata[[tissue]])[!names(metadata[[tissue]]) %in% c("Donor", "Sample", "HardyScale", "IschemicTime", "RIN", "ExonicRate", "PEER1", "PEER2", "Age", "Ancestry", "Sex", "BMI", "Smoking")]
  table[table$Tissue==tissue,"Diseases"] <- paste0(diseases, collapse=", ")
}

tissue_names <- read.csv("data/public/tissue_abreviation.txt")
table$Tissue <- sapply(table$Tissue, function(tissue) tissue_names$Name[tissue_names$X==tissue])

library(xlsx)
write.xlsx(table, 'output/Supplementary_table_1.xlsx', row.names = F)
