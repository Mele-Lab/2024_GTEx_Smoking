#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to correlate smoking with the PEER values
# @software version: R=4.2.2

setwd(system("pwd", intern = T)) #If in linux
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")

tissue <- "Lung"
print("Reading metadata")
metadata <- readRDS(paste0("tissues/", tissue, "/methylation_metadata.rds"))
metadata$Smoking <- as.factor(metadata$Smoking)
metadata$Ancestry <- as.factor(metadata$Ancestry)
if(sum(metadata$Ancestry=="AMR")<5){
  print("Not enough AMRs")
  metadata <- metadata[metadata$Ancestry!="AMR",]
  metadata$Ancestry <- droplevels(metadata$Ancestry)
  data <- data[,colnames(data) %in% rownames(metadata)]
}
metadata$SEX <- as.factor(metadata$SEX)
if(length(levels(metadata$SEX))==1){
  print("Sexual tissue")
  metadata <- metadata[,-which(names(metadata) == "SEX")]
  individual_variables <- c("Smoking", "Ancestry", "AGE", "BMI", "TRISCHD", "DTHHRDY")
} else{
  individual_variables <- c("Smoking", "Ancestry", "SEX", "AGE", "BMI", "TRISCHD", "DTHHRDY")
}
metadata$DTHHRDY <- as.factor(metadata$DTHHRDY)
rownames(metadata) <- metadata$SUBJID
metadata$SUBJID <- NULL


#Check correlations 

correlations <- c()
p_values <- c()
for(i in 1:20){ #There are 20 peer values
  peer <- paste0("PEER", i)
  t <- cor.test(as.numeric(metadata$Smoking), metadata[,eval(peer)], method="pearson")
  correlations <- c(correlations, t$estimate)
  p_values <- c(p_values, t$p.value)
}
p_adjusted <- p.adjust(p_values)

names(correlations) <- paste0("PEER", 1:20)
names(p_adjusted) <- paste0("PEER", 1:20)


#We want to put black boxes in the significant correlations
columns_to_highlight <- which(p_adjusted<=0.05)
btm <- 1 - (columns_to_highlight / length(correlations))
top <- btm + (1/length(correlations))
to_save <- list("correlations"=correlations, "columns_to_highlight"=columns_to_highlight, "top"=top, "btm"=btm)
saveRDS(to_save, file = "../figures/data/PEER_correlations.rds")

#Plot is in another code
