#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to run residuals on gene expression 
# @software version: R=4.2.2
rm(list=ls())

# Load libraries ####
suppressMessages(library(edgeR))

#Set path 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")

# -------------- #
print(Sys.time())
#-------------- #

tissues <- list.dirs("tissues/", full.names = F)[-1]


sexual_tissues <- c("Prostate", "Testis", "Vagina", "Uterus", "Ovary")
#Do a for loop 
for(tissue in tissues){
  print(tissue)
  metadata <- readRDS(paste0("tissues/", tissue, "/metadata.rds"))
  if(tissue %in% sexual_tissues){
    metadata <- metadata[,names(metadata)!="Sex"]
  }

  #Exclude AMR if not enough samples
  if(sum(metadata$Ancestry=="AMR")==0){
    metadata$Ancestry <- droplevels(metadata$Ancestry)
  }
  
  # Counts
  counts <- readRDS(paste0("tissues/", tissue, "/counts.rds")) #Counts of expressed genes in samples of interest for the given tissue
  
  #Make sure that the sample ids in metadata are for the samples with expression
  metadata <- metadata[metadata$Sample %in% colnames(counts),] #Only needed in lung, because I previously edited the metadata for some test
  metadata <- metadata[!duplicated(metadata$Sample),] 
  
  # Create DGEList object
  dge <- DGEList(counts)
  
  # Calculate normalization factors (does not do the normalization, only computes the factors)
  dge <- calcNormFactors(dge)
  
  # Voom
  v <- voom(dge, design = NULL, normalize="quantile", save.plot=F, plot = F) # samples are treated as replicates
  
  if(tissue %in% sexual_tissues){
    individual_traits <- c("Age", "Ancestry", "BMI")
  } else{
    individual_traits <- c("Age", "Ancestry", "Sex", "BMI")
  }
  covariates <- names(metadata)[!names(metadata) %in% c("Donor", "Sample", individual_traits)]
  covariates <- covariates[!covariates %in% "Smoking"]
  
  #Including demographic traits for other purposes in other scripts:
  covariates <- c(covariates, individual_traits) 
  fml_args_mod <- paste(c(covariates), collapse = " + ")
  mod_2 <- model.matrix( as.formula(paste(" ~  ", paste(fml_args_mod,collapse = " "))), data =  metadata)
  fit_2 <- lmFit(v, mod_2)
  exprs_residuals_2 <- resid(fit_2, v)
  saveRDS(exprs_residuals_2, paste0("tissues/", tissue, "/expression_residuals_demographic_traits_no_smoking.rds"))
}


# -------------- #
print(Sys.time())
#-------------- #