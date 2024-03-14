#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to exclude clinical traits that do not retrieve differentially expressed genes
# @software version: R=4.2.2

#Set path 
setwd(system("pwd", intern = T)) #If in linux, as the code will be run from its parent folder rather than inside /scripts/ 
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
# setwd("..")

# Load libraries ####
suppressMessages(library(edgeR)) #Already includes limma
library(optparse)
options(warn=-1)

# Parsing
parser <- OptionParser()
parser <- add_option(parser, opt_str=c("-t", "--tissue"), type="character",
                     dest="tissue",
                     help="Tissue")

options=parse_args(parser)
tissue=options$tissue
# tissue <- "Lung"

print(tissue)

# Reading metadata
metadata <- readRDS(paste0("tissues/", tissue, "/metadata.rds"))
metadata_og <- metadata #metadata_og will be the object from which we will remove unnecessary covariates and save

metadata$Smoking <- NULL #We are not interested in this variable yet

#Vector of variables we will correct for:
covariates <- c("HardyScale", "IschemicTime", "RIN", "ExonicRate", "PEER1", "PEER2", "Age", "Ancestry", "Sex", "BMI")
#Getting clinical traits that we could potentially correct for:
diseases <- names(metadata)[!names(metadata) %in% c("Donor", "Sample", covariates, "Smoking")]

#If sexual tissue, we should remove the covariate sex
sex_tissues <- c("Vagina", "Uterus", "Ovary", "Prostate", "Testis")
if(tissue %in% sex_tissues){ 
  covariates <- covariates[covariates!="Sex"]
  metadata <- metadata[,names(metadata)!="Sex"] }

#Maybe in a tissue we do not have all potential levels in a variable
if(length(unique(metadata$HardyScale))<5){metadata$HardyScale <- droplevels(metadata$HardyScale)}
if(length(unique(metadata$Ancestry))<3){metadata$Ancestry <- droplevels(metadata$Ancestry)}


# Reading expression data
counts <- readRDS(paste0("tissues/",tissue,"/counts.rds"))

#  Normalising gene expression distributions:
# Create DGEList object
dge <- DGEList(counts)

# Calculate normalization factors (does not do the normalization yet, only computes the factors)
dge <- calcNormFactors(dge)

# Voom
v <- voom(dge, design = NULL, normalize="quantile", save.plot=F, plot = F) 

# Limma function ####
limma_lm <- function(fit, covariate, covariate_data){ #It returns the differentially expressed genes for a particular variable based on a limma fit object
  covariate <<- covariate  #makeContrast does not read the function's environment, so I add covariate to the general environment in my session
  contrast.matrix <- suppressWarnings(makeContrasts(covariate, levels=fit$design)) #Warnings due to change of name from (Intercept) to Intercept
  fitConstrasts <- suppressWarnings(contrasts.fit(fit,contrast.matrix))
  eb = eBayes(fitConstrasts)
  tt.smart.sv <- topTable(eb,adjust.method = "BH",number=Inf)
  return(tt.smart.sv)
}

#Run one model per disease:
for(disease in diseases){

  fml_args_mod <- paste(c(covariates, disease), collapse = " + ") #the formula we want to use

  #Creating model
  mod <- model.matrix( as.formula(paste(" ~  ", paste(fml_args_mod,collapse = " "))), data =  metadata) 
  
  # Limma fit
  fit <- lmFit(v, mod) 
  
  dea_res <- limma_lm(fit, paste0(disease, 1), metadata)
  
  #If we see at least one differentially expressed gene or 5?
  if(sum(dea_res$adj.P.Val<0.05)<=5){
    print(paste(disease, "is going to be excluded"))
    metadata_og <- metadata_og[,-which(colnames(metadata_og)==disease)]
  } else{
    print(paste(disease, "will be kept"))
  }
}

saveRDS(metadata_og, paste0("tissues/", tissue, "/metadata.rds")) #Update the metadata


