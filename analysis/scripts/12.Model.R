#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to run differential methylation analysis
# @software version: R=4.2.2

# Parsing
library(optparse)
parser <- OptionParser()
parser <- add_option(parser, opt_str=c("-t", "--tissue"), type="character",
                     dest="tissue",
                     help="Tissue")
parser <- add_option(parser, opt_str=c("-s", "--subset"), action="store_true", 
                     default=FALSE,
                     dest="subset",
                     help="In case we want to run a model with les samples")
#The following variables are only needed for our downsampling
parser <- add_option(parser, opt_str=c("-d", "--downsampling"), action="store_true", 
                     default=FALSE,
                     dest="downsampling",
                     help="Boolean specifying whether we want to model a downsampling set of samples")
parser <- add_option(parser, opt_str=c("-m", "--metadata"),  type="character",
                     dest="metadata",
                     help="Path for downsampled metadata")
options=parse_args(parser)
tissue=options$tissue
subset=options$subset
downsampling=options$downsampling
metadata=options$metadata
# tissue <- "Lung"

setwd(system("pwd", intern = T)) #If in linux
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
# setwd("..")


print("Reading data")
Sys.time()
data <- readRDS(paste0("tissues/", tissue, "/methylation_data.rds")) #From whole compressed data in 5.6G to compressed 1.4G/1.1Gb only in Lung (the highest number of samples)
Sys.time() #12 minutes to load 15 Gb
beta <- data

print("Reading metadata")
if(subset){
  metadata <- readRDS(paste0("tissues/", tissue, "/methylation_metadata_subset.rds"))
  beta <- beta[,colnames(beta) %in% metadata$SUBJID]
}else if(downsampling){
  print(metadata)
  number <- strsplit(last(strsplit(metadata, "_")[[1]]), "[.]")[[1]][1]
  output_file <- gsub(".rds", "_results.rds", metadata)
  metadata <- readRDS(metadata)
}else{
  metadata <- readRDS(paste0("tissues/", tissue, "/methylation_metadata.rds"))
}
metadata$Smoking <- as.factor(metadata$Smoking)
metadata$Ancestry <- as.factor(metadata$Ancestry)
if(sum(metadata$Ancestry=="AMR")<5){
  print("Not enough AMRs")
  metadata <- metadata[metadata$Ancestry!="AMR",]
  metadata$Ancestry <- droplevels(metadata$Ancestry)
  beta <- beta[,colnames(beta) %in% metadata$SUBJID]
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
if(sum(metadata$Ancestry=="AFR")<5){ #Only in testis
  metadata <- metadata[, !colnames(metadata) %in% c("Ancestry")]
  individual_variables <- individual_variables[!individual_variables %in% "Ancestry"]
}

#Add clinical traits to the metadata
expression_metadata <- readRDS(paste0("tissues/", tissue, "/metadata.rds"))
clinical_traits <- colnames(expression_metadata)[!colnames(expression_metadata) %in% c("Donor", "Sample", "HardyScale", "IschemicTime",
                                     "RIN", "ExonicRate", "PEER1", "PEER2",
                                     "Age", "Ancestry", "Sex", "BMI", "Smoking")]
if(length(clinical_traits)==0){
  print("no clinical traits added")
} else{
  clinical_data <- read.csv("data/public/histological_data.csv")
  if(tissue=="BreastMammaryTissue"){
    tissue_id <- "Breast - Mammary Tissue"
  } else if(tissue=="ColonTransverse"){
    tissue_id <- "Colon - Transverse"
  } else if(tissue=="KidneyCortex"){
    tissue_id <- "Kidney - Cortex"
  } else if(tissue=="MuscleSkeletal"){
    tissue_id <- "Muscle - Skeletal"
  } else{
    tissue_id <- tissue
  }
  clinical_data <- clinical_data[clinical_data$Tissue==tissue_id,]
  clinical_data <- clinical_data[clinical_data$Subject.ID %in% rownames(metadata),]
  clinical_data <- clinical_data[,colnames(clinical_data) %in% c("Subject.ID", clinical_traits)]
  metadata$Subject.ID <- rownames(metadata)
  metadata <- merge(metadata, clinical_data, by = "Subject.ID")
}

probes <- rownames(beta)

metadata_2 <- metadata[,c("PEER1", "PEER2", individual_variables, clinical_traits)]
names(metadata_2)
#I will use limma to create my contrasts, as other build in functions doesn't allow me to modify it
library(limma)

limma_function <- function(fit, x){
  covariate <<- x #makeContrast does not read the function's environment, so I add covariate to the general environment in my session
  contrast.matrix <- suppressWarnings(makeContrasts(covariate, levels=fit$design)) #Warnings due to change of name from (Intercept) to Intercept
  fitConstrasts <- suppressWarnings(contrasts.fit(fit, contrast.matrix)) #Warning due to Intercept name
  eb = eBayes(fitConstrasts)
  tt.smart.sv <- topTable(eb,adjust.method = "BH",number=Inf)
  return(tt.smart.sv)
}

beta <- sapply(beta, as.numeric)
rownames(beta) <- probes
M <- log2(beta/(1-beta)) # M is better for differential although beta should be plotted for interpretation

to_run_residuals <- colnames(metadata_2)[!colnames(metadata_2) %in% c("Smoking")]
mod_no_smoking <- model.matrix( as.formula(paste0("~", paste0(to_run_residuals, collapse="+"))), data =  metadata_2)
mod_2 <- model.matrix( as.formula(paste0("~", paste0(colnames(metadata_2), collapse="+"))), data =  metadata_2)

model_function <- function(mod, mod_no_smoking){
  print("Modelling")
  Sys.time()
  fit_M <- lmFit(M, mod)
  summary(decideTests(fit_M))
  # fit_b <- lmFit(beta, mod)
  
  if(!subset){
    #Saving residuals:
    fit_residuals <- lmFit(M, mod_no_smoking)
    residuals <- resid(fit_residuals, M)
    saveRDS(residuals, paste0("tissues/", tissue, "/methylation_residuals.rds")) #Correcting for covariates and demographic traits for the subset with smoking annotation but without smoking
  }
  
  if(!"Ancestry" %in% names(metadata)){ #Only happens in Testis
    to_run <- c("Smoking1", "Smoking2", "SEX2", "AGE", "BMI","Smoking1-Smoking2")
  } else if(length(levels(metadata$Ancestry))>2){ #Do we have Ancestry AMR?
    to_run <- c("Smoking1", "Smoking2", "AncestryAMR", "AncestryEUR", "SEX2", "AGE", "BMI", "AncestryAMR-AncestryEUR", "Smoking1-Smoking2")
  } else{
    to_run <- c("Smoking1", "Smoking2", "AncestryEUR", "SEX2", "AGE", "BMI", "Smoking1-Smoking2")
  }
  if(!"SEX" %in% colnames(metadata)){
    to_run <- to_run[!to_run=="SEX2"]
  }
  # to_run <- c("Smoking2") #For the subset
  print(to_run)
  res_M <- lapply(to_run, function(x) limma_function(fit_M, x) )
  names(res_M) <- to_run

  to_return <- res_M
  Sys.time()
  return(to_return)
}

res_2 <- model_function(mod_2, mod_no_smoking)

if(subset){
  saveRDS(res_2, paste0("tissues/", tissue, "/DML_results_subset.rds"))
}else if(downsampling){
  saveRDS(res_2, output_file)
}else{
  saveRDS(res_2, paste0("tissues/", tissue, "/DML_results.rds"))
}

print("Using 2 PEERs:")
print(paste0("     Using M: ", sum(res_2$Smoking2$adj.P.Val<0.05)))

