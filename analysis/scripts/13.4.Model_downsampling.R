#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to run differential methylation analysis in downsampling so we don't have to read too many times the heaby beta values
# @software version: R=4.2.2

# Parsing
library(optparse)
parser <- OptionParser()
parser <- add_option(parser, opt_str=c("-t", "--tissue"), type="character",
                     dest="tissue",
                     help="Tissue")

options=parse_args(parser)
tissue=options$tissue

# tissue <- "Lung"
# tissue <- "Prostate"

setwd(system("pwd", intern = T)) #If in linux
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
# setwd("..")


print("Reading data")
Sys.time()
data <- readRDS(paste0("tissues/", tissue, "/methylation_data.rds")) #From whole compressed data in 5.6G to compressed 1.4G/1.1Gb only in Lung (the highest number of samples)
Sys.time() #12 minutes to load 15 Gb
probes <- rownames(data)
data <- sapply(data, as.numeric)
rownames(data) <- probes
M <- log2(data/(1-data)) # M is better for differential although beta should be plotted for interpretation

for(number in c(25)){
# for(number in c(11, 25, 50, 35)){
    print(number)
  for(i in 1:50){
    print("Reading metadata")
    file <- paste0("Downsampling/", number, "/", tissue, "/metadata_downsampling_methylation_", i, ".rds")
    if(!file.exists(file)){break}
    metadata <- readRDS(file)
    output_file <- gsub(".rds", "_results.rds", file)
    
    metadata$Smoking <- as.factor(metadata$Smoking)
    metadata$Ancestry <- as.factor(metadata$Ancestry)
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

    #Add clinical traits to the metadata
    expression_metadata <- readRDS(paste0("tissues/", tissue, "/metadata.rds"))
    clinical_traits <- colnames(expression_metadata)[!colnames(expression_metadata) %in% c("Donor", "Sample", "HardyScale", "IschemicTime",
                                                                                           "RIN", "ExonicRate", "PEER1", "PEER2",
                                                                                           "Age", "Ancestry", "Sex", "BMI", "Smoking")]
    if(length(clinical_traits)==0){
      print("no clinical traits added")
      metadata$Subject.ID <- rownames(metadata)
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
    rownames(metadata) <- metadata$Subject.ID
    
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
    
    mod_2 <- model.matrix( as.formula(paste0("~", paste0(colnames(metadata_2), collapse="+"))), data =  metadata_2)
    M_subset <- M[,colnames(M) %in% rownames(mod_2)]
    M_subset <- M_subset[,match(colnames(M_subset), rownames(mod_2))]
    print("Modelling")
    Sys.time()
    fit_M <- lmFit(M_subset, mod_2)
    
    to_run <- c("Smoking2")
    print(to_run)
    res_M <- lapply(to_run, function(x) limma_function(fit_M, x) )
    names(res_M) <- to_run
    
    saveRDS(res_M, output_file)
    print(paste0("     Using M: ", sum(res_M$Smoking2$adj.P.Val<0.05)))
  }
}





