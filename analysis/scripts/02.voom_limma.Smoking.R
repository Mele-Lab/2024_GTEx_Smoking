#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to run differential expression analysis with smoking
# @software version: R=4.2.2

# Load libraries ####
suppressMessages(library(edgeR)) #Already includes limma
library(dplyr) #function last
library(optparse)
options(warn=-1)

#Set path 
setwd(system("pwd", intern = T)) #If in linux
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
# setwd("..")

# Parsing
parser <- OptionParser()
parser <- add_option(parser, opt_str=c("-t", "--tissue"), type="character",
                     dest="tissue",
                     help="Tissue")
parser <- add_option(parser, opt_str=c("-i", "--interaction"), type="character",
                     dest="interaction", default = NULL,
                     help="Interaction to include in the model (e.g., Sex:Smoking)")
parser <- add_option(parser, opt_str=c("-s", "--subset"), action="store_true", 
                     default=FALSE,
                     dest="subset",
                     help="Boolean specifying whether we want to model a smaller set of samples")
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
interaction=options$interaction
subset=options$subset
downsampling=options$downsampling
metadata=options$metadata
n=options$n

threshold <- 20 #for interactions

# tissue <- "Thyroid"
# tissue <- "Lung"

# interaction <- NULL
# interaction <- "Age:Smoking BMI:Smoking Sex:Smoking Ancestry:Smoking"

print(paste("Running", tissue))

# Reading metadata
if(subset){
  metadata <- readRDS(paste0("tissues/", tissue, "/metadata_subset.rds"))
}else if(downsampling){
  print(metadata)
  number <- strsplit(last(strsplit(metadata, "_")[[1]]), "[.]")[[1]][1]
  output_file <- gsub(".rds", "_results.rds", metadata)
  metadata <- readRDS(metadata)
}else{
  metadata <- readRDS(paste0("tissues/", tissue, "/metadata.rds"))
}

#If sexual tissue, we should remove the covariate sex
sex_tissues <- c("Vagina", "Uterus", "Ovary", "Prostate", "Testis")
if(tissue %in% sex_tissues){ metadata <- metadata[,names(metadata)!="Sex"] }

#Maybe in a tissue we do not have all potential levels in a variable
if(length(unique(metadata$HardyScale))<5){metadata$HardyScale <- droplevels(metadata$HardyScale)}
if(length(unique(metadata$Ancestry))<3){metadata$Ancestry <- droplevels(metadata$Ancestry)}
if(!downsampling){ #Don't exclude samples
  if(sum(metadata$Ancestry=="AMR")<5){ #It only happens when we are running with a subset of the samples
    print("Not enough AMRs")
    metadata <- metadata[metadata$Ancestry!="AMR",]
    metadata$Ancestry <- droplevels(metadata$Ancestry)
  }
}


#Variables we are interested in:
individual_traits <- c("Age", "Ancestry", "Sex", "BMI", "Smoking")
if(tissue %in% sex_tissues){ individual_traits <- individual_traits[individual_traits!="Sex"] }

#Vector of variables we will correct for:
covariates <- c(colnames(metadata)[!colnames(metadata) %in% c("Donor", "Sample", individual_traits)])


if(!is.null(interaction)){
  interaction <- strsplit(interaction, " ")[[1]]
  #Exclude interactions when we don't have enough samples:
  metadata_subset <- metadata[metadata$Smoking!=1,] #We are not interested in ex-smokers for interactions
  metadata_subset$Smoking <- droplevels(metadata_subset$Smoking)
  metadata_subset$Age_bucket <- cut(metadata_subset$Age,breaks = c(0,39,100),labels=c('under40','over40'))
  metadata_subset$BMI_bucket <- cut(metadata_subset$BMI,breaks = c(0,24.99,29.99,100),labels=c('healthy','overweight','obese'))
  
  acronym_interaction <- "Smoking"
  #Do we have enough samples per bucket in the interaction Disease:Age? 20 healthy younger than 40, 20 diseased younger than 40, 20 healthy above 40 and 20 diseased above 40 
  Age_boolean <- sum(table(metadata_subset[[acronym_interaction]],metadata_subset$Age_bucket) >= threshold)==length(levels(metadata_subset[[acronym_interaction]]))*length(levels(metadata_subset$Age_bucket)) 
  if(!Age_boolean){ #If FALSE, we don't test the interaction Age:Disease
    interaction <- interaction[!grepl("Age", interaction)]
  }
  if(sum(grepl("Sex", names(metadata_subset)))>0){ #Do we consider sex?
    Sex_boolean <- sum(table(metadata_subset[[acronym_interaction]],metadata_subset$Sex) >= threshold)==length(levels(metadata_subset[[acronym_interaction]]))*length(levels(metadata_subset$Sex)) #Do we have enough samples per bucket in the interaction Disease:Sex? 
    if(!Sex_boolean){ #If FALSE, we don't test the interaction Sex:Disease
      interaction <- interaction[!grepl("Sex", interaction)]
    }
  }else{
    interaction <- interaction[!grepl("Sex", interaction)]
  }
  BMI_boolean <- sum(table(metadata_subset[[acronym_interaction]],metadata_subset$BMI_bucket) >= threshold)==length(levels(metadata_subset[[acronym_interaction]]))*length(levels(metadata_subset$BMI_bucket)) #Do we have enough samples per bucket in the interaction Disease:Sex? 
  if(!BMI_boolean){ #If FALSE, we don't test the interaction BMI:Disease
    interaction <- interaction[!grepl("BMI", interaction)]
  }
  #For ancestry only interactions with AncestryEUR-AFR as we have very few AMR samples
  metadata_subset_1 <- metadata_subset[metadata_subset$Ancestry!="AMR",] 
  metadata_subset_1$Ancestry <- droplevels(metadata_subset_1$Ancestry)
  Ancestry_boolean <- sum(table(metadata_subset_1[[acronym_interaction]],metadata_subset_1$Ancestry) >= threshold)==length(levels(metadata_subset_1[[acronym_interaction]]))*length(levels(metadata_subset_1$Ancestry)) #Do we have enough samples per bucket in the interaction Disease:Sex? 
  if(!Ancestry_boolean){ #If Ancestry_boolean==FALSE, we don't test the interaction Ancestry:Smoking
    interaction <- interaction[!grepl("Ancestry", interaction)]
  } 
  print("Interactions to be included in the model:")
  print(interaction)
  if(length(interaction)==0){
    print("As no interactions are going to be included the code ends here")
    stop()
  }
}

# Functions ####
limma_lm2 <- function(fit, covariate, covariate_data){
   # print(covariate)
  if(covariate %in% nonEstimable(fit$design)){
    print("Covariate non estimable")
    print(covariate)
    return(0)
  }
  v.contrast <- rep(0,ncol(fit$design))
  if(grepl(":", covariate)){ #If interaction
    v.contrast[ which( colnames(fit$design) == covariate) ] <- 1
  } else if(grepl("-", covariate)){ #If comparison not in the default
    first <- strsplit(covariate, "-")[[1]][1]
    second <- strsplit(covariate, "-")[[1]][2]
    v.contrast[ which( colnames(fit$design) == first) ] <- 1
    v.contrast[ which( colnames(fit$design) == second) ] <- -1
  } else{
    v.contrast[ which( colnames(fit$design) == covariate) ] <- 1
  }  
  contrast.matrix <- cbind( "C1" = v.contrast)
  fitConstrasts <- contrasts.fit(fit,contrast.matrix)
  eb = eBayes(fitConstrasts)
  tt.smart.sv <- topTable(eb,adjust.method = "BH",number=Inf)
  return(tt.smart.sv)
}


# Reading expression data
counts <- readRDS(paste0("tissues/", tissue, "/counts.rds"))

#Change order of samples in count to match the order in metadata (relevant for voom limma) in case we use a different metadata (i.e., downsampling):
# identical(metadata$Sample, colnames(counts)) #In most cases this will already be the case
counts <- counts[, match(metadata$Sample, colnames(counts))] 

#  Normalising gene expression distributions:
# Create DGEList object
dge <- DGEList(counts)

# Calculate normalization factors (does not do the normalization yet, only computes the factors)
dge <- calcNormFactors(dge)

# Voom
v <- voom(dge, design = NULL, normalize="quantile", save.plot=F, plot = F) 


#Removing interactions with Sex for sexual tissues
if(tissue %in% sex_tissues){interaction <- interaction[!grepl("Sex", interaction)]}

#Preparing model
fml_args_mod <- paste(c(covariates, individual_traits, interaction), collapse = " + ")
print(fml_args_mod)

#Creating model
mod <- model.matrix( as.formula(paste(" ~  ", paste(fml_args_mod,collapse = " "))), data =  metadata)

# Limma fit
fit <- lmFit(v, mod)

#Limma test
if(!is.null(interaction)){
  to_run <- colnames(fit)[grepl(":", colnames(fit))]
  to_run <- to_run[grepl("Smoking2", to_run)]
  to_run <- to_run[!grepl("AMR", to_run)]
} else{
  individual_traits <- colnames(fit)[grepl("Age|Ancestry|BMI|Sex|Smoking", colnames(fit))]
  to_run <- c(individual_traits, "AncestryAMR-AncestryEUR", "Smoking1-Smoking2")
}

dea_res <- lapply(to_run, function(phenotype) limma_lm2(fit, phenotype, metadata ) )
names(dea_res) <- to_run

# Save table with results
if(!is.null(interaction)){
  saveRDS(dea_res, paste0("tissues/", tissue, "/interactions.rds"))
} else{
  if(subset){
    saveRDS(dea_res, paste0("tissues/", tissue, "/voom_limma_results_subset.rds"))
  }else if(downsampling){
    saveRDS(dea_res, output_file)
  }else{
    saveRDS(dea_res, paste0("tissues/", tissue, "/voom_limma_results.rds"))
  }
}

