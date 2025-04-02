#!/usr/bin/env Rscript
# @Author: Rogério Eduardo Ramos Ribeiro
# @E-mail: rogerio.e.ramos.ribeiro@gmail.com
# @Description: Code to rerun Horvath and Vadim models on the epigenic dataset 
# @software version: R=4.2.2

# setwd(system("pwd", intern = T)) #If in linux
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio

## Load libraries
suppressPackageStartupMessages(library(tidyverse))
library(readr)
library(lmtest)
library(data.table)
library(impute)
library(ggpubr)

# Function to remove everything after the second "-"
remove_after_second_dash <- function(x) {
  parts <- unlist(strsplit(x, "-"))
  if (length(parts) > 2) {
    return(paste(parts[1:2], collapse = "-"))
  } else {
    return(x)
}

# Load normalization functions (Horvarth clock only)
source("19.helper_normalization.R")

# Horvath ----
# Functions
trafo= function(x,adult.age=20) { x=(x+1)/(1+adult.age); y=ifelse(x<=1, log( x),x-1);y }
anti.trafo= function(x,adult.age=20) { ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }
probeAnnotation21kdatMethUsed=read.csv("public/clockData/probeAnnotation21kdatMethUsed.csv")
probeAnnotation27k=read.csv("public/clockData/datMiniAnnotation27k.csv")
datClock=read.csv("public/clockData/AdditionalFile3.csv")

# Set tissues to run analysis

tissues <- c("Lung", "ColonTransverse", "BreastMammaryTissue", "MuscleSkeletal", "Ovary", "Prostate", "Testis", "WholeBlood")

sex_tissues <- c("Ovary", "Prostate", "Testis")

## Load methylation and metadata ----
predicted_age_horvath_list <- list()
for (tissue in tissues) {
  # Load data
  data <- fread(paste0("data_mn5/tissues/", tissue, "/methylation_data.csv"))
  colnames(data) <- unname(sapply(colnames(data), remove_after_second_dash))
  
  metadata <- readRDS(paste0("data_mn5/tissues/", tissue, "/methylation_metadata_subset.rds"))
  
  # Filter for samples with age and smoking data
  data <- data %>% 
    select("probe", all_of(metadata %>% pull("SUBJID")))
  
  # Keep only probes in the 21K
  data21k <- data %>%
    filter(probe %in% probeAnnotation21kdatMethUsed$Name)
  
  # Impute missing genes in the 21K probe array with NA values
  match1 <- match(probeAnnotation21kdatMethUsed$Name, data21k$probe)
  missing_probes <- probeAnnotation21kdatMethUsed$Name[is.na(match1)]
  
  missing_probes_data <- matrix(data = NA, nrow = length(missing_probes), ncol = ncol(data21k))
  missing_probes_data[,1] <- missing_probes
  colnames(missing_probes_data) <- colnames(data21k)
  
  missing_probes_data <- as.data.frame(missing_probes_data)
  data21kImputed <- rbind(data21k, missing_probes_data) %>%
    as_tibble() %>%
    mutate_at(2:ncol(data21k), as.numeric)
  
  ## Run pipeline ----
  
  ### Step 1: Impute missing data
  datMethUsed <- t(data21kImputed[,-1])
  colnames(datMethUsed) <- as.character(data21kImputed$probe)
  
  dimnames1 <- dimnames(datMethUsed)
  datMethUsed <- data.frame(t(impute.knn(t(datMethUsed))$data))
  dimnames(datMethUsed) <- dimnames1
  
  data21kImputed <- datMethUsed
  
  ### Step 2: Normalize data
  data21kImputed <- BMIQcalibration(
    datM = data21kImputed,
    goldstandard.beta = probeAnnotation21kdatMethUsed$goldstandard2,
    plots = FALSE
  )
  
  ### Step 3: Predict age with Horvath clock
  selectCpGsClock <- is.element(dimnames(data21kImputed)[[2]], as.character(datClock$CpGmarker[-1]))
  
  datMethClock0 <- data.frame(data21kImputed[, selectCpGsClock])
  datMethClock <- data.frame(datMethClock0[as.character(datClock$CpGmarker[-1])])
  
  predictedAge <- as.numeric(anti.trafo(datClock$CoefficientTraining[1] + as.matrix(datMethClock) %*% as.numeric(datClock$CoefficientTraining[-1])))
  
  predicted_age_horvath_list[[tissue]] <- data.frame("SUBJECT" = row.names(datMethClock0), 
                                                     "meth_age" = predictedAge)
  
  rm(data)
}


saveRDS(predicted_age_horvath_list, "../figures/data/predicted_age_ht.rds")

### Compare Smokers vs never smokers----
# Are Smokers predicted with an higher error compared with never-smokers? 

# Are Smokers with a higher error compared with never-smokers, if we adjust for 
residuals_regression <- data.frame("tissue" = c(), "age" = c(), "residual" = c(), "smoking_status" = c())
horvath_coef_table <- c()
variables <- c("Ancestry", "SEX", "AGE", "BMI", "TRISCHD", "DTHHRDY",
               "Smoking")

for (tissue in tissues){
  ages_pred <- predicted_age_horvath_list[[tissue]]
  ages_real <- readRDS(paste0("data_mn5/tissues/", tissue, "/methylation_metadata_subset.rds"))
  
  ages <- merge(ages_pred, ages_real, by.x = "SUBJECT" , by.y = "SUBJID")
  ages$error <- ages$meth_age - ages$AGE
  
  if (tissue %in% sex_tissues){
    variables_2 <- variables
    ages <- ages %>% 
      select(all_of(variables_2), error) %>%
    mutate_at(.vars = c("Ancestry", "DTHHRDY", "Smoking"), as.factor)
    
  }else{
    variables_2 <- variables
    ages <- ages %>% 
      select(all_of(variables_2), error) %>%
      mutate_at(.vars = c("SEX", "Ancestry", "DTHHRDY", "Smoking"), as.factor)
    
  }
  
  # Regression to run differences
  formula <- as.formula(paste("error ~", paste(variables_2, collapse = " + ")))
  # Regression to get residuals
  variables_resids <- variables_2[variables_2 != "Smoking"]
  formula_resids <- as.formula(paste("error ~", paste(variables_resids, collapse = " + ")))
  
  regression_horvath <- lm(formula, data = ages)
  
  
  horvath_coef_table <- rbind(horvath_coef_table, 
                              c(tissue, "Never vs Smoking",
                                summary(regression_horvath)$coefficients["Smoking2","Pr(>|t|)"], 
                                summary(regression_horvath)$coefficients["Smoking2","Estimate"],
                                confint(regression_horvath)["Smoking2", "2.5 %"],
                                confint(regression_horvath)["Smoking2", "97.5 %"]
                                )
  )
  horvath_coef_table <- rbind(horvath_coef_table, 
                              c(tissue, "Never vs Ex-Smoker",
                                summary(regression_horvath)$coefficients["Smoking1","Pr(>|t|)"],
                                summary(regression_horvath)$coefficients["Smoking1","Estimate"],
                                confint(regression_horvath)["Smoking1", "2.5 %"],
                                confint(regression_horvath)["Smoking1", "97.5 %"]
                                )
  )
  regression_horvath <- lm(formula, data = ages)
  
  #Model comparing Smokers and ex-smokers
  ages$Smoking <- factor(ages$Smoking, levels = c("1", "0", "2"))
  regression_horvath <- lm(formula, data = ages)
  horvath_coef_table <- rbind(horvath_coef_table, 
                              c(tissue, "Ex-Smoker vs Smoker",
                                summary(regression_horvath)$coefficients["Smoking2","Pr(>|t|)"],
                                summary(regression_horvath)$coefficients["Smoking2","Estimate"],
                                confint(regression_horvath)["Smoking2", "2.5 %"],
                                confint(regression_horvath)["Smoking2", "97.5 %"]
                                )
  )
  
  
  
  # Save residuals from the regression to plot
  regression_horvath_2 <- lm(formula_resids, data = ages)
  residuals <- resid(regression_horvath_2)
  residuals_regression <- rbind(residuals_regression, 
                                data.frame("tissue" = tissue, 
                                           "age" = ages$AGE,
                                           "error" = residuals,
                                           "smoking_status" = ages$Smoking)
  )
}


colnames(horvath_coef_table) <- c("tissue", "comparison", "p.value", "estimate", "lCI", "hCI")
horvath_coef_table <- horvath_coef_table %>%
  as.data.frame() %>% 
  mutate(p.value = as.numeric(p.value)) %>% 
  mutate(FDR = p.adjust(p.value, method = "fdr"))


data_to_export <- list()
data_to_export$resids <- residuals_regression
data_to_export$comparison <- horvath_coef_table

saveRDS(data_to_export, "../figures/data/hovarth_molecular_clocks.rds")



# Vadim clocks ----
## Load molecular clocks coefficients for Vadim clocks----
trafo= function(x,adult.age=20) { x=(x+1)/(1+adult.age); y=ifelse(x<=1, log( x),x-1);y }
anti.trafo= function(x,adult.age=20) { ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }
probeAnnotation21kdatMethUsed=read.csv("public/clockData/probeAnnotation21kdatMethUsed.csv")
probeAnnotation27k=read.csv("public/clockData/datMiniAnnotation27k.csv")
datClock=read.csv("public/clockData/AdditionalFile3.csv")


clocks <- list()
# CauAge has too many missing samples
#clocks$CauAge <- read.csv("public/clockData/YingCausAge.csv")
clocks$DamAge <- read.csv("public/clockData/YingDamAge.csv")
clocks$AdaptAge <- read.csv("public/clockData/YingAdaptAge.csv")

# Check for missing data across molecular clocks ----
#for (clock in names(clocks)){
#  cat("Clock ", clock, "has ", sum(!clocks[[clock]]$term[-1] %in% data[["Lung"]]$probe), 
#      "missing probes from ", nrow(clocks[[clock]]), "\n")
#}

clocks_to_run <- clocks

# Keep data only from the 450K array (same array used to train the Vadim Clocks)
array450k <- read.csv("public/clockData/humanmethylation450_15017482_v1-2.csv")

predicted_ages_vadim <- list()

for (tissue in tissues){
  # Load data
  data <- fread(paste0("data_mn5/tissues/", tissue, "/methylation_data.csv"))
  colnames(data) <- unname(sapply(colnames(data), remove_after_second_dash))
  metadata <- readRDS(paste0("data_mn5/tissues/", tissue, "/methylation_metadata_subset.rds"))
  
  
  data <- data %>% 
    select("probe", all_of(metadata %>% pull("SUBJID")))

  
  data450k <- data %>%  filter(probe %in% array450k$Name)
  
  
  # Run pipeline ----
  # Due to missing genes from the 450K probe array, we need to input them.
  # Imputation is based on average methylation values across the samples 
  # This is based on Horvath methodology to deal with missing data in the online cal
  # Due to computational limits, I only used the the probes in the model as background 
  # Using the 450 could potencial improve something, but I not sure the best way to 
  # do it is. Horvath uses all the probes in the array (21k)
  
  ## Adapt Age----
  
  # Add missing probes with value "NA"
  missing_probes <- clocks_to_run$AdaptAge$term[!clocks_to_run$AdaptAge$term %in% data450k$probe]
  missing_probes <- missing_probes[-1] # removing Intercept
  
  
  tissue_data <- data450k
  
  missing_probes_data <- matrix(data = NA, nrow = length(missing_probes), ncol = ncol(tissue_data))
  missing_probes_data[,1] <- missing_probes
  colnames(missing_probes_data) <-  colnames(tissue_data)
  
  missing_probes_data <- as.data.frame(missing_probes_data)
  
  # Filter data (keeping only the clock features)
  #tissue_data_cpg_filtered <- tissue_data %>% 
  #  filter(probe %in% clocks_to_run$AdaptAge$term[-1])
  
  # Add NA data
  tissue_data <- rbind(tissue_data, missing_probes_data)
  tissue_data <- tissue_data %>% 
    as_tibble() %>%
    mutate_at(2:ncol(tissue_data), as.numeric)
  
  datMethUsed <- t(tissue_data[,-1])
  colnames(datMethUsed)=as.character(tissue_data$probe)
  
  #Note There another faster imputation method in the Horvath source code
  # but this is fine
  
  dimnames1=dimnames(datMethUsed)
  datMethUsed= data.frame(t(impute.knn(t(datMethUsed))$data))
  dimnames(datMethUsed)=dimnames1
  
  # Predict age in these tissues
  AdaptAge <- clocks_to_run$AdaptAge
  
  datMethClock <-  data.frame(datMethUsed[as.character(AdaptAge$term[-1])])
  predictedAge <- as.numeric(
    AdaptAge$estimate[1]+as.matrix(datMethClock)%*% as.numeric(AdaptAge$estimate[-1]))
  
  predictedAdaptAge <- data.frame("SUBJECT" =row.names(datMethClock), 
                               "Adapt_age" = predictedAge)

  ## DamAge Age----
  
  # Add missing probes with value "NA"
  missing_probes <- clocks_to_run$DamAge$term[!clocks_to_run$DamAge$term %in% data450k]
  missing_probes <- missing_probes[-1] # removing Intercept
  
  tissue_data <- data450k
  
  missing_probes_data <- matrix(data = NA, nrow = length(missing_probes), ncol = ncol(tissue_data))
  missing_probes_data[,1] <- missing_probes
  colnames(missing_probes_data) <-  colnames(tissue_data)
  
  missing_probes_data <- as.data.frame(missing_probes_data)
  
  # Filter data (keeping only the clock features)
  tissue_data_cpg_filtered <- tissue_data %>% 
    filter(probe %in% clocks_to_run$DamAge$term[-1])
  
  # Add NA data
  tissue_data <- rbind(tissue_data_cpg_filtered, missing_probes_data)
  tissue_data <- tissue_data %>% 
    as_tibble() %>%
    mutate_at(2:ncol(tissue_data), as.numeric)
  
  datMethUsed <- t(tissue_data[,-1])
  colnames(datMethUsed)=as.character(tissue_data$probe)
  
  #Note There another faster imputation method in the Horvath source code
  # but this is fine
  
  dimnames1=dimnames(datMethUsed)
  datMethUsed= data.frame(t(impute.knn(t(datMethUsed))$data))
  dimnames(datMethUsed)=dimnames1
  
  # Predict age in these tissues
  DamAge <- clocks_to_run$DamAge
  
  datMethClock <-  data.frame(datMethUsed[as.character(DamAge$term[-1])])
  predictedAge <- as.numeric(
    DamAge$estimate[1]+as.matrix(datMethClock)%*% as.numeric(DamAge$estimate[-1]))
  
  predictedDamAGE <- data.frame("SUBJECT" =row.names(datMethClock), 
                               "Dam_age" = predictedAge)
  
  
  predicted_age <- merge(predictedDamAGE, predictedAdaptAge, by ="SUBJECT") 
  predicted_ages_vadim[[tissue]] <- predicted_age  
}

saveRDS(predicted_ages_vadim, "../figures/data/vadim_age_pred.rds")

### Comparison between non-smoker and smokers ----
#Are Smokers predicted with an higher error compared with never-smokers? 
#Do smokers have a higher error compared with never-smokers, if we adjust for 
# Other factors? 
variables <- c("Ancestry", "SEX", "AGE", "BMI", "TRISCHD", "DTHHRDY",
               "Smoking")

AdaptAge_coef_table <- c()
DamAge_coef_table <- c()

residuals_regression_vadim <- data.frame()

for (tissue in tissues){
  ages_pred <- predicted_ages_vadim[[tissue]]
  ages_real <- readRDS(paste0("data_mn5/tissues/", tissue, "/methylation_metadata_subset.rds"))
  
  ages <- merge(ages_pred, ages_real, by.x = "SUBJECT" , by.y = "SUBJID")
  
  ages$errorAdaptAge <- ages$Adapt_age - ages$AGE
  ages$errorDamAge <- ages$Dam_age - ages$AGE
  
  if (tissue %in% sex_tissues){
    variables_2 <- variables
    ages <- ages %>% 
      select(all_of(variables_2), errorAdaptAge, errorDamAge) %>%
      mutate_at(.vars = c("Ancestry", "DTHHRDY", "Smoking"), as.factor)
    
  }else{
    variables_2 <- variables
    ages <- ages %>% 
      select(all_of(variables_2),  errorAdaptAge, errorDamAge) %>%
      mutate_at(.vars = c("SEX", "Ancestry", "DTHHRDY", "Smoking"), as.factor)
  }
  
  # For AdaptAge
  # Regression to run differences
  formula <- as.formula(paste("errorAdaptAge ~", paste(variables_2, collapse = " + ")))
  # Regression to get residuals
  variables_resids <- variables_2[variables_2 != "Smoking"]
  formula_resids <- as.formula(paste("errorAdaptAge ~", paste(variables_resids, collapse = " + ")))
  
  regression_adapt <- lm(formula, data = ages)
  
  
  AdaptAge_coef_table <- rbind(AdaptAge_coef_table, 
                              c(tissue, "Never vs Smoking",
                                summary(regression_adapt)$coefficients["Smoking2","Estimate"],
                                summary(regression_adapt)$coefficients["Smoking2","Pr(>|t|)"],
                                confint(regression_adapt)["Smoking2", "2.5 %"],
                                confint(regression_adapt)["Smoking2", "97.5 %"]
                                )
  )
  AdaptAge_coef_table <- rbind(AdaptAge_coef_table, 
                              c(tissue, "Never vs Ex-Smoker",
                                summary(regression_adapt)$coefficients["Smoking1","Estimate"],
                                summary(regression_adapt)$coefficients["Smoking1","Pr(>|t|)"], 
                                confint(regression_adapt)["Smoking1", "2.5 %"],
                                confint(regression_adapt)["Smoking1", "97.5 %"])
  )
  
  #Model comparing Smokers and ex-smokers
  ages$Smoking <- factor(ages$Smoking, levels = c("1", "0", "2"))
  regression_adapt <- lm(formula, data = ages)
  AdaptAge_coef_table <- rbind(AdaptAge_coef_table, 
                              c(tissue, "Ex-Smoker vs Smoker",
                                summary(regression_adapt)$coefficients["Smoking2","Estimate"],
                                summary(regression_adapt)$coefficients["Smoking2","Pr(>|t|)"],
                                confint(regression_adapt)["Smoking2", "2.5 %"],
                                confint(regression_adapt)["Smoking2", "97.5 %"]
                              ))
  

  # Save residuals from the regression to plot
  regression_adaptAge <- lm(formula_resids, data = ages)
  residuals <- resid(regression_adaptAge)
  residuals_regression_vadim <- rbind(residuals_regression_vadim, 
                                data.frame("tissue" = tissue, 
                                           "age" = ages$AGE,
                                           "error" = residuals,
                                           "smoking_status" = ages$Smoking, 
                                           "clock" = "AdaptAge")
  )
  
  # For DamAge
  ages$Smoking <- factor(ages$Smoking, levels = c("0", "1", "2"))
  
  # Regression to run differences
  formula <- as.formula(paste("errorDamAge ~", paste(variables_2, collapse = " + ")))
  # Regression to get residuals
  variables_resids <- variables_2[variables_2 != "Smoking"]
  formula_resids <- as.formula(paste("errorDamAge ~", paste(variables_resids, collapse = " + ")))
  
  regression_dam <- lm(formula, data = ages)
  
  DamAge_coef_table <- rbind(DamAge_coef_table, 
                               c(tissue, "Never vs Smoking",
                                 summary(regression_dam)$coefficients["Smoking2","Estimate"],
                                 summary(regression_dam)$coefficients["Smoking2","Pr(>|t|)"],
                                 confint(regression_dam)["Smoking2", "2.5 %"],
                                 confint(regression_dam)["Smoking2", "97.5 %"])
  )
  DamAge_coef_table <- rbind(DamAge_coef_table, 
                               c(tissue, "Never vs Ex-Smoker",
                                 summary(regression_dam)$coefficients["Smoking1","Estimate"],
                                 summary(regression_dam)$coefficients["Smoking1","Pr(>|t|)"],
                                 confint(regression_dam)["Smoking1", "2.5 %"],
                                 confint(regression_dam)["Smoking1", "97.5 %"])
  )
  
  #Model comparing Smokers and ex-smokers
  ages$Smoking <- factor(ages$Smoking, levels = c("1", "0", "2"))
  regression_dam <- lm(formula, data = ages)
  DamAge_coef_table <- rbind(DamAge_coef_table, 
                               c(tissue, "Ex-Smoker vs Smoker",
                                 summary(regression_dam)$coefficients["Smoking2","Estimate"],
                                 summary(regression_dam)$coefficients["Smoking2","Pr(>|t|)"],
                                 confint(regression_dam)["Smoking2", "2.5 %"],
                                 confint(regression_dam)["Smoking2", "97.5 %"])
                                 )
  
  
  # Save residuals from the regression to plot
  regression_DamAge <- lm(formula_resids, data = ages)
  residuals <- resid(regression_DamAge)
  residuals_regression_vadim <- rbind(residuals_regression_vadim, 
                                      data.frame("tissue" = tissue, 
                                                 "age" = ages$AGE,
                                                 "error" = residuals,
                                                 "smoking_status" = ages$Smoking, 
                                                 "clock" = "DamAge_coef_table")
  )
}


colnames(DamAge_coef_table) <- c("tissue", "comparison", "estimate", "p.value", "lCI", "hCI")
colnames(AdaptAge_coef_table) <- c("tissue", "comparison","estimate", "p.value", "lCI", "hCI")

DamAge_coef_table <- DamAge_coef_table %>%
  as.data.frame() %>% 
  mutate(p.value = as.numeric(p.value)) %>% 
  mutate(FDR = p.adjust(p.value, method = "fdr"))

AdaptAge_coef_table <- AdaptAge_coef_table %>%
  as.data.frame() %>% 
  mutate(p.value = as.numeric(p.value)) %>% 
  mutate(FDR = p.adjust(p.value, method = "fdr"))

data_to_export <- list()
data_to_export$adapt <- list()
data_to_export$damn <- list()



data_to_export$resids <- residuals_regression_vadim
data_to_export$adapt$comparison <- AdaptAge_coef_table
data_to_export$damn$comparison <- DamAge_coef_table

saveRDS(data_to_export, "../figures/data/figure_clock_vadim.rds")
