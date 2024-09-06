#!/usr/bin/env Rscript
# @Author: Rogério Eduardo Ramos Ribeiro
# @E-mail: rogerio.e.ramos.ribeiro@gmail.com
# @Description: Code to rerun Horvath and Vadim models on the epigenic dataset 
# @software version: R=4.2.2

# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio

## Load libraries
suppressPackageStartupMessages(library(tidyverse))
library(readr)
library(lmtest)
library(data.table)
library(impute)
library(ggpubr)

# Load normalization functions (Horvarth clock only)
source("19.helper_normalization.R")

# Horvath ----
# Functions
trafo= function(x,adult.age=20) { x=(x+1)/(1+adult.age); y=ifelse(x<=1, log( x),x-1);y }
anti.trafo= function(x,adult.age=20) { ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }
probeAnnotation21kdatMethUsed=read.csv("data/public/clockData/probeAnnotation21kdatMethUsed.csv")
probeAnnotation27k=read.csv("data/public/clockData/datMiniAnnotation27k.csv")
datClock=read.csv("data/public/clockData/AdditionalFile3.csv")

# Set tissues to run analysis

#tissues <- c("Lung", "ColonTransverse", "BreastMammaryTissue", "MuscleSkeletal", "Ovary", "Prostate", "Testis", "WholeBloods")
#sex_tissues <- c("Ovary", "Prostate", "Testis")
tissues <- "Lung"

## Load methylation and metadata ----
data <- list()
metadata <- list()

for (tissue in tissues){
    data[[tissue]] <- readRDS(paste0("tissues/", tissue, "/methylation_data.rds"))
    metadata[[tissue]] <- readRDS(paste0("tissues/", tissue, "/methylation_metadata_subset.rds"))
}


#Process the data
# Keep only subject-ID part
data <- lapply(tissues, function(tissue){
  tissue_data <- data[[tissue]]
  colnames(tissue_data) <- gsub("-[A-Za-z0-9]{4,5}-SM-.*", "\\1", colnames(tissue_data))
  
  tissue_data
})

names(data) <- tissues

# Filter for samples with age and smoking data ----
data <- lapply(tissues, function(tissue) data[[tissue]] %>% 
                 select("probe", all_of(metadata[[tissue]] %>% pull("SUBJID")))
               )

names(data) <- tissues

# Keep only probes in the 21K (we have some missing values, which we will set to NA)
data21k <- lapply(tissues, function(tissue) data[[tissue]] %>% 
                    filter(probe %in% probeAnnotation21kdatMethUsed$Name)
)

names(data21k) <- tissues

# Due to missing genes from the 21K probe array, we need to input them with NA values.
# Imputation is based on average values across samples
data21kImputed <- lapply(tissues, function(tissue){
  tissue_data <- data21k[[tissue]] 
  
  match1 <- match(probeAnnotation21kdatMethUsed$Name , tissue_data$probe)
  missing_probes <- probeAnnotation21kdatMethUsed$Name[is.na(match1)]
  missing_probes_data <- matrix(data = NA, nrow = length(missing_probes), ncol = ncol(tissue_data))
  missing_probes_data[,1] <- missing_probes
  colnames(missing_probes_data) <-  colnames(tissue_data)
  
  missing_probes_data <- as.data.frame(missing_probes_data)
  
  tissue_data <- rbind(tissue_data, missing_probes_data)
  tissue_data <- tissue_data %>% 
    as_tibble() %>%
    mutate_at(2:ncol(tissue_data), as.numeric)
})

names(data21kImputed) <- tissues

## Run pipeline ----

### Step 1: Impute missing data ----
data21kImputed <- lapply(tissues, function(tissue){
  tissue_data <- data21kImputed[[tissue]]
  
  datMethUsed <- t(tissue_data[,-1])
  colnames(datMethUsed)=as.character(tissue_data$probe)
  
  
  #Note There another faster imputation method in the Horvath source code
  dimnames1=dimnames(datMethUsed)
  datMethUsed= data.frame(t(impute.knn(t(datMethUsed))$data))
  dimnames(datMethUsed)=dimnames1
  
  datMethUsed
})

names(data21kImputed) <- tissues


### Step2: Normalize data ----
data21kImputed <- lapply(tissues, function(tissue) BMIQcalibration(
    datM = data21kImputed[[tissue]],
    goldstandard.beta = probeAnnotation21kdatMethUsed$goldstandard2,
    plots = FALSE
  )
)

names(data21kImputed) <- tissues

### Step 3: Predict age with Horvath clock ----
selectCpGsClock <- is.element(dimnames(data21kImputed$Lung)[[2]], as.character(datClock$CpGmarker[-1]))

predicted_age_horvath <- lapply(tissues, function(tissue){
  datMethClock0 <- data.frame(data21kImputed[[tissue]][,selectCpGsClock])
  datMethClock <-  data.frame(datMethClock0[as.character(datClock$CpGmarker[-1])])
  predictedAge <- as.numeric(anti.trafo(datClock$CoefficientTraining[1]+as.matrix(datMethClock)%*% as.numeric(datClock$CoefficientTraining[-1])))
  
  data_to_return <- data.frame("SUBJECT" =row.names(datMethClock0), 
                               "meth_age" = predictedAge)
  
  data_to_return
})


names(predicted_age_horvath) <- tissues

### Compare Smokers vs never smokers----
# Are Smokers predicted with an higher error compared with never-smokers? 

# Are Smokers with a higher error compared with never-smokers, if we adjust for 
# Other factors? 
age_error <- data.frame("tissue" = c(), "error" = c(), "smoking_status" = c())
residuals_regression <- data.frame("tissue" = c(), "age" = c(), "residual" = c(), "smoking_status" = c())

variables <- c("Ancestry", "SEX", "AGE", "BMI", "TRISCHD", "DTHHRDY",
                          "Smoking")

horvath_coef_table <- c()

for (tissue in tissues){
  ages_pred <- predicted_age_horvath[[tissue]]
  ages_real <- metadata[[tissue]]
  
  ages <- merge(ages_pred, ages_real, by.x = "SUBJECT" , by.y = "SUBJID")
  
  ages$error <- ages$meth_age - ages$AGE
  
  # Divide into smokers and never smokers
  error_smokers <- ages %>% filter(Smoking == 0) %>% pull(error)
  error_ns <- ages %>% filter(Smoking == 2) %>% pull(error)
  
  #Add to the error dataframe 
  age_error <- rbind(age_error, data.frame("tissue" = tissue, "error" = error_smokers, "smoking_status" = "Smoker"))
  age_error <- rbind(age_error, data.frame("tissue" = tissue, "error" = error_ns,"smoking_status" = "Never Smoker"))
  

  # Correct for other factors such as BMI, sex and others 
  ages.s.ns <- ages %>% 
    filter(Smoking != 1) %>% 
    mutate_at(.vars = c("Smoking", "DTHHRDY", "SEX", "Ancestry"), as.factor)

  #if (tissue %in% sex_tissues){
  #  variables_2 <- variables[variables != "SEX"]
  #  formula <- as.formula(paste("error ~", paste(variables_2, collapse = " + ")))
  #  variables3 <- variables_2[variables_2 != "Smoking"] 
  #  formula3 <- as.formula(paste("error ~", paste(variables3, collapse = " + ")))
    
  #}else{
   formula <- as.formula(paste("error ~", paste(variables, collapse = " + ")))
   variables3 <- variables[variables != "Smoking"] 
   formula3 <- as.formula(paste("error ~", paste(variables3, collapse = " + ")))
    
  #}
  regression_horvath <- lm(formula, data = ages.s.ns)
  
  # Save residuals from the regression to plot
  regression_horvath_2 <- lm(formula3, data = ages.s.ns)
  residuals <- resid(regression_horvath_2)
  residuals_regression <- rbind(residuals_regression, 
                                data.frame("tissue" = tissue, 
                                           "age" = ages.s.ns$AGE,
                                           "error" = residuals,
                                           "smoking_status" = ages.s.ns$Smoking)
                                )

  # Save the p.value of the coefficients in a table
  horvath_coef_table <- rbind(horvath_coef_table, 
                             c(tissue, 
                                summary(regression_horvath)$coefficients["Smoking2","Pr(>|t|)"])
  )
}


colnames(horvath_coef_table) <- c("tissue", "p.value")
horvath_coef_table <- horvath_coef_table %>%
  as.data.frame() %>% 
  mutate(p.value = as.numeric(p.value)) %>% 
  mutate(FDR = p.adjust(p.value, method = "fdr"))


data_to_export <- list()
data_to_export$Lung_residuals_regression <- residuals_regression %>% filter(tissue == "Lung")

saveRDS(to_plot_f, "../figures/data/figure_clock_hovath.rds")

# Vadim clocks ----
## Load molecular clocks coefficients for Vadim clocks----
trafo= function(x,adult.age=20) { x=(x+1)/(1+adult.age); y=ifelse(x<=1, log( x),x-1);y }
anti.trafo= function(x,adult.age=20) { ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }
probeAnnotation21kdatMethUsed=read.csv("data/public/clockData/probeAnnotation21kdatMethUsed.csv")
probeAnnotation27k=read.csv("data/public/clockData/datMiniAnnotation27k.csv")
datClock=read.csv("data/public/clockData/AdditionalFile3.csv")


clocks <- list()
clocks$CauAge <- read.csv("data/public/clockData/YingCausAge.csv")
clocks$DamAge <- read.csv("data/public/clockData/YingDamAge.csv")
clocks$AdaptAge <- read.csv("data/public/clockData/YingAdaptAge.csv")

# Check for missing data across molecular clocks ----
for (clock in names(clocks)){
  cat("Clock ", clock, "has ", sum(!clocks[[clock]]$term[-1] %in% data[["Lung"]]$probe), 
      "missing probes from ", nrow(clocks[[clock]]), "\n")
}

# CauAge has too many missing samples

clocks_to_run <- clocks
clocks_to_run$CauAge <- NULL

# Keep data only from the 450K array (same array used to train the Vadim Clocks)
array450k <- read.csv("data/public/clockData/humanmethylation450_15017482_v1-2.csv")

data450k <- lapply(tissues, function(tissue) data[[tissue]] %>% 
                    filter(probe %in% array450k$Name)
)

names(data450k) <- tissues

# Run pipeline ----
# Due to missing genes from the 450K probe array, we need to input them.
# Imputation is based on average methylation values across the samples 
# This is based on Horvath methodology to deal with missing data in the online cal
# Due to computational limits, I only used the the probes in the model as background 
# Using the 450 could potencial improve something, but I not sure the best way to 
# do it is. Horvath uses all the probes in the array (21k)


## Adapt Age----

# Add missing probes with value "NA"
example_data <- data450k$Lung
missing_probes <- clocks_to_run$AdaptAge$term[!clocks_to_run$AdaptAge$term %in% data450k$Lung$probe]
missing_probes <- missing_probes[-1] # removing Intercept


predictedAdaptAge <- lapply(tissues, function(tissue){
  tissue_data <- data450k[[tissue]] 

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
  
  data_to_return <- data.frame("SUBJECT" =row.names(datMethClock), 
                               "Adapt_age" = predictedAge)
  
  data_to_return
})


names(predictedAdaptAge) <- tissues

## DamAge Age----

# Add missing probes with value "NA"
example_data <- data450k$Lung
missing_probes <- clocks_to_run$DamAge$term[!clocks_to_run$DamAge$term %in% data450k$Lung$probe]
missing_probes <- missing_probes[-1] # removing Intercept

predictedDamAge <- lapply(tissues, function(tissue){
  tissue_data <- data450k[[tissue]] 
  
  
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
  
  data_to_return <- data.frame("SUBJECT" =row.names(datMethClock), 
                               "Dam_age" = predictedAge)
  
  data_to_return
})

names(predictedDamAge) <- tissues

### Comparison between non-smoker and smokers ----
#Are Smokers predicted with an higher error compared with never-smokers? 
#Do smokers have a higher error compared with never-smokers, if we adjust for 
# Other factors? 

age_error_adaptAge <- data.frame("tissue" = c(), "error" = c(), "smoking_status" = c())
age_error_DamAge <- data.frame("tissue" = c(), "error" = c(), "smoking_status" = c())

variables <- c("Ancestry", "SEX", "AGE", "BMI", "TRISCHD", "DTHHRDY",
               "Smoking")

AdaptAge_coef_table <- c()
DamAge_coef_table <- c()


for (tissue in tissues){
  ages_predAdapt <- predictedAdaptAge[[tissue]]
  ages_predDam <- predictedDamAge[[tissue]]
  ages_pred <- merge(ages_predAdapt, ages_predDam, by = "SUBJECT")
  
  ages_real <- metadata[[tissue]]
  
  ages <- merge(ages_pred, ages_real, by.x = "SUBJECT" , by.y = "SUBJID")
  
  ages$errorAdaptAge <- ages$Adapt_age - ages$AGE
  ages$errorDamAge <- ages$Dam_age - ages$AGE
  
  ## For AdaptAge
  # Divide into smokers and never smokers
  error_smokers <- ages %>% filter(Smoking == 0) %>% pull(errorAdaptAge)
  error_ns <- ages %>% filter(Smoking == 2) %>% pull(errorAdaptAge)
  
  #Add to the error dataframe 
  age_error_adaptAge <- rbind(age_error_adaptAge, data.frame("tissue" = tissue, "error" = error_smokers, "smoking_status" = "Smoker"))
  age_error_adaptAge <- rbind(age_error_adaptAge, data.frame("tissue" = tissue, "error" = error_ns,"smoking_status" = "Never Smoker"))
  
  # Correct for other factors such as BMI, sex and others 
  ages.s.ns <- ages %>% 
    filter(Smoking != 1) %>% 
    mutate_at(.vars = c("Smoking", "DTHHRDY", "SEX"), as.factor)
  
  #if (tissue %in% sex_tissues){
  #  variables_2 <- variables[variables != "SEX"]
  #  formula <- as.formula(paste("errorAdaptAge  ~", paste(variables_2, collapse = " + ")))
    
  #}else{
  formula <- as.formula(paste("errorAdaptAge  ~", paste(variables, collapse = " + "))) 
  #}
  regression <- lm(formula, data = ages.s.ns)
  
  # Save the p.value of the coefficients in a table
  AdaptAge_coef_table <- rbind(AdaptAge_coef_table, c(tissue, 
                                                    summary(regression)$coefficients["Smoking2","Pr(>|t|)"],
                                                    summary(regression)$coefficients["Smoking2","Estimate"])
  )
  
  ## For DamAge 
  # Divide into smokers and never smokers
  error_smokers <- ages %>% filter(Smoking == 0) %>% pull(errorDamAge)
  error_ns <- ages %>% filter(Smoking == 2) %>% pull(errorDamAge)
  
  #Add to the error dataframe 
  age_error_DamAge <- rbind(age_error_DamAge, data.frame("tissue" = tissue, "error" = error_smokers, "smoking_status" = "Smoker"))
  age_error_DamAge <- rbind(age_error_DamAge, data.frame("tissue" = tissue, "error" = error_ns,"smoking_status" = "Never Smoker"))
  

  # Correct for other factors such as BMI, sex and others 
  ages.s.ns <- ages %>% 
    filter(Smoking != 1) %>% 
    mutate_at(.vars = c("Smoking", "DTHHRDY", "SEX"), as.factor)
  
  #if (tissue %in% sex_tissues){
  #  variables_2 <- variables[variables != "SEX"]
  #  formula <- as.formula(paste("errorDamAge  ~", paste(variables_2, collapse = " + ")))
    
  #}else{
  formula <- as.formula(paste("errorDamAge ~", paste(variables, collapse = " + ")))
  #}
  
  regression <- lm(formula, data = ages.s.ns)
  
  # Save the p.value of the coefficients in a table
  DamAge_coef_table <- rbind(DamAge_coef_table, c(tissue, 
                                                  summary(regression)$coefficients["Smoking2","Pr(>|t|)"],
                                                  summary(regression)$coefficients["Smoking2","Estimate"])
  )
}

colnames(DamAge_coef_table) <- c("tissue", "p.value", "estimate")
colnames(AdaptAge_coef_table) <- c("tissue", "p.value", "estimate")

DamAge_coef_table <- DamAge_coef_table %>%
  as.data.frame() %>% 
  mutate(p.value = as.numeric(p.value)) %>% 
  mutate(FDR = p.adjust(p.value, method = "fdr"))

AdaptAge_coef_table <- AdaptAge_coef_table %>%
  as.data.frame() %>% 
  mutate(p.value = as.numeric(p.value)) %>% 
  mutate(FDR = p.adjust(p.value, method = "fdr"))