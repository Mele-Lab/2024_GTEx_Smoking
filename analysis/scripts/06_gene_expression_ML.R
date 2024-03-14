#!/usr/bin/env Rscript
# @Author: Rogério Eduardo Ramos Ribeiro
# @E-mail: zroger499@gmail.com
# @Description: Code to run ML in gene expression (Lung + Thyroid + Esophagus mucosa)
# @software version: R=4.2.2

####The strategy is to build models and validate using 5CV on smokers and never smokers and then using them on ex-smokers###

Sys.time()

# Load libraries

suppressPackageStartupMessages(library(tidyverse))
library(caret)
library(lightgbm)
library(splitTools)

#Make working dir 
setwd(system("pwd", intern = T)) #If in linux

#Set SEED 
set.seed("42")

Sys.time()

# Load data 
expression_data <- read.delim("data/expression/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz", skip = 2)
metadata <- readRDS(file = "data/metadata/metadata.rds")
biotype <- readRDS(file = "data/metadata/biotype.rds")

#Subset expression dataset for pcoding and lincRNA

selected.genes <- biotype %>% 
  filter(biotype == "protein_coding" | biotype == "lincRNA") %>% 
  pull(ensembl)

expression_data <- expression_data %>% 
  filter(Name %in% selected.genes)


cat("The expression dataset has ", nrow(expression_data)," genes\n")

##Analysis function (saves everything in intermediary files)
out.folder <- "output/ML_expression/"
if (!dir.exists(out.folder)){
  dir.create(out.folder, recursive = T)
}

tissues <- c("Lung", "Thyroid", "EsophagusMucosa")


# 5CV cross validation
perform_cross_validation <- function(folds, tissue_expression, labels, covars.n.ns){
  
  metrics <- c() 
  metrics$Accuracy <- c()
  metrics$Specificity <- c()
  metrics$Sensitivity <- c()
  metrics$Precision <- c()
  metrics$F1 <- c()
  
  
  cor_covars <- c()
  for (j in 1:length(folds)){
    fold <- folds[[j]]
    fold_name <- paste0("fold_", j)
    
    dtrain.fold <- lgb.Dataset(data = tissue_expression[fold, ], label = labels[fold], params = list(feature_pre_filter=FALSE))
    dtest.fold <- lgb.Dataset.create.valid(dtrain.fold, data = tissue_expression[-fold, ], label = labels[-fold], params = list(feature_pre_filter=FALSE))
    
    model <- lgb.train(
      list(
        objective  = "binary",
        metric     = "binary_error",
        nthread = 8,
        num_class = 1L),
      data = dtrain.fold,
      verbose= -1)
    
    
    my_preds <- predict(model, tissue_expression[-fold, ])
    label.test  <- labels[-fold]
    my_preds_rnd = round(my_preds)
    confusion_matrix <- confusionMatrix(data = as.factor(my_preds_rnd), reference = as.factor(label.test))
    
    metrics$Accuracy[j]<- confusion_matrix$overall[1]
    metrics$Specificity[j] <- confusion_matrix$byClass[2]
    metrics$Sensitivity[j] <-  confusion_matrix$byClass[1]
    metrics$Precision[j] <-confusion_matrix$byClass[5]
    metrics$F1[j] <- confusion_matrix$byClass[7]
    
    
    #Correlation with covariates in the test set
    for (i in 1:ncol(covars.n.ns)){
      data <- covars.n.ns[-fold, i]
      var <- colnames(data)[1]
      data = unlist(unname(data))
      test <-  cor.test(as.numeric(data), my_preds, method = "spearman")
      correlation <-test$estimate
      corraltion_pvalue <-test$p.value
      cor_covars <- rbind(cor_covars, c(correlation, corraltion_pvalue, var, j))
    }
    
  }
  
  colnames(cor_covars) <- c("rho", "pvalue", "var", "interation")
  cor_covars <- cor_covars %>% 
    as_tibble() %>% 
    mutate_at(.vars = c(1,2,4), as.numeric) %>% 
    mutate(padjust = p.adjust(pvalue))
  
  res <- c()
  res$metrics <- metrics
  res$cor_covars <- cor_covars
  res
}


# Main analysis function
run_ML <- function(exp, metadata, tissues, out.folder, cv = 5){
  #Run Machine learning analysis
  
  #Parse colnames
  colnames(exp) <- gsub("\\.SM.*", "", colnames(exp))
  colnames(exp) <- gsub("\\.", "-", colnames(exp))
  row.names(exp) <- exp$Name
  
  for (tissue in tissues){
    cat("Starting analysis on ", tissue, "\n")
    #create out folder 
    tissue_folder <- paste0(out.folder, "/", tissue, "/")
    if (!dir.exists(tissue_folder)){
      dir.create(tissue_folder)
    }
    
    tissue_samples <- metadata[[tissue]]
    
    #Subset smoker samples
    smk_never_smoker_samples <- tissue_samples %>% 
      filter(Smoking == 0 | Smoking == 2)
    
    #Subset expression matrix in smoker_samples
    exp_tissue <- exp %>% 
      select(all_of(smk_never_smoker_samples$Sample))
    
    cat("Starting ML analysis on ", tissue, " with ", ncol(exp_tissue), " samples (NS + S) and ", nrow(exp_tissue), " genes\n")
    
    # Log input data
    exp_tissue.log <- log2(exp_tissue + 1)
    
    # filter genes by median expression
    gene.median.tpm = apply(exp_tissue.log,1, median)
    exp_tissue.log.sel = exp_tissue.log[names(which(gene.median.tpm > 1)),] #Filter for median > 1 across samples
    
    # Transpose the matrix
    exp_tissue.log.sel = t(exp_tissue.log.sel)
    
    # Extract labels
    labels <- as.character(smk_never_smoker_samples$Smoking)
    labels <- as.numeric(as.factor(droplevels.factor(labels))) - 1 #1 are smokers, 0 are never smokers
    
    print(paste0("Input dataset variables has dimensions ", dim(exp_tissue.log.sel)))
    print("Running sample check ", all(smk_never_smoker_samples$Sample == row.names(exp_tissue.log.sel)))
    
    folds <- create_folds(labels, k = cv)
    
    #Create a covariate table (for correlations with prediction scores)
    covars.n.ns <-  smk_never_smoker_samples %>% 
      select(-c(Donor, Sample)) %>% 
      separate_rows(Ancestry, sep = ", ") %>% 
      mutate(dummy = 1) %>% 
      pivot_wider(names_from = Ancestry, values_from = dummy, values_fill = 0) %>% 
      mutate(HardyScale = paste0("Hardy", HardyScale)) %>%
      separate_rows(HardyScale, sep = ", ") %>% 
      mutate(dummy = 1) %>% 
      pivot_wider(names_from = HardyScale, values_from = dummy, values_fill = 0) %>% 
      mutate(Smoking = as.numeric(as.factor(Smoking)) -1)
    
    
    cat("Start running cross validation on Smokers and Never smokers")
    cv_result <- perform_cross_validation(folds, exp_tissue.log.sel,labels, covars.n.ns)
    
    #Save metrics and correlation with covars
    saveRDS(cv_result$metrics,paste0(tissue_folder, "5cv_metrics.rds"))
    saveRDS(cv_result$cor_covars,paste0(tissue_folder, "5cv_cor_covars.rds"))
    
    
    #Train the model on the complete dataset
    dtrain.complete <- lgb.Dataset(data = exp_tissue.log.sel, label = labels, params = list(feature_pre_filter=FALSE))
    model_smokers_neversmokers.complete <- lgb.train(
      list(
        objective  = "binary",
        metric     = "binary_error",
        nthread = 8,
        num_class = 1L),
      data = dtrain.complete)
    
    lgb.save(model_smokers_neversmokers.complete, paste0(tissue_folder, "model.smokers_neversmokers.complete"))
    
    #Get feature importance
    feature_importance <- lgb.importance(model_smokers_neversmokers.complete, percentage = TRUE)
    saveRDS(feature_importance,paste0(tissue_folder, "feature_importance_complete_model"))
    
    # Predict ex-smokers using the complete model
    # select only ex-smokers for classification
    ex_samples <- tissue_samples %>% 
      filter(Smoking == 1)
      
    exp_tissue_ex <- exp %>% 
      select(all_of(ex_samples$Sample))
  
    # create expression table 
    selected.genes =  colnames(exp_tissue.log.sel)
    exp_tissue.log.sel.ex = log2(1 + exp_tissue_ex[as.character(selected.genes),])
    exp_tissue.log.sel.ex = t(exp_tissue.log.sel.ex)
    
    cat("Predicting Smk status of ", nrow(exp_tissue.log.sel.ex), " Ex smokers\n")
    
    #Predict
    my_preds.exSmokers <- predict(model_smokers_neversmokers.complete, exp_tissue.log.sel.ex)
    
    #Save smoking predictions
    ex_smoker_percent <- data.frame("pred" = my_preds.exSmokers, "samples" = row.names( exp_tissue.log.sel.ex))
    saveRDS(ex_smoker_percent, paste0(tissue_folder, "exSmoker_lung_prediction.rds"))
  }
}


run_ML(expression_data, metadata, tissues, out.folder, cv = 5)
Sys.time()
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
