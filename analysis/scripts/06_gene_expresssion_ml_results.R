suppressPackageStartupMessages(library(tidyverse))
library(lightgbm)

biotype <- read.delim("data/public/gencode.v26.GRCh38.genes.bed")


## --------------------------------------------------------------------------------------------------------------------------------------------------------
get_metric_tables <- function(metrics){
  metric_table <- data.frame()
  for (metric in names(metrics)){
    mean_metric <- round(mean(metrics[[metric]]), 3)
    sd_metric <- round(sd(metrics[[metric]]), 3)
    metric_table <- rbind(metric_table, c(metric, mean_metric, sd_metric))
    print(paste0("Metric ", metric, " had an average of ", mean_metric, " +/- ", sd_metric))
  }  
  colnames(metric_table) <- c("metrics","avg", "sd")
  metric_table
}


#Average metrics (accuracy, F1 score, etc)
lung_metrics <- readRDS(file = "output/ML_expression/Lung/5cv_metrics.rds")
lung_metrics_avg <- get_metric_tables(lung_metrics) 
lung_metrics_avg$tissue <- "Lung"

#Covariate correlation in NS and S 
covar.cor.lung <- readRDS(file = "output/ML_expression/Lung/5cv_cor_covars.rds")

#Feature importance 
fetureImportance_lung <- readRDS("output/ML_expression/Lung/feature_importance_complete_model") %>% 
  merge(biotype, by.x = "Feature", by.y = "ensembl") %>% 
  arrange(-Gain)

#Ex-smokers classification
ex_smokers_classification.lung <- readRDS(file = "output/ML_expression/Lung/exSmoker_lung_prediction.rds")
ex_smokers_classification.lung$predRound <- round(ex_smokers_classification.lung$pred)
ex_smokers_classification.lung$tissue <- "Lung"


## --------------------------------------------------------------------------------------------------------------------------------------------------------
#Average metrics (accuracy, F1 score, etc)
thyroid_metrics <- readRDS(file = "output/ML_expression/Thyroid/5cv_metrics.rds")
thyroid_metrics_avg <- get_metric_tables(thyroid_metrics)
thyroid_metrics_avg$tissue <- "Thyroid"

#Covariate correlation in NS and S 
covar.cor.thyroid <- readRDS(file = "output/ML_expression/Thyroid/5cv_cor_covars.rds")

#Feature importance 
fetureImportance_thyroid <- readRDS("output/ML_expression/Thyroid/feature_importance_complete_model") %>% 
  merge(biotype, by.x = "Feature", by.y = "ensembl") %>% 
  arrange(-Gain)

#Ex-smokers classification
ex_smokers_classification.thyroid <- readRDS(file = "output/ML_expression/thyroid/exSmoker_lung_prediction.rds")
ex_smokers_classification.thyroid$predRound <- round(ex_smokers_classification.thyroid$pred)
ex_smokers_classification.thyroid$tissue <- "Thyroid"


## --------------------------------------------------------------------------------------------------------------------------------------------------------
#Average metrics (accuracy, F1 score, etc)
esophagus_muc_metrics <- readRDS(file = "output/ML_expression/EsophagusMucosa/5cv_metrics.rds")
esophagus_muc_metrics_avg <- get_metric_tables(esophagus_muc_metrics)
esophagus_muc_metrics_avg$tissue <- "Esophagus Mucosa"

#Covariate correlation in NS and S 
covar.cor.esophagus_muc <- readRDS(file = "output/ML_expression/EsophagusMucosa/5cv_cor_covars.rds")

#Feature importance 
fetureImportance_esophagus_muc <- readRDS("output/ML_expression/EsophagusMucosa/feature_importance_complete_model") %>% 
  merge(biotype, by.x = "Feature", by.y = "ensembl") %>% 
  arrange(-Gain)

#Ex-smokers classification
ex_smokers_classification.esophagus_muc <- readRDS(file = "output/ML_expression/EsophagusMucosa/exSmoker_lung_prediction.rds")
ex_smokers_classification.esophagus_muc$predRound <- round(ex_smokers_classification.esophagus_muc$pred)
ex_smokers_classification.esophagus_muc$tissue <- "Esophagus Mucosa"


## --------------------------------------------------------------------------------------------------------------------------------------------------------
metrics <- rbind(lung_metrics_avg, rbind(thyroid_metrics_avg, esophagus_muc_metrics_avg))
ex_smoker_classification <- rbind(ex_smokers_classification.lung, rbind(ex_smokers_classification.thyroid, ex_smokers_classification.esophagus_muc))


## --------------------------------------------------------------------------------------------------------------------------------------------------------
ml_results_summary <- list()
ml_results_summary$metrics <- metrics
ml_results_summary$lung_feature_importance <- fetureImportance_lung
ml_results_summary$thyroid_feature_importance <- fetureImportance_thyroid
ml_results_summary$esophagus_mucosa_feature_importance <- fetureImportance_esophagus_muc
ml_results_summary$ex_smokers_classification <- ex_smoker_classification

saveRDS(ml_results_summary, file = "../figures/data/ML_expression.rds")

