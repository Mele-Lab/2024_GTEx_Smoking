####################################################################
# Libraries
#####################################################################
library(lightgbm)
library(tidyverse)
library(caret)
library(matrixStats)
library(splitTools)

setwd(system("pwd", intern = T)) #If in linux
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio

# for lightgbm do not crash
Sys.setenv(KMP_DUPLICATE_LIB_OK = TRUE)

data_to_export <- list()

#Load data
smoking_annot <- readRDS("output/metadata_smoking.rds")
smoking_annot <- smoking_annot %>% mutate(Smoking = ifelse(Smoking == 0, "NeverSmoker", ifelse(Smoking == 1, "ExSmoker", "Smoker")))
row.names(smoking_annot) <- smoking_annot$Sample

# Load array annotation 
array_annotation <- read.csv(file = "public/methylation_epic_v1.0b5.csv")

#Seed to use
SEED = 3456

# Load the methylation EPIC array 
## Methylation table avaiable at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213478
epic_array_lung <- data.table::fread("public/GSE213478_methylation_DNAm_noob_final_BMIQ_all_tissues_987.txt.gz", sep=",")
epic_array <- epic_array %>% column_to_rownames("V1")


# Load metadata and subset lung samples
metadata <- read.csv("public/eGTExDNA_Pierce_Jan18.09-11-2021.tsv", sep = "\t")
metadata_lung_samples <- metadata %>% filter(Tissue.Site.Detail == "Lung") %>% pull( Sample.ID.for.data.sharing.and.public.release)

# Subset lung samples 
epic_array_lung <- epic_array %>% select(all_of(metadata_lung_samples))


# Tranform lung sample names 
colnames(epic_array_lung) <- str_extract(colnames(epic_array_lung), "GTEX-.{4,5}")
colnames(epic_array_lung) <- gsub("-$", "", colnames(epic_array_lung))


# Transpose the matrix 
epic_array_lung <- t(epic_array_lung)



#Lung ids and smoking annotation 
smoking_labels_smokers_ns <- smoking_annot %>% 
  filter(Smoking == "Smoker" | Smoking == "NeverSmoker") %>% 
  pull("SUBJID")

smoking_labels_smokers_exs <- smoking_annot %>% 
  filter(Smoking == "ExSmoker") %>% 
  pull("SUBJID")


#Filter for only smoking and never smoking samples
epic_array_lung_ns_s <- epic_array_lung[as.character(smoking_labels_smokers_ns), ]
epic_array_lung_exs <- epic_array_lung[as.character(smoking_labels_smokers_exs), ]


#Reoder the labels
cat("All the labels in the correct order:",  all(smoking_labels_smokers_ns == row.names(epic_array_lung_ns_s)))

# make the labels numeric
lung_smoking_labels <- as.numeric(as.factor(smoking_annot %>% filter(Smoking == "Smoker" | Smoking == "NeverSmoker") %>% pull("Smoking"))) - 1


#####################################
##Perform 5-fold cross validation ###
####################################

perform_cross_validation <- function(folds){
  metrics <- c() 
  metrics$Accuracy <- c()
  metrics$Specificity <- c()
  metrics$Sensitivity <- c()
  metrics$Precision <- c()
  metrics$F1 <- c()
  
  cat("Running folds .. ")
  for (j in 1:length(folds)){
    cat(j, ".. ")
    fold <- folds[[j]]
    fold_name <- paste0("fold_", j)
    dtrain.fold <- lgb.Dataset(data = epic_array_lung_ns_s[fold, ], label = lung_smoking_labels[fold], params = list(feature_pre_filter=FALSE))
    dtest.fold <- lgb.Dataset.create.valid(dtrain.fold, data = epic_array_lung_ns_s[-fold, ], label = lung_smoking_labels[-fold], params = list(feature_pre_filter=FALSE))
    
    model <- lgb.train(
      list(
        objective  = "binary",
        metric     = "binary_error",
        nthread = 2,
        num_class = 1L),
      data = dtrain.fold,
      verbose= -1)
    
    my_preds <- predict(model, epic_array_lung_ns_s[-fold, ])
    label.test  <- lung_smoking_labels[-fold]
    my_preds_rnd = round(my_preds)
    
    
    confusion_matrix <- confusionMatrix(data = as.factor(my_preds_rnd), reference = as.factor(label.test))
    
    metrics$Accuracy[j]<- confusion_matrix$overall[1]
    metrics$Specificity[j] <- confusion_matrix$byClass[2]
    metrics$Sensitivity[j] <-  confusion_matrix$byClass[1]
    metrics$Precision[j] <-confusion_matrix$byClass[5]
    metrics$F1[j] <- confusion_matrix$byClass[7]
    
  }
  
  return (metrics)
}


epic_array_lung_ns_s <- as.matrix(epic_array_lung_ns_s)

set.seed(SEED)
n_fold = 5
folds <- create_folds(lung_smoking_labels, k = n_fold)
metrics <- perform_cross_validation(folds)

### Average across validation scores 

metric_table <- data.frame()
for (metric in names(metrics)){
  mean_metric <- round(mean(metrics[[metric]]), 3)
  sd_metric <- round(sd(metrics[[metric]]), 3)
  metric_table <- rbind(metric_table, c(metric, mean_metric, sd_metric))
  print(paste0("Metric ", metric, " had an average of ", mean_metric, " +/- ", sd_metric))
}

colnames(metric_table) <- c("metric", "mean", "sd")
metric_table[1,1] <- "Accuracy"
metric_table <- metric_table %>% 
  mutate(mean = as.numeric(mean)) %>% 
  mutate(sd = as.numeric(sd)) %>% 
  filter(metric != "Recall") %>%
  mutate(metric = factor(metric, levels = metric[length(metric):1]))


safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")


# To predict ex-smoker retrain the whole model with the complete dataframe (smokers + ex-smokers)
dtrain.complete <- lgb.Dataset(data = epic_array_lung_ns_s, label = lung_smoking_labels, params = list(feature_pre_filter=FALSE))
model_smokers_neversmokers.complete <- lgb.train(
  list(
    objective  = "binary",
    metric     = "binary_error",
    nthread = 8,
    num_class = 1L),
  data = dtrain.complete)


## Check feature importance
###Feature importance
tree_importance.complete <- lgb.importance(model_smokers_neversmokers.complete, percentage = TRUE)
tree_importance.complete.genes <- merge(tree_importance.complete, array_annotation, by.x = "Feature", by.y = "Name")

tree_importance.complete.genes.top15 <- tree_importance.complete.genes %>% 
  arrange(Gain) %>% 
  slice_tail(n = 15) %>% 
  mutate(Gene_name = ifelse(GencodeBasicV12_NAME == "", Feature, GencodeBasicV12_NAME)) %>%
  mutate(Feature = factor(Feature, levels = Feature[1:length(Feature)]))

tree_importance.complete.genes.top15$Gene_name <- gsub(";.*", "", tree_importance.complete.genes.top15$Gene_name)
tree_importance.complete.genes.top15$Gene_name <- ifelse(tree_importance.complete.genes.top15$Gene_name == tree_importance.complete.genes.top15$Feature,
                                                         tree_importance.complete.genes.top15$Gene_name, 
                                                         paste0(tree_importance.complete.genes.top15$Feature, " - ", tree_importance.complete.genes.top15$Gene_name))



data_to_export[["feature_importance"]] <- tree_importance.complete.genes.top15


# Classify ex-smokers
## predictions for Ex Smokers using the smokers_neversmokers model
my_preds.exSmokers <- predict(model_smokers_neversmokers.complete, as.matrix(epic_array_lung_exs))
my_preds_rnd.exSmokers = round(my_preds.exSmokers)

#Make a dataframe with the prediction 
ex_smoker_pred_prediction <- data.frame(data = row.names(epic_array_lung_exs), pred = my_preds.exSmokers)

# add the classification labels to sample ids and create data frame
exsmokers_preds = as.data.frame(cbind(rownames(epic_array_lung_exs), my_preds_rnd.exSmokers))
names(exsmokers_preds) = c("sample","SmokePrediction")
exsmokers_preds$ex_smoker_pred <- ifelse(exsmokers_preds$SmokePrediction == 1, "Smoker", "Never Smoker")


ex_smoker_pred_count <- as.data.frame(table(exsmokers_preds$SmokePrediction))
ex_smoker_pred_count$Var1 <- c("Never Smoker", "Smoker")
ex_smoker_pred_count$tissue <- "Lung"


data_to_export$ex_smoker_predictions <- exsmokers_preds

saveRDS(data_to_export, file = "../figures/data/methylation_ml.rds")