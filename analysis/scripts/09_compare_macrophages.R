library(BayesPrism)
library(tidyverse)
library(betareg)
library(multcomp)

#Load color pallete to use
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")


ted <- readRDS(file = "output/scRNA_seq/bayesPrism.lung.signature_genes.rds")
metadata_lung <- readRDS(file = "output/metadata/metadata.rds")$Lung

theta <- get.fraction (bp=ted,
            which.theta="final",
            state.or.type="type")



# Main computing function

#' @param data Data to run the PMI. Includes dependent and independent  variables
#' @param celltype Celltype to test
#' @param formula Formula to give to the independent Probablistic model
da_beta_regression <- function(data, celltype, formula){
    res <- c()
    subsetData <- data %>% filter(Celltype == celltype)
    median_all <- median(subsetData$value)
   
    #Summarise the median differences
    medians <- subsetData %>% 
        group_by(Smoking) %>% 
        summarise(median = median(value))

    #Extract each group median
    median_smoker <- medians %>% filter(Smoking == 2) %>% pull(median)
    median_non_smoker <- medians %>% filter(Smoking == 0) %>% pull(median)
    median_ex_smoker <- medians %>% filter(Smoking == 1) %>% pull(median)

    #compare medians
    median_dif_smoker_non_smoker <- median_smoker - median_non_smoker
    median_fc_smoker_non_smoker <- median_smoker/median_non_smoker

    median_dif_ex_smokers_non_smoker <- median_ex_smoker - median_non_smoker
    median_fc_ex_smokers_non_smoker <- median_ex_smoker/median_non_smoker

    median_dif_smoker_ex_smoker <- median_smoker - median_ex_smoker
    median_fc_smoker_ex_smoker <- median_smoker/median_ex_smoker

    #Run the beta regression model 
    model <- betareg(design, data=subsetData, link = "loglog") 

    #Hypothesis testing with glht
    mult <- summary(glht(model, linfct = c("Smoking2 = 0", "Smoking1 = 0", "Smoking2-Smoking1 = 0")))
    mult$test$pvalues

    res <- rbind(res, c(name,  
                        median_all, median_non_smoker, median_ex_smoker, median_smoker,
                        median_dif_smoker_non_smoker, median_fc_smoker_non_smoker, mult$test$pvalues[1],
                        median_dif_ex_smokers_non_smoker, median_fc_ex_smokers_non_smoker, mult$test$pvalues[2],
                        median_dif_smoker_ex_smoker, median_fc_smoker_ex_smoker, mult$test$pvalues[3]))


  colnames(res) <- c("celltype", 
                     "median","median_never_smoker", "median_ex_smoker", "median_smoker",
                     "median_difference_NS_VS_S", "fc_median__NS_VS_S", "pvalue_NS_VS_S",
                     "median_difference_NS_VS_EX", "fc_median__NS_VS_EX", "pvalue_NS_VS_EX",
                     "median_difference_EX_VS_S", "fc_median__EX_VS_S", "pvalue_EX_VS_S")

  res %>%
    as_tibble() 
}

# Parse and format the data
theta.longer <- reshape2::melt(theta)
colnames(theta.longer) <- c("Sample", "Celltype", "value")

# Parse and format the metadata
annot <- metadata_lung %>% 
  as_tibble()


data <- theta.longer %>% 
  merge(annot, by = "Sample") %>% 
  mutate(Smoking = as.factor(Smoking)) %>%
  as_tibble()


design <- formula(value ~ Smoking + HardyScale + IschemicTime + Age + Ancestry + Sex + BMI + atelectasis+ emphysema + fibrosis + pneumonia)

DA <- da_beta_regression(data, celltype= "Macrophage",design)

saveRDS(DA, file = "../figure/data/bayesPrism.rds")