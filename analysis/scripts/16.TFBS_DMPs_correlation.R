#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Do DMPs in TFBS correlate more with expression changes?
# @software version: R=4.2.2

#Set path 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")

data_path <- "~/Documents/mn4/scratch/bsc83/bsc83535/GTEx/v9/Oliva/"
TFBS_cpgs <- read.csv(paste0(data_path, "chip_seq_processed.csv")) #Lung specific
colnames(TFBS_cpgs)[8] <- c("TF")
TFBS_cpgs <- TFBS_cpgs[,c(5,8)]
library(dplyr)
TFBS_cpgs <- TFBS_cpgs %>% distinct()

cor <- readRDS("tissues/Lung/Correlations_DMP_DEG_new.rds")
cor$TF <- "no"
cor$TF[cor$probe %in% TFBS_cpgs$name_ann] <- "yes"

results_DML <- readRDS(paste0("tissues/Lung/DML_results.rds"))
smoking <- results_DML$Smoking2
smoking_hyper <- rownames(smoking[smoking$adj.P.Val<0.05 & smoking$logFC>0,])
smoking_hypo <- rownames(smoking[smoking$adj.P.Val<0.05 & smoking$logFC<0,])


#Fisher's
# signif <- cor[cor$p.adj<0.05,]
cor$signif <- cor$p.adj<0.05
# TFBS_cpgs <- TFBS_cpgs[TFBS_cpgs$name_ann %in% cor$probe,]  #background is our set of pairs
cor$signif <- factor(cor$signif, levels=c(T, F))
cor$TF <- factor(cor$TF, levels=c("yes", "no"))

hypo <- cor[cor$probe %in% smoking_hypo,] 
hyper <- cor[cor$probe %in% smoking_hyper,] 


#The DMPs associated with DEGs that are in a TFBS are more likely for its methylation level to correlate with the expression level of the gene?
table(hypo$TF, hypo$signif)
fisher.test(hypo$TF, hypo$signif)

table(hyper$TF, hyper$signif)
fisher.test(hyper$TF, hyper$signif)



#Wilcoxon

# library(ggplot2)
# ggplot(hypo) + geom_boxplot(aes(TF, abs(cor)))
# t <- wilcox.test(abs(hypo$cor[hypo$TF=="no"]), abs(hypo$cor[hypo$TF=="yes"])) 
# wilcox.test(abs(hypo$cor[hypo$TF=="no"]), abs(hypo$cor[hypo$TF=="yes"]), alternative = "less") #x is not less than y
# #less means that the alternative hypothesis is that x is smaller than y, which is our hypothesis
# 
# ggplot(hyper) + geom_boxplot(aes(TF, abs(cor)))
# wilcox.test(abs(hyper$cor[hyper$TF=="no"]), abs(hyper$cor[hyper$TF=="yes"]), alternative = "less") #x is not less than y

#With wilcoxon of rho values in correlation we get the same results