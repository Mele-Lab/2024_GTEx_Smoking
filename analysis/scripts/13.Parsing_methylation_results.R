#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Methylation results analysis
# @software version: R=4.2.2

#Set path 
# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")


tissues <- c("BreastMammaryTissue", "ColonTransverse", "Lung", "MuscleSkeletal", "Ovary", "Prostate", "Testis", "WholeBlood")
tissues <- c("ColonTransverse", "Lung")

#Never vs Smokers:

table <- data.frame('CpG'='1', 'logFC'='1', 'AveExpr'='1', 't'='1', 'P.Value'='1', 'adj.P.Val'='1', 'B'='1', 'tissue'='1')
for(tissue in tissues){
  print(tissue)
  results <- readRDS(paste0("tissues/", tissue, "/DML_results.rds"))
  signif <- results$Smoking2[results$Smoking2$adj.P.Val<0.05,]
  if(nrow(signif)>0){
    table <- rbind(table, cbind('CpG'=rownames(signif), signif, tissue))
  }
}
table <- table[-1,]
table <- table[,c("CpG", "logFC", "adj.P.Val", "tissue")]
rownames(table) <- NULL

#Add reversability and save:
table_2 <- data.frame('CpG'='1', 'logFC'='1', 'AveExpr'='1', 't'='1', 'P.Value'='1', 'adj.P.Val'='1', 'B'='1', 'tissue'='1', 'Comparison'='1')
table$reversability <- "partially reversible"
for(tissue in unique(table$tissue)){
  print(tissue)
  results <- readRDS(paste0("tissues/", tissue, "/DML_results.rds"))
  subset <- table[table$tissue==tissue,] #DMPs
  ex_vs_smoker <- results$`Smoking1-Smoking2`[results$`Smoking1-Smoking2`$adj.P.Val<0.05,]
  if(nrow(ex_vs_smoker>0)){
    table_2 <- rbind(table_2, cbind('CpG'=rownames(ex_vs_smoker), ex_vs_smoker, tissue, "Comparison"="Ex vs Smoker"))
  }
  never_vs_ex <- results$Smoking1[results$Smoking1$adj.P.Val<0.05,]
  if(nrow(never_vs_ex>0)){
    table_2 <- rbind(table_2, cbind('CpG'=rownames(never_vs_ex), never_vs_ex, tissue, "Comparison"="Never vs Ex"))
  }
  rev <- subset$CpG[subset$CpG %in% rownames(ex_vs_smoker)]
  rev <- rev[!rev %in% rownames(never_vs_ex)]
  
  non_rev <-  subset$CpG[subset$CpG %in% rownames(never_vs_ex)]
  non_rev <- non_rev[!non_rev %in% rownames(ex_vs_smoker)]
  
  table$reversability[table$tissue==tissue & table$CpG %in% rev] <- "reversible"
  table$reversability[table$tissue==tissue & table$CpG %in% non_rev] <- "non-reversible"
}
table_2 <- table_2[-1,]
table_2 <- table_2[,c("CpG", "logFC", "adj.P.Val", "tissue", "Comparison")]
colnames(table_2)[colnames(table_2)=="tissue"] <- "Tissue"

library(writexl)
write_xlsx(table, "output/Supplementary_table_9.xlsx")
write_xlsx(table_2, "output/Supplementary_table_15.xlsx")
