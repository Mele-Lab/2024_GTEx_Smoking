#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to study the common effects of smoking and other demographic traits
# @software version: R=4.2.2

Sys.time()
#Set path 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")

#Loading functions to do enrichment:
gene_annotation <- read.csv("data/public/gene_annotation.csv")

library(WebGestaltR)
library(clusterProfiler)
library(org.Hs.eg.db)

# Reading data
tissues <- list.dirs("tissues/", full.names = F)[-1]
tissue_info <- read.csv("data/public/tissue_abreviation.txt")

tissues <- c("ArteryAorta", "Pancreas", "SkinSunExposedLowerleg", "Lung", "EsophagusMucosa")
results <- list()
for(tissue in tissues){
  print(tissue)
  dea <- readRDS(paste0("tissues/", tissue, "/voom_limma_results.rds"))
  smoking <- dea$Smoking2[dea$Smoking2$adj.P.Val<0.05,]
  age <- dea$Age[dea$Age$adj.P.Val<0.05,]
  up <- rownames(smoking)[smoking$logFC>0][rownames(smoking)[smoking$logFC>0] %in% rownames(age)[age$logFC>0]]
  down <- rownames(smoking)[smoking$logFC<0][rownames(smoking)[smoking$logFC<0] %in% rownames(age)[age$logFC<0]]
  # common <- c(up, down)
  
  bg <- union(rownames(smoking), rownames(age))
  
  #Do functional enrichments
  up <- gsub("\\..*", "", up)
  down <- gsub("\\..*", "", down)
  bg <- gsub("\\..*", "", bg)
  
  test <- enrichGO(gene= up, universe = bg,
           keyType = "ENSEMBL", OrgDb = org.Hs.eg.db,
           ont = "BP")
  
  if(nrow(test@result[test@result$p.adjust<0.05,])>0){
    results[[paste0(tissue, "_up")]] <- cbind(test@result[test@result$p.adjust<0.05,], Direction="up", Tissue=tissue)
  }
  
  test_down <- enrichGO(gene= down, universe = bg,
                        keyType = "ENSEMBL", OrgDb = org.Hs.eg.db,
                        ont = "BP")
  if(nrow(test_down@result[test_down@result$p.adjust<0.05,])>0){
    results[[paste0(tissue, "_down")]] <- cbind(test_down@result[test_down@result$p.adjust<0.05,], Direction="down", Tissue=tissue)
  }

}
to_share <- do.call(rbind, results)
to_share <- to_share[,c("Tissue", "Direction", "ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust")]
rownames(to_share) <- NULL


library("xlsx")
write.xlsx(to_share, "output/Supplementary_table_X.xlsx", 
           col.names = TRUE, append = FALSE, row.names = F)
