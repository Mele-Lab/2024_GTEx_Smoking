#!/usr/bin/env Rscript
# @Author: Rogério Eduardo Ramos Ribeiro
# @E-mail: zroger499@gmail.com
# @Description: Classify each gene in terms of reversability (R reversible, NR non reversible and PR Partial reversible)
# @software version: R=4.2.2

# Load libraries 
library(tidyverse)

#Set working dir 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load data
degs_analysis <- readRDS(file = "output/differential_expression_results.rds")

## Get differential expressed genes ----
get_degs <- function(DEGS, comparison, fdr  = 0.05){
  res <- c()
  for (tissue in names(DEGS)){
    degs <- DEGS[[tissue]][[comparison]] %>% 
      filter(adj.P.Val < fdr) %>% 
      mutate(ensembl = gsub("\\.\\d+", "", row.names(.))) %>% 
      mutate(tissue = tissue)
    res <- rbind(res, degs)
  }
  res
}


NS_S_DEGS <- get_degs(degs_analysis, comparison ="SmokingSMOKER-NEVER")
EX_S_DEGS <- get_degs(degs_analysis, comparison ="SmokingEX-SMOKER")
NS_EX_DEGS <- get_degs(degs_analysis, comparison = "SmokingEX-NEVER")

## Get reversible genes and non reversible genes per tissue ----

### reversible_genes - DEGS in NS vs S and EX vs S

NS_S_DEGS <- NS_S_DEGS %>% 
  mutate(direction = ifelse(logFC > 0, "+", "-")) %>%
  mutate(name = paste0(ensembl, "_", direction, "_", tissue))

colnames(NS_S_DEGS) <- paste0(colnames(NS_S_DEGS), "_", "NS_VS_S")

EX_S_DEGS <- EX_S_DEGS %>% 
  mutate(direction = ifelse(logFC > 0, "-", "+")) %>% # logFC here is upregulated in EX
  mutate(name = paste0(ensembl, "_", direction, "_", tissue))

colnames(EX_S_DEGS) <- paste0(colnames(EX_S_DEGS), "_", "EX_VS_S")


reversible_genes <- merge(NS_S_DEGS, EX_S_DEGS, by.x = "name_NS_VS_S", by.y = "name_EX_VS_S")

### Non Reversible genes ----

NS_EX_DEGS <- NS_EX_DEGS %>% 
  mutate(direction = ifelse(logFC > 0, "+", "-")) %>% 
  mutate(name = paste0(ensembl, "_", direction, "_", tissue))

colnames(NS_EX_DEGS) <- paste0(colnames(NS_EX_DEGS), "_", "NS_VS_EX")

non_reversible_genes <- merge(NS_S_DEGS, NS_EX_DEGS, by.x = "name_NS_VS_S", by.y = "name_NS_VS_EX")

### Summarise the information in a table ----

#Some genes are classified as both classified as reversible and never reversible
clear_pr_genes <- intersect(gsub("_[+|-]_", "", non_reversible_genes$name_NS_VS_S), gsub("_[+|-]_", "", reversible_genes$name_NS_VS_S))

#Annotate the complete gene table

DEGS <- get_degs(degs_analysis, comparison ="SmokingSMOKER-NEVER") %>% mutate(name = paste0(ensembl, tissue))
DEGS$reversability <- NA

#Add reversible genes
DEGS$reversability <- ifelse(DEGS$name %in% gsub("_[+|-]_", "", reversible_genes$name_NS_VS_S), "reversible", NA)
#Add non reversible genes 
DEGS$reversability <- ifelse(DEGS$name %in% gsub("_[+|-]_", "", non_reversible_genes$name_NS_VS_S), "non-reversible", DEGS$reversability )
#Add partial reversible genes 
DEGS$reversability <- ifelse(DEGS$name %in% clear_pr_genes, "partially reversible", DEGS$reversability)
DEGS$reversability <- ifelse(is.na(DEGS$reversability), "partially reversible", DEGS$reversability)


### Get the paired logFC between genes in the NS vs EX and EX vs S (For PR genes) ----
pr_genes <- DEGS %>% 
  filter(reversability == "partially reversible") %>% 
  pull(name)


NS_EX_pr <- get_degs(degs_analysis, comparison ="SmokingEX-NEVER", fdr = 1) %>%
  rename(logFC_NS_EX = logFC) %>%
  mutate(gene_tissue = paste0(ensembl, tissue)) %>% 
  filter(gene_tissue %in% pr_genes)


EX_S_pr <- get_degs(degs_analysis, comparison = "SmokingEX-SMOKER", fdr = 1) %>% 
  rename(logFC_EX_VS_S = logFC) %>%
  mutate(gene_tissue = paste0(ensembl, tissue)) %>%
  filter(gene_tissue %in% pr_genes)


EX_comparisons_logFC <- merge(NS_EX_pr, EX_S_pr, by = "gene_tissue") %>% 
  select(logFC_NS_EX, logFC_EX_VS_S, gene_name.x, tissue.x) %>% 
  rename(gene_name = gene_name.x) %>%
  rename(tissue = tissue.x)

## Percent of PR genes close to Never smokers 
sum(abs(EX_comparisons_logFC$logFC_NS_EX) < abs(EX_comparisons_logFC$logFC_EX_VS_S)) / nrow(EX_comparisons_logFC)

#Wilcoxon test between logFC
test <- wilcox.test(abs(EX_comparisons_logFC$logFC_NS_EX), abs(EX_comparisons_logFC$logFC_EX_VS_S), paired = T)
test$p.value

median(abs(EX_comparisons_logFC$logFC_NS_EX)) - median(abs(EX_comparisons_logFC$logFC_EX_VS_S))

## What is the precent of genes with higher logFC in S vs EX?
100 *sum(abs(EX_comparisons_logFC$logFC_EX_VS_S) > abs(EX_comparisons_logFC$logFC_NS_EX)) / nrow(EX_comparisons_logFC)

### Percentage of reversible genes in each tissue 
percent_reverse_status <- c()
for (tiss in names(table(DEGS$tissue))){
  for (status in c("reversible", "partially reversible", "non-reversible")){
    percent_group <-  nrow(DEGS %>% filter(reversability == status & tissue == tiss)) /  nrow(DEGS %>% filter(tissue == tiss))
    percent_reverse_status <- rbind(percent_reverse_status, c(tiss, status, percent_group*100))
  }
}
colnames(percent_reverse_status) <- c("tissue", "status", "percent")
percent_reverse_status <- data.frame(percent_reverse_status) %>% mutate(percent = as.numeric(percent))

DEGS.per.tissue <- DEGS %>% 
  group_by(tissue) %>% 
  summarise(n = n())

percent_reverse_status <- merge(percent_reverse_status, DEGS.per.tissue, by = "tissue")


## Overall percentage of reversible and non-reversible gene 
100 * sum(DEGS$reversability == "reversible") / nrow(DEGS)
100 * sum(DEGS$reversability == "non-reversible") / nrow(DEGS)
100 * sum(DEGS$reversability == "partially reversible") / nrow(DEGS)


### Examples of Reversible, PR and non reversible genes

metadata  <- readRDS(file = "output/metadata/metadata.rds")
expression <- fread(file = "public/expression/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz") #Expression table in TPM


get_gene_expression <- function(exp, gene, metadata, tissue){
  #Get data to plot expression levels of genes (logs +1 the TPM)
  
  #Get samples names from metadata 
  samples.smk <- metadata[[tissue]] %>% filter(Smoking == 2) %>% pull(Sample)
  samples.ns <- metadata[[tissue]] %>% filter(Smoking == 0) %>% pull(Sample)
  samples.ex <- metadata[[tissue]] %>% filter(Smoking == 1) %>% pull(Sample)
  
  #Modify sample names in the expression dataset
  colnames(exp) <- gsub("\\.SM.*", "", colnames(exp))
  colnames(exp) <- gsub("\\.", "-", colnames(exp))
  exp$Name <- gsub("\\.\\d+", "", exp$Name)
  
  exp_gene.smk <- log(as.numeric(exp %>% filter(Name == gene) %>% select(all_of(samples.smk))) + 1)
  exp_gene.Nonsmk <- log(as.numeric(exp %>% filter(Name == gene) %>% select(all_of(samples.ns))) + 1) 
  exp_gene.exsmk <- log(as.numeric(exp %>% filter(Name == gene) %>% select(all_of(samples.ex))) + 1)
  
  tpm <- c(exp_gene.smk, exp_gene.Nonsmk, exp_gene.exsmk)
  label <- c(rep("Smoker", length(exp_gene.smk)), rep("Never Smoker", length(exp_gene.Nonsmk)), rep("Ex Smoker", length(exp_gene.exsmk)))
  
  expression_data <- data.frame(exp = tpm, label = label)
  
  expression_data %>% 
    mutate(label = factor(label, levels = c("Never Smoker", "Ex Smoker", "Smoker"))) %>% 
    mutate(gene_name = gene) %>% 
    mutate(tissue = tissue)
}


# Example of genes
reversible_gene_example <- get_gene_expression(expression, "ENSG00000140465", metadata, "Lung") #	ENSG00000140465 (CYP1A1)
pr_gene_example <- get_gene_expression(expression, "ENSG00000154165", metadata, "Lung") #ENSG00000154165 (GPR15)
non_reversible_genes_example <-  get_gene_expression(expression, "ENSG00000105696", metadata, "Thyroid") #ENSG00000135898 (TMEM59L)


# Save the information 
output <- list()
output$reversability <- DEGS
output$logFC <- EX_comparisons_logFC
output$reversible_example <- reversible_gene_example
output$pr_gene_example <- pr_gene_example
output$non_reversible_genes_example <- non_reversible_genes_example

saveRDS(output, file = "../figures/data/reversible_genes.rds")