#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to correlate gene expression and DNA methylation
# @software version: R=4.2.2

#Loading libraries
library(dplyr)

#Set path 
setwd(system("pwd", intern = T)) #If in linux
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
# setwd("..")

gene_probe <- readRDS("output/gene_probe_pairs.rds")
gene_probe <- distinct(gene_probe)

#Get only DMPs and DEGs:
print("Running lung")
results_DML <- readRDS("tissues/Lung/DML_results.rds")
smoking2 <- results_DML$Smoking2
smoking2 <- rownames(smoking2[smoking2$adj.P.Val<0.05,])

dim(gene_probe)
gene_probe <- gene_probe[gene_probe$Name %in% smoking2,]
dim(gene_probe)

degs <-  readRDS("tissues/Lung/voom_limma_results.rds")
degs <- rownames(degs$Smoking2)[degs$Smoking2$adj.P.Val<0.05]

#From ensembl id to gene symbol
gene_annotation <- read.delim("data/public/gencode.v26.GRCh38.genes.bed", header=F)[,c(6,7)]
colnames(gene_annotation) <- c("gene", "symbol")
#These were genes duplicated, I changed their names to their correct one
gene_annotation$symbol[gene_annotation$gene=="ENSG00000253972.5"] <- "MAL2-AS1" #I found this a posteriori
gene_annotation$symbol[gene_annotation$gene=="ENSG00000283992.1"] <- "SLURP2" #Insted of LYNX1
gene_annotation$symbol[gene_annotation$gene=="ENSG00000235271.5"] <- "GC22P027299"
gene_annotation$symbol[gene_annotation$gene=="ENSG00000229694.6"] <- "C9orf73"
gene_annotation$symbol[gene_annotation$gene=="ENSG00000228741.2"] <- "GC13P024553" 

degs <- sapply(degs, function(gene) gene_annotation$symbol[gene_annotation$gene==gene])

# table(gene_probe$UCSC_RefGene_Name %in% degs) #10966/(10966+63954)*100

dim(gene_probe)
gene_probe <- gene_probe[gene_probe$UCSC_RefGene_Name %in% degs,]
dim(gene_probe)



#Reading methylationresiduals
beta <- readRDS(paste0("tissues/Lung/methylation_residuals.rds"))
Sys.time() 

#Reading expression residuals in the lung:
expression <- readRDS(paste0("tissues/Lung/expression_residuals_demographic_traits_no_smoking.rds"))

rownames(expression) <- sapply(rownames(expression), function(gene) gene_annotation$symbol[gene_annotation$gene==gene])
#From sample id to donor id
colnames(expression) <- sapply(colnames(expression), function(id) paste0(strsplit(id, "-")[[1]][-3], collapse="-"))
#Subset expression data to match the donors in DNA methylation data
expression <- expression[,colnames(expression) %in% colnames(beta)]
beta <- beta[,colnames(beta) %in% colnames(expression)] #I guess we have beta values for donors with no expression data because I only keep expression data for donors with smoking annotation
identical(colnames(beta), colnames(expression)) #same order

  
Sys.time() 
print("Starting to compute correlations")

#Compute correlations

#Now correlate beta values for each position and expression for its associated gene and compute Pearson's correlation
correlation_function <- function(gene, probe){
  if(!gene %in% rownames(expression)){ return(NA) }
  exp <- expression[rownames(expression)==gene,]
  if(!probe %in% rownames(beta)){ return(NA) }
  bet <- beta[rownames(beta)==probe,]
  test <- cor.test(exp, as.numeric(bet))
  to_return <- list("p.val"= test$p.value, "cor"=test$estimate, "gene"=gene, "probe"=probe)
  return(to_return)
}


output <- mapply(function(x,y) correlation_function(x, y), gene_probe[,2], gene_probe[,1] )
#Cleaning output
output_2 <- as.data.frame(t(output))
output_2$cor <- unlist(output_2$cor)
output_2$p.val <- unlist(output_2$p.val)
output_2$gene <- unlist(output_2$gene)
output_2$probe <- unlist(output_2$probe)
output_2 <- output_2[!is.na(output_2$p.val),]   #Removing NAs
output_2$p.adj <- p.adjust(output_2$p.val, method = "BH")

print("Finished computing correlations for lung")
Sys.time()

saveRDS(output_2, paste0("tissues/Lung/Correlations_DMP_DEG_new.rds"))



print("Running correlations for colon")
gene_probe <- readRDS("output/gene_probe_pairs.rds")
gene_probe <- distinct(gene_probe)

results_DML <- readRDS("tissues/ColonTransverse/DML_results.rds")
smoking2 <- results_DML$Smoking2
smoking2 <- rownames(smoking2[smoking2$adj.P.Val<0.05,])

dim(gene_probe)
gene_probe <- gene_probe[gene_probe$Name %in% smoking2,]
dim(gene_probe)

#Keeping only DEGs
degs <-  readRDS("tissues/ColonTransverse/voom_limma_results.rds")
degs <- rownames(degs$Smoking2)[degs$Smoking2$adj.P.Val<0.05]
degs <- sapply(degs, function(gene) gene_annotation$symbol[gene_annotation$gene==gene])

dim(gene_probe)
gene_probe <- gene_probe[gene_probe$UCSC_RefGene_Name %in% degs,]
dim(gene_probe)

#Reading methylation residuals
beta <- readRDS(paste0("tissues/ColonTransverse/methylation_residuals.rds"))
#Reading expression residuals in the lung:
expression <- readRDS(paste0("tissues/ColonTransverse/expression_residuals_demographic_traits_no_smoking.rds"))

rownames(expression) <- sapply(rownames(expression), function(gene) gene_annotation$symbol[gene_annotation$gene==gene])
#From sample id to donor id
colnames(expression) <- sapply(colnames(expression), function(id) paste0(strsplit(id, "-")[[1]][-3], collapse="-"))
#Subset expression data to match the donors in DNA methylation data
expression <- expression[,colnames(expression) %in% colnames(beta)]
beta <- beta[,colnames(beta) %in% colnames(expression)] #I guess we have beta values for donors with no expression data because I only keep expression data for donors with smoking annotation
identical(colnames(beta), colnames(expression)) #same order

Sys.time() 
print("Starting to compute correlations")
output <- mapply(function(x,y) correlation_function(x, y), gene_probe[,2], gene_probe[,1] )
#Cleaning output
output_2 <- as.data.frame(t(output))
output_2$cor <- unlist(output_2$cor)
output_2$p.val <- unlist(output_2$p.val)
output_2$gene <- unlist(output_2$gene)
output_2$probe <- unlist(output_2$probe)
output_2 <- output_2[!is.na(output_2$p.val),]   #Removing NAs
output_2$p.adj <- p.adjust(output_2$p.val, method = "BH")

print("Finished computing correlations for colon")
Sys.time()

saveRDS(output_2, paste0("tissues/ColonTransverse/Correlations_DMP_DEG_new.rds"))
Sys.time() 
