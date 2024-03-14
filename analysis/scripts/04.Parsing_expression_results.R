#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to parse expression results to share in Zenodo
# @software version: R=4.2.2

Sys.time()
#Set path 
setwd(system("pwd", intern = T)) #If in linux
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
# setwd("..")

# Tissues ----
tissues <- list.dirs("tissues/", full.names = F)[-1]
tissues <- tissues[tissues!="KidneyCortex"]
sex_tissues <- c("Uterus", "Vagina", "Ovary", "Prostate", "Testis")

# Differential expression analyses: results tables ----
dea_res <- lapply(tissues, function(tissue) readRDS(paste0("tissues/", tissue, "/voom_limma_results.rds")))
names(dea_res) <- tissues

#Residuals just to share:
residuals <- lapply(tissues, function(tissue) readRDS(paste0("tissues/", tissue, "/expression_residuals_demographic_traits_no_smoking.rds")))
names(residuals) <- tissues
saveRDS(residuals, "output/residuals_to_share.rds")

# Hier.part ----
hier_part <- lapply(tissues, function(tissue) readRDS(paste0("tissues/", tissue, "/hier.part.rds"))) #Residuals when correcting for all covariates but not Age, Ancestry, BMI, Sex, Smoking
names(hier_part) <- tissues

print("finished reading data")

for(tissue in tissues){
  print(tissue)
  #Renaming variables to a more clear name
  names(dea_res[[tissue]])[names(dea_res[[tissue]])=="AncestryAMR"] <- "AncestryAMR-AFR"
  names(dea_res[[tissue]])[names(dea_res[[tissue]])=="AncestryEUR"] <- "AncestryEUR-AFR"
  names(dea_res[[tissue]])[names(dea_res[[tissue]])=="AncestryAMR-AncestryEUR"] <- "AncestryAMR-EUR"
  names(dea_res[[tissue]])[names(dea_res[[tissue]])=="Sex2"] <- "Sex" #The reference level was male
  names(dea_res[[tissue]])[names(dea_res[[tissue]])=="Smoking1"] <- "SmokingEX-NEVER"
  names(dea_res[[tissue]])[names(dea_res[[tissue]])=="Smoking2"] <- "SmokingSMOKER-NEVER"
  names(dea_res[[tissue]])[names(dea_res[[tissue]])=="Smoking1-Smoking2"] <- "SmokingEX-SMOKER"

  for(trait in names(dea_res[[tissue]])){
    if(grepl("Ancestry", trait)){
      trait_name <- "Ancestry"
    } else if(grepl("Smoking", trait)){
      trait_name <- "Smoking"
    } else{
      trait_name <- trait
    }
    if(tissue %in% sex_tissues & trait == "Sex"){
      next
    }else{
      dea_res[[tissue]][[trait]]$R2 <- NA
      for(gene in rownames(dea_res[[tissue]][[trait]])){
        if(dea_res[[tissue]][[trait]][gene, "adj.P.Val"] < 0.05){
          dea_res[[tissue]][[trait]][gene, "R2"] <- hier_part[[tissue]][gene, paste0(trait_name, "_abs")]
        }
      }
    }
  }
}

Sys.time()
saveRDS(dea_res, "output/differential_expression_results_uncomplete.rds")
# dea_res <- readRDS("output/differential_expression_results_uncomplete.rds")

#Add gene_name
print("Adding gene name")
gene_annotation <- read.delim("data/public/gencode.v26.GRCh38.genes.bed", header=F)[,c(6,7)]
colnames(gene_annotation) <- c("gene", "symbol")

for(tissue in tissues){
  print(tissue)
  Sys.time()
  for(trait in names(dea_res[[tissue]])){
    if(tissue %in% sex_tissues & trait == "Sex"){
      print(paste0(tissue, ": ", trait))
      next
    }else{
      dea_res[[tissue]][[trait]]$gene_name <- sapply(rownames(dea_res[[tissue]][[trait]]), function(ensembl) gene_annotation$symbol[gene_annotation$gene==ensembl])
      dea_res[[tissue]][[trait]] <- dea_res[[tissue]][[trait]][, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "gene_name", "R2")] #Changing order
    }
  }
}
print("saving")
Sys.time()

saveRDS(dea_res, "output/differential_expression_results.rds")

#Parsing metadata as well:
metadata <- lapply(tissues, function(tissue) readRDS(paste0("tissues/", tissue, "/metadata.rds"))) #Residuals when correcting for all covariates but not Age, Ancestry, BMI, Sex, Smoking
names(metadata) <- tissues
saveRDS(metadata, file="output/metadata.rds")
