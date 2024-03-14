#!/usr/bin/env Rscript
# @Author: Rogério Eduardo Ramos Ribeiro
# @E-mail: zroger499@gmail.com
# @Description: Code to run enrichement analysis on a per tissue basis for up and down genes in seperate
# @software version: R=4.2.2


# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
suppressPackageStartupMessages(library(tidyverse))

#Set work dir 
setwd(system("pwd", intern = T)) #If in linux
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio

#Load degs 
#DEGS <- readRDS(file = "../../DGE/differential_expression_results.rds")
DEGS <- readRDS(file ="output/differential_expression_results.rds")

#Create ouput folder 
#out.folder = "../../enrichement/"
out.folder = "output/DEGS_enrichement/"

if (!dir.exists(out.folder)){
  dir.create(out.folder, recursive = T)
}


#Main analysis function
perform_enrichment <- function(data, comparison){
  res <- c()  
  for (tissue in names(data)){ 
    tissue_data <- data[[tissue]][[comparison]]
    deg_tissue_list <- tissue_data %>% filter(adj.P.Val < 0.05)
    deg_tissue_list.up <- gsub("\\.\\d+", "", row.names(tissue_data %>% filter(adj.P.Val < 0.05 & logFC > 0)))
    deg_tissue_list.down <- gsub("\\.\\d+", "", row.names(tissue_data %>% filter(adj.P.Val < 0.05 & logFC < 0)))

    out <- tryCatch({
      deg_tissue_list_entrez.up <- bitr(deg_tissue_list.up, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)["ENTREZID"]
    },error=function(cond) {
      deg_tissue_list_entrez.up <- c()
    })
    
    out <- tryCatch({
      deg_tissue_list_entrez.down <- bitr(deg_tissue_list.down, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)["ENTREZID"]
    },error=function(cond) {
      deg_tissue_list_entrez.down <- c()
    })
    
    
    background <- gsub("\\.\\d+", "", row.names(tissue_data))
    background_entrez <- bitr(background, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)["ENTREZID"]
    
    # Enrichment analysis
    if (length(deg_tissue_list.up) > 5){
      go.bp.up <- enrichGO(gene = deg_tissue_list.up, keyType = "ENSEMBL", OrgDb='org.Hs.eg.db', pvalueCutoff = 0.05, ont = "BP", pAdjustMethod = "BH", universe = background)
      go.mf.up <- enrichGO(gene = deg_tissue_list.up, keyType = "ENSEMBL", OrgDb='org.Hs.eg.db', pvalueCutoff = 0.05, ont = "MF", pAdjustMethod = "BH", universe = background)
      go.cc.up <- enrichGO(gene = deg_tissue_list.up, keyType = "ENSEMBL", OrgDb='org.Hs.eg.db', pvalueCutoff = 0.05, ont = "CC", pAdjustMethod = "BH", universe = background)
      res[[tissue]][["BP_up"]] <- go.bp.up
      res[[tissue]][["MF_up"]] <- go.mf.up
      res[[tissue]][["CC_up"]] <- go.cc.up
    }
    
    if (length(deg_tissue_list.down) > 5){
      go.bp.down <- enrichGO(gene = deg_tissue_list.down, keyType = "ENSEMBL", OrgDb='org.Hs.eg.db', pvalueCutoff = 0.05, ont = "BP", pAdjustMethod = "BH", universe = background)
      go.mf.down <- enrichGO(gene = deg_tissue_list.down, keyType = "ENSEMBL", OrgDb='org.Hs.eg.db', pvalueCutoff = 0.05, ont = "MF", pAdjustMethod = "BH", universe = background)
      go.cc.down <- enrichGO(gene = deg_tissue_list.down, keyType = "ENSEMBL", OrgDb='org.Hs.eg.db', pvalueCutoff = 0.05, ont = "CC", pAdjustMethod = "BH", universe = background)
      res[[tissue]][["BP_down"]] <- go.bp.down
      res[[tissue]][["MF_down"]] <- go.mf.down
      res[[tissue]][["CC_down"]] <- go.cc.down
    }
    
    if (length(deg_tissue_list_entrez.up$ENTREZID) > 5){
      kegg.up <- enrichKEGG(gene = deg_tissue_list_entrez.up$ENTREZID,  organism = "hsa", pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_entrez$ENTREZID)
      do.up <- enrichDO(gene = deg_tissue_list_entrez.up$ENTREZID, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_entrez$ENTREZID)
      dgn.up <- enrichDGN(gene = deg_tissue_list_entrez.up$ENTREZID, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_entrez$ENTREZID)
      res[[tissue]][["kegg_up"]] <- kegg.up
      res[[tissue]][["do_up"]] <- do.up
      res[[tissue]][["dng_up"]] <- dgn.up
    }
    
    if (length(deg_tissue_list_entrez.down$ENTREZID) > 5){
      kegg.down <- enrichKEGG(gene = deg_tissue_list_entrez.down$ENTREZID,  organism = "hsa", pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_entrez$ENTREZID)
      do.down <- enrichDO(gene = deg_tissue_list_entrez.down$ENTREZID, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_entrez$ENTREZID)
      dgn.down <- enrichDGN(gene = deg_tissue_list_entrez.down$ENTREZID, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_entrez$ENTREZID)
      
      res[[tissue]][["kegg_down"]] <- kegg.down
      res[[tissue]][["do_down"]] <- do.down
      res[[tissue]][["dng_down"]] <- dgn.down
      
    }
    
  }
  
  return (res)
}


############### ------------------- ####################

#"SmokingEX-NEVER" is the comparison between 0 and 1 (never vs ex-smokers), 
#"SmokingSMOKER-NEVER" is the comparison between 0 and 2 (never vs smokers) and 
#"SmokingEX-SMOKER" is the comparison between 1 and 2 (ex-smokers vs smokers). 

#########################################################


#Perform enrichment analysis for each pairwise comparison
enrichement_NS_vs_S <- perform_enrichment(DEGS, "SmokingSMOKER-NEVER")
#enrichement_EX_vs_S <- perform_enrichment(DEGS, "SmokingEX-SMOKER")
#enrichement_NS_vs_EX <- perform_enrichment(DEGS, "SmokingEX-NEVER")


# Save results 
saveRDS(enrichement_NS_vs_S, file = paste0(out.folder, "never_vs_smoker_enrichement.rds"))
#saveRDS(enrichement_EX_vs_S, paste0(out.folder, "ex_vs_smoker_enrichement.rds"))
#saveRDS(enrichement_NS_vs_EX, paste0(out.folder, "never_vs_ex_enrichement.rds"))



### Prepare data to run ORSUM (Never vs Smoker)
save_enriched_terms <- function(enrichment, slot, FDR, minCount, folder){
  if (!dir.exists(folder)){
    dir.create(folder, recursive = T)
  }
  for (tissue in names(enrichment)){
    if (!is.null(enrichment[[tissue]]) & !is.null(enrichment[[tissue]][[slot]])){
      sig <- enrichment[[tissue]][[slot]]@result %>% filter(p.adjust < FDR & Count > minCount) %>% pull(ID)
      if (length(sig) > 0){
        filename <- paste0(tissue, "_", slot)
        full_filename <- paste0(folder, "/", filename, ".csv")
        write.table(sig, file = full_filename, quote = F, row.names = F, col.names = F)
      }
    }
  }
}



save_enriched_terms(enrichement_NS_vs_S, "BP_up", 0.05, 5, paste0(out.folder, "enriched_terms.go.bp.upregulated"))
save_enriched_terms(enrichement_NS_vs_S, "BP_down", 0.05, 5, paste0(out.folder, "enriched_terms.go.bp.downregulated"))


# Generate script to run ORSUM
generate_orsum_command <- function(dir, ontology_file, ontology, direction, server_dir_files ,additional_options = ""){
  command <- "orsum.py --gmt "
  command <- paste(command, ontology_file)
  file_list <- list.files(dir)
  file_list_filtered <- file_list[grepl(ontology, file_list)]
  
  
  command <- paste0(command, " --files")
  for (i in 1:length(file_list_filtered)){
    command <- paste0(command, ' \'', dir, "/", file_list_filtered[i],'\'')
  }
  
  
  command <- paste0(command, " --fileAliases ")
  for (i in 1:length(file_list_filtered)){
    tissue <-  str_split(file_list_filtered[i], "_")[[1]][1]
    command <- paste(command, tissue, sep = " ")
  }
  
  command = paste0(command, " ", additional_options)
  return (command)
}


gmt_file <- "data/gmt/hsapiens.GO_BP.name.gmt"
#output_orsum = "../../enrichement/"  
output_orsum <- "output/DEGS_enrichement/"



command_up.bp <- generate_orsum_command(dir = paste0(out.folder, "enriched_terms.go.bp.upregulated"), 
                                        gmt_file, "BP", "up", 
                                        additional_options = paste0("--outputFolder ", output_orsum, "BP.up", " --minTermSize 1"))


command_up.down <-  generate_orsum_command(dir = paste0(out.folder, "enriched_terms.go.bp.downregulated"), 
                                           gmt_file, "BP", "down", 
                                           additional_options = paste0("--outputFolder ", output_orsum, "BP.down", " --minTermSize 1"))


#Write commands to a file
fileConn<-file("scripts/05_run_orsum.sh")

writeLines(c("#!/bin/bash", "\n", 
             paste0("mkdir -p ", output_orsum, "BP.up"), 
             paste0("mkdir -p ", output_orsum, "BP.down"),
             "\n", "\n",
             command_up.bp,command_up.down,
             "orsum.py -v > orsum.version"), 
             fileConn)

close(fileConn)


# Write session info in R 
writeLines(capture.output(sessionInfo()), "enrichement_sessionInfo.txt")
