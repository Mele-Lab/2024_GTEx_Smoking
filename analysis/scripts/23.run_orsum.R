#! /usr/bin/env Rscript
# @Author: Rogério Eduardo Ramos Ribeiro
# @E-mail: rogerio.e.ramos.ribeiro
# @Description: Plot the GO terms from MOFA enrichement analysis with permutations
# @software version: R=4.4.2

library(MOFAdata)
library(DOSE)
library(GO.db)
library(tidyverse)

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888", 
                             "black")

out.dir <- "MOFA/terms_per_view_lung"
if (!dir.exists(out.dir)) {dir.create(out.dir)}

res_enrich <- readRDS(file = "MOFA/enrichement_MOFA.rds")

names(res_enrich$lung) <- c("Gene Expression", "Promoters", "Enhancers", "Gene Body")
factor_n <- 2 
enriched_terms <- data.frame()
for (view in names(res_enrich$lung)){
  for (direction in c("positive", "negative")){
    write.table(res_enrich$lung[[view]][[direction]]$sigPathways[[factor_n]], paste0(out.dir, "/", direction, "_GO_terms_", str_replace(view, " ", ""), ".csv"), quote = F, row.names = F, col.names = F)
    enriched_terms <- rbind(enriched_terms, 
                            data.frame("view" = view, "direction" = direction, "terms" = res_enrich$lung[[view]][[direction]]$sigPathways[[factor_n]]))
  }
}

# Run ORSUM
#orsum.py --gmt  public/gmt/hsapiens.GO_BP.name.gmt --files "MOFA/terms_per_view_lung/positive_GO_terms_Enhancers.csv" "MOFA/terms_per_view_lung/positive_GO_terms_GeneBody.csv" "MOFA/terms_per_view_lung/positive_GO_terms_GeneExpression.csv" "MOFA/terms_per_view_lung/positive_GO_terms_Promoters.csv"  "MOFA/terms_per_view_lung/negative_GO_terms_Enhancers.csv" "MOFA/terms_per_view_lung/negative_GO_terms_GeneBody.csv" "MOFA/terms_per_view_lung/negative_GO_terms_GeneExpression.csv" "MOFA/terms_per_view_lung/negative_GO_terms_Promoters.csv" --fileAliases Enhancers_positive GeneBody_positive GeneExpression_positive Promoters_positive Enhancers_negative GeneBody_negative GeneExpression_negative Promoters_negative --outputFolder MOFA/ORSUM/all --minTermSize 1
