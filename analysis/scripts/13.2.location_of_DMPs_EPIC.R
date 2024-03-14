#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Enrichment on genomic location for differentially methylated positions/loci
# @software version: R=4.2.2

library(dplyr)
library(tidyr)
library(tidyverse)

#Set path 
# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")

Sys.time()
data_path <- "~/Documents/mn4/scratch/bsc83/bsc83535/GTEx/v9/Oliva/"
# annotation <- read.csv(paste0(first_dir, "/scratch/bsc83/bsc83535/GTEx/v9/Oliva/GPL21145_MethylationEPIC_15073387_v-1-0.csv"))
# annotation <- annotation[7:nrow(annotation),] #removing header
# row_odd <- seq_len(nrow(annotation)) %% 2 #dummy variable for odd rows
# annotation_row_odd <- annotation[row_odd == 1, ]
# annotation_row_even <- annotation[row_odd == 0, ]
# annotation_row_odd <- annotation_row_odd[1:867927,] #There are empty rows at the bottom because of the controls
# annotation_row_even <- annotation_row_even[1:867927,] #There are extra information at the bottom because of the controls
# final_annotation <- cbind(annotation_row_odd, annotation_row_even)
# #There are 7 extra columns at the end with no column name
# colnames(final_annotation) <- final_annotation[1,]
# final_annotation <- final_annotation[-1,]

#Annotate nearest genes to promoter and enhancer probes with the gene missing. 
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# library(bumphunter)
# genes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
# to_annotate <- final_annotation[final_annotation$UCSC_RefGene_Name=="",]
# to_annotate <- to_annotate[to_annotate$Phantom5_Enhancers!="" | to_annotate$Regulatory_Feature_Group=="Promoter_Associated",]
# ## annotate promoter/enhancer probes that do not have a gene assigned
# to_annotate$CHR <- paste0('chr',to_annotate$CHR)
# to_annotate$start <- to_annotate$MAPINFO
# to_annotate$end <- as.numeric(to_annotate$MAPINFO) + 2
# rownames(to_annotate) <- to_annotate$IlmnID
# tab <- matchGenes(to_annotate,genes)
# to_annotate$UCSC_RefGene_Name <- tab$name
# ### merge data ----
# to_annotate$start <- NULL
# to_annotate$end <- NULL
# to_annotate$CHR <- gsub('chr','',to_annotate$CHR)
#
# final_annotation <- final_annotation[!final_annotation$Name %in% to_annotate$Name,]
# final_annotation <- rbind(final_annotation, to_annotate) #63945
# write.csv(final_annotation, paste0(first_dir, "/scratch/bsc83/bsc83535/GTEx/v9/Oliva/GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv"))

#Subsetting for the tested ones:
# tissue <- "Lung"
# res <- readRDS(paste0("tissues/", tissue , "/DML_results.rds"))[["Smoking2"]]
# annotation <- annotation[annotation$Name %in% rownames(res),]
# write.csv(annotation, paste0(data_path, "GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv"))

annotation <- read.csv(paste0(data_path, "GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv"))
Sys.time()

tissue <- "Lung"
# tissue <- "ColonTransverse"
res <- readRDS(paste0("tissues/", tissue , "/DML_results.rds"))[["Smoking2"]]

#We will divide all cpgs into promoter, enhancer, gene body or intergenic:

#We first classify as promoters the CpG annotated as "Promoter_Associated" and TSS. If a CpG is associated to two genes, one as promoter and another as gene body, we keep the CpG only as associated to the gene as promoter and exclude the other gene (very few cases) 
library(dplyr)
test <- annotation %>% separate_rows(UCSC_RefGene_Name, UCSC_RefGene_Group, sep = ';')
promoter <- distinct(test[test$Regulatory_Feature_Group=="Promoter_Associated" | grepl("TSS200|TSS1500", test$UCSC_RefGene_Group), c("Name", "UCSC_RefGene_Name")])
promoter_cpg <- unique(promoter$Name)

#We then classify the remaining cpgs as enhancers:
test <- test[!test$Name %in% promoter_cpg,]
enhancer <- distinct(test[test$Phantom5_Enhancers!="", c("Name", "UCSC_RefGene_Name")])
enhancer_cpg <- unique(enhancer$Name)

#We then assign the other cpg as either gene body or intergenic:
test <- test[!test$Name %in% enhancer_cpg,]
body <- test[grepl("Body|1stExon|5URT|3UTR|ExonBnd", test$UCSC_RefGene_Group), c("Name", "UCSC_RefGene_Name")]
body_cpg <- unique(body$Name)

test <- test[!test$Name %in% body_cpg,]
intergenic_cpg <- test$Name

identical(sum(length(promoter_cpg), length(enhancer_cpg), length(body_cpg), length(intergenic_cpg)), nrow(annotation))

#Save file with probe id and gene
to_save <- rbind(cbind(promoter, category="promoter"), cbind(enhancer, category="enhancer"), cbind(body, category="gene_body"))
to_save <- distinct(to_save)
saveRDS(to_save, "output/gene_probe_pairs.rds")

fisher_function <- function(variable, direction){
  print(variable)
  if(variable=="promoter"){
    type_df <- promoter_cpg
    other_type <- annotation[!annotation$IlmnID %in% promoter_cpg,"IlmnID"]
  } else if(variable=="enhancer"){
    type_df <- enhancer_cpg
    other_type <- annotation[!annotation$IlmnID %in% enhancer_cpg,"IlmnID"]
  } else if(variable=="gene_body"){
    type_df <- body_cpg
    other_type <- annotation[!annotation$IlmnID %in% body_cpg,"IlmnID"]
  } else if(variable=="intergenic"){
    type_df <- intergenic_cpg
    other_type <- annotation[!annotation$IlmnID %in% intergenic_cpg,"IlmnID"]
  }else if(variable=="island"){
    type_df <- annotation[annotation$Relation_to_UCSC_CpG_Island=="Island","IlmnID"]
    other_type <- annotation[!annotation$Relation_to_UCSC_CpG_Island=="Island","IlmnID"]
    # type_df <- annotation[annotation$Relation_to_UCSC_CpG_Island=="Island" & !annotation$IlmnID %in% promoter_cpg,"IlmnID"]
    # other_type <- annotation[!annotation$Relation_to_UCSC_CpG_Island=="Island" & !annotation$IlmnID %in% promoter_cpg,"IlmnID"]
  } else if(variable=="shelf"){
    type_df <- annotation[annotation$Relation_to_UCSC_CpG_Island %in% c("S_Shelf", "N_Shelf"),"IlmnID"]
    other_type <- annotation[!annotation$Relation_to_UCSC_CpG_Island %in% c("S_Shelf", "N_Shelf"),"IlmnID"]
  } else if(variable=="shore"){
    type_df <- annotation[annotation$Relation_to_UCSC_CpG_Island %in% c("S_Shore", "N_Shore"),"IlmnID"]
    other_type <- annotation[!annotation$Relation_to_UCSC_CpG_Island %in% c("S_Shore", "N_Shore"),"IlmnID"]
  } else if(variable=="open_sea"){
    type_df <- annotation[!annotation$Relation_to_UCSC_CpG_Island %in% c("Island", "S_Shelf", "N_Shelf", "S_Shore", "N_Shore"),"IlmnID"]
    other_type <- annotation[annotation$Relation_to_UCSC_CpG_Island %in% c("Island", "S_Shelf", "N_Shelf", "S_Shore", "N_Shore"),"IlmnID"]
  } 
  
  if(direction=="hypo"){
    type_diff <- length(type_df[type_df %in% rownames(res[res$adj.P.Val<0.05 & res$logFC<0,])])
    type_notdiff <- length(type_df) - type_diff
    other_type_diff <- length(other_type[other_type %in% rownames(res[res$adj.P.Val<0.05 & res$logFC<0,])])
    other_type_notdiff <- length(other_type) - other_type_diff
  } else if(direction=="hyper"){
    type_diff <- length(type_df[type_df %in% rownames(res[res$adj.P.Val<0.05 & res$logFC>0,])])
    type_notdiff <- length(type_df) - type_diff
    other_type_diff <- length(other_type[other_type %in% rownames(res[res$adj.P.Val<0.05 & res$logFC>0,])])
    other_type_notdiff <- length(other_type) - other_type_diff
  }
  
  ### test significance
  m <- matrix(c(type_diff, type_notdiff, other_type_diff, other_type_notdiff), 2,2, byrow = T)
  f <- fisher.test(m)
  return(list("f" = f, "m" = type_diff))
}

families <- c("promoter", "enhancer", "gene_body", "intergenic", "island", "shelf", "shore", "open_sea")
fisher_results <- lapply(families, function(region) fisher_function(region, "hypo"))
names(fisher_results) <- families
saveRDS(fisher_results, paste0('output/final_enrichment_hypo_', tissue,'.rds'))

fisher_results <- lapply(families, function(region) fisher_function(region, "hyper"))
names(fisher_results) <- families
saveRDS(fisher_results, paste0('output/final_enrichment_hyper_', tissue,'.rds'))


#Smoking-Age-DMPs
age <- readRDS(paste0("tissues/", tissue , "/DML_results.rds"))[["AGE"]]

fisher_function_age <- function(variable, direction){
  print(variable)
  if(variable=="promoter"){
    type_df <- promoter_cpg
    other_type <- annotation[!annotation$IlmnID %in% promoter_cpg,"IlmnID"]
  } else if(variable=="enhancer"){
    type_df <- enhancer_cpg
    other_type <- annotation[!annotation$IlmnID %in% enhancer_cpg,"IlmnID"]
  } else if(variable=="gene_body"){
    type_df <- body_cpg
    other_type <- annotation[!annotation$IlmnID %in% body_cpg,"IlmnID"]
  } else if(variable=="intergenic"){
    type_df <- intergenic_cpg
    other_type <- annotation[!annotation$IlmnID %in% intergenic_cpg,"IlmnID"]
  }else if(variable=="island"){
    type_df <- annotation[annotation$Relation_to_UCSC_CpG_Island=="Island","IlmnID"]
    other_type <- annotation[!annotation$Relation_to_UCSC_CpG_Island=="Island","IlmnID"]
  } else if(variable=="shelf"){
    type_df <- annotation[annotation$Relation_to_UCSC_CpG_Island %in% c("S_Shelf", "N_Shelf"),"IlmnID"]
    other_type <- annotation[!annotation$Relation_to_UCSC_CpG_Island %in% c("S_Shelf", "N_Shelf"),"IlmnID"]
  } else if(variable=="shore"){
    type_df <- annotation[annotation$Relation_to_UCSC_CpG_Island %in% c("S_Shore", "N_Shore"),"IlmnID"]
    other_type <- annotation[!annotation$Relation_to_UCSC_CpG_Island %in% c("S_Shore", "N_Shore"),"IlmnID"]
  } else if(variable=="open_sea"){
    type_df <- annotation[!annotation$Relation_to_UCSC_CpG_Island %in% c("Island", "S_Shelf", "N_Shelf", "S_Shore", "N_Shore"),"IlmnID"]
    other_type <- annotation[annotation$Relation_to_UCSC_CpG_Island %in% c("Island", "S_Shelf", "N_Shelf", "S_Shore", "N_Shore"),"IlmnID"]
  } 
  
  if(direction=="hypo"){
    probes_hypo_hypo <- rownames(res[res$adj.P.Val<0.05 & res$logFC<0,])[rownames(res[res$adj.P.Val<0.05 & res$logFC<0,]) %in% rownames(age[age$adj.P.Val<0.05 & age$logFC<0,])]
    type_diff <- length(type_df[type_df %in% probes_hypo_hypo])
    type_notdiff <- length(type_df) - type_diff
    other_type_diff <- length(other_type[other_type %in% probes_hypo_hypo])
    other_type_notdiff <- length(other_type) - other_type_diff
  } else if(direction=="hyper"){
    probes_hyper_hyper <- rownames(res[res$adj.P.Val<0.05 & res$logFC>0,])[rownames(res[res$adj.P.Val<0.05 & res$logFC>0,]) %in% rownames(age[age$adj.P.Val<0.05 & age$logFC>0,])]
    type_diff <- length(type_df[type_df %in% probes_hyper_hyper])
    type_notdiff <- length(type_df) - type_diff
    other_type_diff <- length(other_type[other_type %in% probes_hyper_hyper])
    other_type_notdiff <- length(other_type) - other_type_diff
  }
  
  ### test significance
  m <- matrix(c(type_diff, type_notdiff, other_type_diff, other_type_notdiff), 2,2, byrow = T)
  f <- fisher.test(m)
  return(list("f" = f, "m" = type_diff))
}


fisher_results <- lapply(families, function(region) fisher_function_age(region, "hypo"))
names(fisher_results) <- families
saveRDS(fisher_results, 'output/enrichment_anno_age_hypo.rds')
 
fisher_results <- lapply(families, function(region) fisher_function_age(region, "hyper"))
names(fisher_results) <- families
saveRDS(fisher_results, 'output/enrichment_anno_age_hyper.rds')
