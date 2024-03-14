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
# annotation <- read.csv(paste0(data_path, "GPL21145_MethylationEPIC_15073387_v-1-0.csv"))
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
# 
# #Keep only probes that we consider in the analysis
# results_DML <- readRDS(paste0("tissues/Lung/DML_results.rds"))
# 
# final_annotation <- final_annotation[final_annotation$IlmnID %in% rownames(results_DML$Smoking2),]
# write.csv(final_annotation, paste0(data_path, "GPL21145_MethylationEPIC_15073387_v-1-0_processed_jose.csv"))


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

annotation <- read.csv(paste0(data_path, "GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv"))
Sys.time()

tissue <- "Lung"
# tissue <- "ColonTransverse"
results_DML <- readRDS(paste0("tissues/", tissue , "/DML_results.rds"))

chromhmm <- read.delim('~/Documents/mn4/projects/bsc83/Projects/GTEx_v8/Methylation/Data/EpiMap/Lung.BSS01190_18_CALLS_segments.bed.gz', sep='\t', header=F) #One donor of lung from EpiMap to run ChromHMM #for colon we should use another file
# chromhmm <- read.delim('~/Documents/mn4/projects/bsc83/Projects/GTEx_v8/Methylation/Data/EpiMap/ColonTransverse.BSS01848_18_CALLS_segments.bed.gz', sep='\t', header=F) 

### overlap chromhmm and cpgs tested #### 
ann_bed <- annotation[!is.na(annotation$MAPINFO) & !is.na(annotation$CHR),] %>%
  dplyr::select(chrom=CHR, start=MAPINFO, end=MAPINFO, name=IlmnID) %>% distinct()
ann_bed$chrom <- paste0('chr',ann_bed$chrom)
head(ann_bed)
ann_bed$start <- ann_bed$start-1 #Bed format is 0-based for start and 1-based for end

library(valr)
chrom_df <- chromhmm[,c(1:4)]
colnames(chrom_df) <- c('chrom','start','end','region')
chromhmm_cpgs <- bed_intersect(ann_bed, chrom_df, suffix = c("_ann", "_chromhmm"))

#From 18 states to 14:
chromhmm_cpgs$region_chromhmm[chromhmm_cpgs$region_chromhmm %in% c("TssFlnkD", "TssFlnk", "TssFlnkU")] <- "Flanking TSS"
chromhmm_cpgs$region_chromhmm[chromhmm_cpgs$region_chromhmm %in% c("TssBiv")] <- "Bivalent TSS"
chromhmm_cpgs$region_chromhmm[chromhmm_cpgs$region_chromhmm %in% c("TssA")] <- "Active TSS"

chromhmm_cpgs$region_chromhmm[chromhmm_cpgs$region_chromhmm %in% c("EnhA2", "EnhA1")] <- "Active enhancer"
chromhmm_cpgs$region_chromhmm[chromhmm_cpgs$region_chromhmm %in% c("EnhWk")] <- "Weak enhancer"
chromhmm_cpgs$region_chromhmm[chromhmm_cpgs$region_chromhmm %in% c("EnhBiv")] <- "Bivalent enhancer"
chromhmm_cpgs$region_chromhmm[chromhmm_cpgs$region_chromhmm %in% c("EnhG1", "EnhG2")] <- "Genic enhancer"

chromhmm_cpgs$region_chromhmm[chromhmm_cpgs$region_chromhmm %in% c("ReprPCWk")] <- "Weak repressed polycomb"
chromhmm_cpgs$region_chromhmm[chromhmm_cpgs$region_chromhmm %in% c("ReprPC")] <- "Repressed polycomb"

chromhmm_cpgs$region_chromhmm[chromhmm_cpgs$region_chromhmm %in% c("Quies")] <- "Quiescent"
chromhmm_cpgs$region_chromhmm[chromhmm_cpgs$region_chromhmm %in% c("Het")] <- "Heterochromatin"
chromhmm_cpgs$region_chromhmm[chromhmm_cpgs$region_chromhmm %in% c("ZNF/Rpts")] <- "ZNF genes & repeats"

chromhmm_cpgs$region_chromhmm[chromhmm_cpgs$region_chromhmm %in% c("TxWk")] <- "Weak transcription"
chromhmm_cpgs$region_chromhmm[chromhmm_cpgs$region_chromhmm %in% c("Tx")] <- "Strong transcription"

#Save chromHMM in lung:
write_csv(chromhmm_cpgs, "data/public/lung_chromhmm.csv")

#What chromatin states have our EPIC enhancers?
# annotation_promoter <- annotation$IlmnID[annotation$Regulatory_Feature_Group=="Promoter_Associated"]
# enhancers <- annotation$IlmnID[annotation$Phantom5_Enhancers!=""]
# enhancers <- enhancers[!enhancers%in% annotation_promoter]
# chromhmm_enhacers <- chromhmm_cpgs[chromhmm_cpgs$name_ann %in% enhancers,]
# table(chromhmm_enhacers$region_chromhmm)/table(chromhmm_cpgs$region_chromhmm)
# table(chromhmm_enhacers$region_chromhmm)/length(enhancers)
#35% are active enhancers

#Are smoking- and age-DEGs enriched in CpG island because of bivalent states?
# smoking <- rownames(results_DML$Smoking2[results_DML$Smoking2$adj.P.Val<0.05,])
# age <- rownames(results_DML$AGE[results_DML$AGE$adj.P.Val<0.05,])
# common <- age[age %in% smoking]
# #From common, do a fisher test, are bivalent enriched in cpg? If we remove the bivalent, where are the other CpGs and do we still see enrichment in CpG?
# common_annotation <- chromhmm_cpgs[chromhmm_cpgs$name_ann %in% common,]
# common_annotation_cpg <- annotation[annotation$Name %in% common,]
# bivalent <- common_annotation$name_ann[common_annotation$region_chromhmm %in% c("Bivalent TSS", "Bivalent enhancer", "Repressed polycomb")]
# cpg <- common_annotation_cpg$Name[common_annotation_cpg$Relation_to_UCSC_CpG_Island == "Island"]
# c <- sum(cpg %in% bivalent)
# matrix <- matrix(c(c,
#                  length(bivalent)-c,
#                  length(cpg)-c,
#                  length(common)-length(bivalent)-length(cpg)), nrow=2)
# fisher.test(matrix)
  
#### enrichment

my_fisher <- function(type, direction){
  print(type)
  res <- results_DML[["Smoking2"]]
  # res <- results_DML[["AGE"]]
  chrom_tissue <- as.data.frame(chromhmm_cpgs)
  type_df <- chrom_tissue[chrom_tissue$region_chromhmm == type & chrom_tissue$name_ann %in% rownames(res),]
  other_type <- chrom_tissue[chrom_tissue$region_chromhmm != type & chrom_tissue$name_ann %in% rownames(res),]
  if(direction=="hypo"){
    type_diff <- nrow(type_df[type_df$name_ann %in% rownames(res[res$adj.P.Val<0.05 & res$logFC<0,]),])
    type_notdiff <- nrow(type_df) - type_diff
    other_type_diff <- nrow(other_type[other_type$name_ann %in% rownames(res[res$adj.P.Val<0.05 & res$logFC<0,]),])
    other_type_notdiff <- nrow(other_type) - other_type_diff
  } else if(direction=="hyper"){
    type_diff <- nrow(type_df[type_df$name_ann %in% rownames(res[res$adj.P.Val<0.05 & res$logFC>0,]),])
    type_notdiff <- nrow(type_df) - type_diff
    other_type_diff <- nrow(other_type[other_type$name_ann %in% rownames(res[res$adj.P.Val<0.05 & res$logFC>0,]),])
    other_type_notdiff <- nrow(other_type) - other_type_diff
  }

  ### test significance
  m <- matrix(c(type_diff, type_notdiff, other_type_diff, other_type_notdiff), 2,2, byrow = T)
  print(m)
  f <- fisher.test(m)
  # print(f)
  return(list("f" = f, "m" = type_diff))
}

# Two-tailed Fisher test
families <- as.vector(unique(chromhmm_cpgs$region_chromhmm))
fisher_results <- lapply(families, function(region) my_fisher(region, "hypo"))
names(fisher_results) <- families
saveRDS(fisher_results, paste0("output/enrichment_chromhmm_hypo_", tissue,".rds"))
# saveRDS(fisher_results, paste0("output/enrichment_chromhmm_hypo_", tissue,"_only_age.rds"))

fisher_results <- lapply(families, function(region) my_fisher(region, "hyper"))
names(fisher_results) <- families
saveRDS(fisher_results, paste0("output/enrichment_chromhmm_hyper_", tissue,".rds"))
# saveRDS(fisher_results, paste0("output/enrichment_chromhmm_hyper_", tissue,"_only_age.rds"))




#Enrichments for Smoking + Age DMPs:

my_fisher_age <- function(type, direction){ #same function as before with a few modifications
  print(type)
  res <- results_DML[["Smoking2"]]
  age <- results_DML[["AGE"]]
  chrom_tissue <- as.data.frame(chromhmm_cpgs)
  type_df <- chrom_tissue[chrom_tissue$region_chromhmm == type & chrom_tissue$name_ann %in% rownames(res),] #rownames(res) are all tested positions
  other_type <- chrom_tissue[chrom_tissue$region_chromhmm != type & chrom_tissue$name_ann %in% rownames(res),]
  if(direction=="hypo"){
    probes_hypo_hypo <- rownames(res[res$adj.P.Val<0.05 & res$logFC<0,])[rownames(res[res$adj.P.Val<0.05 & res$logFC<0,]) %in% rownames(age[age$adj.P.Val<0.05 & age$logFC<0,])]
    type_diff <- nrow(type_df[type_df$name_ann %in% probes_hypo_hypo,])
    type_notdiff <- nrow(type_df) - type_diff
    other_type_diff <- nrow(other_type[other_type$name_ann %in% probes_hypo_hypo,])
    other_type_notdiff <- nrow(other_type) - other_type_diff
  } else if(direction=="hyper"){
    probes_hyper_hyper <- rownames(res[res$adj.P.Val<0.05 & res$logFC>0,])[rownames(res[res$adj.P.Val<0.05 & res$logFC>0,]) %in% rownames(age[age$adj.P.Val<0.05 & age$logFC>0,])]
    type_diff <- nrow(type_df[type_df$name_ann %in% probes_hyper_hyper,])
    type_notdiff <- nrow(type_df) - type_diff
    other_type_diff <- nrow(other_type[other_type$name_ann %in% probes_hyper_hyper,])
    other_type_notdiff <- nrow(other_type) - other_type_diff
  }
  
  ### test significance
  m <- matrix(c(type_diff, type_notdiff, other_type_diff, other_type_notdiff), 2,2, byrow = T)
  f <- fisher.test(m)
  # print(f)
  return(list("f" = f, "m" = type_diff))
}
families <- as.vector(unique(chromhmm_cpgs$region_chromhmm))
fisher_results <- lapply(families, function(region) my_fisher_age(region, "hypo"))
names(fisher_results) <- families
saveRDS(fisher_results, 'output/enrichment_chromhmm_hypo_age.rds')

fisher_results <- lapply(families, function(region) my_fisher_age(region, "hyper"))
names(fisher_results) <- families
saveRDS(fisher_results, 'output/enrichment_chromhmm_hyper_age.rds')

