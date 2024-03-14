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
annotation <- read.csv(paste0(data_path, "GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv"))
Sys.time()

#Christiansen uses EPIC 800k positions
library(readxl)
Christiansen <- read_excel("data/literature/Christiansen.xlsx", sheet = 2)
names(Christiansen) <- Christiansen[2,]
Christiansen <- Christiansen[3:nrow(Christiansen),]
trues <- lengths(regmatches(Christiansen$`Effect (+ hypermethylation, - hypomethylation)`, gregexpr("-", Christiansen$`Effect (+ hypermethylation, - hypomethylation)`))) > 4
Christiansen_up <- Christiansen$IlmnID[!trues]
Christiansen_down <- Christiansen$IlmnID[trues]

#Joehannes uses 450k
# blood <- read_excel("data/literature/Blood_Sup_2.xlsx", sheet = 3)
# colnames(blood) <- blood[2,]
# blood <- blood[-c(1,2),]
# # blood <- blood[blood$`Probe ID` %in% background,] #Only 94% of the probes in 450K are in EPIC
# blood_up <- blood[blood$Effect>0,]
# blood_down <- blood[blood$Effect<0,]
# blood_up <- blood_up$`Probe ID`
# blood_down <- blood_down$`Probe ID`

# blood_up <- c(blood_up, Christiansen_up)
# blood_down <- c(blood_down, Christiansen_down)


# annotation <- annotation[annotation$Methyl450_Loci==T,]

#The only female in PBMC: BSS01419
chromhmm <- read.delim('~/Documents/mn4/projects/bsc83/Projects/GTEx_v8/Methylation/Data/EpiMap/PBMC.Blood.BSS01419_18_CALLS_segments.bed.gz', sep='\t', header=F)

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
write_csv(chromhmm_cpgs, "data/public/blood_chromhmm_all.csv")
chrom_tissue <- as.data.frame(chromhmm_cpgs)

#### enrichment

my_fisher <- function(type, direction){
  print(type)
  # res <- results_DML[["Smoking2"]]
  # res <- results_DML[["AGE"]]
  type_df <- chrom_tissue[chrom_tissue$region_chromhmm == type & chrom_tissue$name_ann %in% chromhmm_cpgs$name_ann,] #Our background will be only the ones we also test
  other_type <- chrom_tissue[chrom_tissue$region_chromhmm != type & chrom_tissue$name_ann %in% chromhmm_cpgs$name_ann,]
  if(direction=="hypo"){
    type_diff <- nrow(type_df[type_df$name_ann %in% Christiansen_down,])
    type_notdiff <- nrow(type_df) - type_diff
    other_type_diff <- nrow(other_type[other_type$name_ann %in% Christiansen_down,])
    other_type_notdiff <- nrow(other_type) - other_type_diff
  } else if(direction=="hyper"){
    type_diff <- nrow(type_df[type_df$name_ann %in% Christiansen_up,])
    type_notdiff <- nrow(type_df) - type_diff
    other_type_diff <- nrow(other_type[other_type$name_ann %in% Christiansen_up,])
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
saveRDS(fisher_results, paste0("output/enrichment_chromhmm_hypo_blood.rds"))

fisher_results <- lapply(families, function(region) my_fisher(region, "hyper"))
names(fisher_results) <- families
saveRDS(fisher_results, paste0("output/enrichment_chromhmm_hyper_blood.rds"))




# library(missMethyl)
# 
# smoking_hypo <- blood_up
# smoking_hyper <- blood_down
# smoking_hypo <- Christiansen_up
# smoking_hyper <- Christiansen_down
# 
# results <- data.frame("ONTOLOGY"="1", "TERM"="1", "N"="1", "DE"="1", "P.DE"="1", "FDR"="1", "region"="1", "direction"="1")
# 
# for(region in unique(chromhmm_cpgs$region_chromhmm)){
#   print(region)
#   subset <- chromhmm_cpgs$name_ann[chromhmm_cpgs$region_chromhmm==region]
#   subset_hypo <- subset[subset %in% smoking_hypo]
#   subset_hyper <- subset[subset %in% smoking_hyper]
#   go_hypo <- gometh(sig.cpg=subset_hypo, all.cpg = subset, collection="GO", array.type="EPIC") # Background: analysed CpGs per region
#   # go_hypo <- gometh(sig.cpg=subset_hypo, all.cpg = subset, collection="GO", array.type="450K") # Background: analysed CpGs per region
#   go_hypo <- go_hypo[go_hypo$ONTOLOGY=="BP",]
#   print(paste(region, "hypo", sum(go_hypo$FDR<0.05)))
#   go_hyper <- gometh(sig.cpg=subset_hyper, all.cpg = subset, collection="GO", array.type="EPIC")
#   # go_hyper <- gometh(sig.cpg=subset_hyper, all.cpg = subset, collection="GO", array.type="450K")
#   go_hyper <- go_hyper[go_hyper$ONTOLOGY=="BP",]
#   print(paste(region, "hyper", sum(go_hyper$FDR<0.05)))
#   
#   if(sum(go_hyper$FDR<0.05)!=0){results <- rbind(results, cbind(go_hyper[go_hyper$FDR<0.05,], "region"=region, "direction"="hyper"))}
#   if(sum(go_hypo$FDR<0.05)!=0){results <- rbind(results, cbind(go_hypo[go_hypo$FDR<0.05,], "region"=region, "direction"="hypo"))}
# }
# 
# results <- results[results$FDR<0.05,]
# results <- results[order(results$region, results$direction, results$FDR),]
# 
