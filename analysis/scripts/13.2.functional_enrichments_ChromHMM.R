#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Functional enrichment differentially methylated positions/loci
# @software version: R=4.2.2

library(dplyr)
library(tidyr)
library(tidyverse)
library(missMethyl)
library(ggplot2)

#Set path 
# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")


#Functional enrichments

#Reading annotation
Sys.time()
data_path <- "~/Documents/mn4/scratch/bsc83/bsc83535/GTEx/v9/Oliva/"
annotation <- read.csv(paste0(data_path, "GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv"))
Sys.time()

#Reading data
results <- readRDS("tissues/Lung/DML_results.rds")
smoking2 <- results$Smoking2
signif <- smoking2[smoking2$adj.P.Val<0.05,]
# results <- readRDS("tissues/ColonTransverse/DML_results.rds")
# smoking2 <- results$Smoking2
# signif <- smoking2[smoking2$adj.P.Val<0.05,]

chromhmm <- read.delim('~/Documents/mn4/projects/bsc83/Projects/GTEx_v8/Methylation/Data/EpiMap/Lung.BSS01190_18_CALLS_segments.bed.gz', sep='\t', header=F) #One donor of lung from EpiMap to run ChromHMM #for colon we should use another file

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



#Functions to plot

plot_go <- function(go, direction){
  go <- go[go$ONTOLOGY=="BP",] #Only BP
  sig <- sum(go$FDR<0.05)
  print(sig)
  if(sig>20){
    max<-20
  } else{max<-sig}
  topgo <- topGSA(go, max)
  ggplot(data = topgo, aes(x=DE/N, y = factor(TERM, levels=rev(TERM)), 
                           color = FDR, size = DE)) + 
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    theme_bw() + ylab("") + xlab("Gene Ratio") + 
    ggtitle(paste0("GO ", direction, " BP"))
}


#Enrichment separated by genomic location:
smoking_hypo <- rownames(signif[signif$logFC<0,])
smoking_hyper <- rownames(signif[signif$logFC>0,])

results <- data.frame("ONTOLOGY"="1", "TERM"="1", "N"="1", "DE"="1", "P.DE"="1", "FDR"="1", "region"="1", "direction"="1")

for(region in unique(chromhmm_cpgs$region_chromhmm)){
  print(region)
  subset <- chromhmm_cpgs$name_ann[chromhmm_cpgs$region_chromhmm==region]
  subset_hypo <- subset[subset %in% smoking_hypo]
  subset_hyper <- subset[subset %in% smoking_hyper]
  go_hypo <- gometh(sig.cpg=subset_hypo, all.cpg = subset, collection="GO", array.type="EPIC") # Background: analysed CpGs per region
  # go_hypo <- gometh(sig.cpg=subset_hypo, all.cpg = annotation$IlmnID, collection="GO", array.type="EPIC") # Background: analysed CpGs in the array
  # go_hypo <- gometh(sig.cpg=subset_hypo, all.cpg = annotation_subset$IlmnID, collection="GO", array.type="EPIC") # Background: DMPs
  go_hypo <- go_hypo[go_hypo$ONTOLOGY=="BP",]
  print(paste(region, "hypo", sum(go_hypo$FDR<0.05)))
  go_hyper <- gometh(sig.cpg=subset_hyper, all.cpg = subset, collection="GO", array.type="EPIC") #Output with genes takes a lot
  # go_hyper <- gometh(sig.cpg=subset_hyper, all.cpg = annotation$IlmnID, collection="GO", array.type="EPIC") #Output with genes takes a lot
  # go_hyper <- gometh(sig.cpg=subset_hyper, all.cpg = annotation_subset$IlmnID, collection="GO", array.type="EPIC") #Output with genes takes a lot
  go_hyper <- go_hyper[go_hyper$ONTOLOGY=="BP",]
  print(paste(region, "hyper", sum(go_hyper$FDR<0.05)))
  
  if(sum(go_hyper$FDR<0.05)!=0){results <- rbind(results, cbind(go_hyper[go_hyper$FDR<0.05,], "region"=region, "direction"="hyper"))}
  if(sum(go_hypo$FDR<0.05)!=0){results <- rbind(results, cbind(go_hypo[go_hypo$FDR<0.05,], "region"=region, "direction"="hypo"))}
}

results <- results[results$FDR<0.05,]
results <- results[order(results$region, results$direction, results$FDR),]
results$Trait <- "Smoking"
results_smoking <- results


#Enrichments for age-smoking-DMPs
res <- readRDS("tissues/Lung/DML_results.rds")
age <- res$AGE[res$AGE$adj.P.Val<0.05,]

smoking_hypo <- rownames(signif[signif$logFC<0,])[rownames(signif[signif$logFC<0,]) %in% rownames(age[age$logFC<0,])]
smoking_hyper <- rownames(signif[signif$logFC>0,])[rownames(signif[signif$logFC>0,]) %in% rownames(age[age$logFC>0,])]

function_to_run_enrichments <- function(smoking_hypo, smoking_hyper, trait){
  results <- data.frame("ONTOLOGY"="1", "TERM"="1", "N"="1", "DE"="1", "P.DE"="1", "FDR"="1", "region"="1", "direction"="1")
  
  for(region in unique(chromhmm_cpgs$region_chromhmm)){
    print(region)
    subset <- chromhmm_cpgs$name_ann[chromhmm_cpgs$region_chromhmm==region]
    subset_hypo <- subset[subset %in% smoking_hypo]
    subset_hyper <- subset[subset %in% smoking_hyper]
    go_hyper <- gometh(sig.cpg=subset_hyper, all.cpg = subset, collection="GO", array.type="EPIC") #Output with genes takes a lot
    go_hyper <- go_hyper[go_hyper$ONTOLOGY=="BP",]
    print(paste(region, "hyper", sum(go_hyper$FDR<0.05)))
    if(sum(go_hyper$FDR<0.05)!=0){results <- rbind(results, cbind(go_hyper[go_hyper$FDR<0.05,], "region"=region, "direction"="hyper"))}
    go_hypo <- try(gometh(sig.cpg=subset_hypo, all.cpg = subset, collection="GO", array.type="EPIC"), silent=T) # Background: analysed CpGs per region
    if(class(go_hypo)!="try-error"){
      go_hypo <- go_hypo[go_hypo$ONTOLOGY=="BP",]
      print(paste(region, "hypo", sum(go_hypo$FDR<0.05)))
      if(sum(go_hypo$FDR<0.05)!=0){results <- rbind(results, cbind(go_hypo[go_hypo$FDR<0.05,], "region"=region, "direction"="hypo"))}
    }
  }
  
  results <- results[results$FDR<0.05,]
  results <- results[order(results$region, results$direction, results$FDR),]
  results$Trait <- "Smoking-age"
  return(results)
}

results_age <- function_to_run_enrichments(smoking_hypo, smoking_hyper, trait="Smoking-age")



#If we remove the CpGs sitting in Polycomb TFBS do we still see development in hypermethylation?

#Removing CpGs sitting in TFBS from Polycomb members:
Sys.time()
tfbs <- read.csv("~/Documents/mn4/scratch/bsc83/bsc83535/GTEx/v9/Oliva/chip_seq_processed.csv") #Lung specific
colnames(tfbs)[8] <- c("TF")
tfbs <- tfbs[,c(5,8)]
tfbs <- tfbs %>% distinct()
Sys.time()
polycomb_genes <- c(#these first are part of PRC1
    "CBX2", "CBX4", "CBX7", "CBX8",
    "BMI1", "PCGF4", "RNF51", #these refer to the same
    "MBLR", "PCGF6", "RING6", "RNF134", #these refer to the same
    "MEL18", "PCGF2", #these refer to the same
    "NSPC", "PCGF1", "RNF68", #these refer to the same
    "RNF1", "RING1A", "RING1", #these refer to the same
    "RNF2", "RING1B", "RING2", #these refer to the same
    "RYBP", "YAF2",
    #these next are part of PRC2
    "EZH1", "EZH2", "SUZ12", "EED", "JARID2",
    "RbAp46", "RBBP7", #these refer to the same
    "RbAp48", "RBBP4", #these refer to the same
    "PCL13",
    #others
    "YY1",
    "PHC1", "PHC2", "KDM2B", "CBX", "REST"
  )
# non_polycomb_cpgs <- tfbs$name_ann[!tfbs$TF %in% polycomb_genes]
polycomb_cpgs <- unique(tfbs$name_ann[tfbs$TF %in% polycomb_genes])
non_polycomb_cpgs <- rownames(res$Smoking2)[!rownames(res$Smoking2) %in% polycomb_cpgs]

smoking_hypo <- rownames(signif[signif$logFC<0,])
smoking_hyper <- rownames(signif[signif$logFC>0,])

function_to_run_enrichments_no_poly <- function(smoking_hypo, smoking_hyper, trait){
  results <- data.frame("ONTOLOGY"="1", "TERM"="1", "N"="1", "DE"="1", "P.DE"="1", "FDR"="1", "region"="1", "direction"="1")
  
  for(region in unique(chromhmm_cpgs$region_chromhmm)){
    print(region)
    subset <- chromhmm_cpgs$name_ann[chromhmm_cpgs$region_chromhmm==region]
    subset <- subset[subset %in% non_polycomb_cpgs]
    subset_hypo <- subset[subset %in% smoking_hypo]
    subset_hyper <- subset[subset %in% smoking_hyper]
    if(length(subset_hypo)>0){
      go_hypo <- try(gometh(sig.cpg=subset_hypo, all.cpg = subset, collection="GO", array.type="EPIC"), silent=T) # Background: analysed CpGs per region
      if(class(go_hypo)!="try-error"){
        go_hypo <- go_hypo[go_hypo$ONTOLOGY=="BP",]
        print(paste(region, "hypo", sum(go_hypo$FDR<0.05)))
        if(sum(go_hypo$FDR<0.05)!=0){results <- rbind(results, cbind(go_hypo[go_hypo$FDR<0.05,], "region"=region, "direction"="hypo"))}
      }
    }
    if(length(subset_hyper)>0){
      go_hyper <- gometh(sig.cpg=subset_hyper, all.cpg = subset, collection="GO", array.type="EPIC") #Output with genes takes a lot
      go_hyper <- go_hyper[go_hyper$ONTOLOGY=="BP",]
      print(paste(region, "hyper", sum(go_hyper$FDR<0.05)))
      if(sum(go_hyper$FDR<0.05)!=0){results <- rbind(results, cbind(go_hyper[go_hyper$FDR<0.05,], "region"=region, "direction"="hyper"))}
      
    }
  }
  
  results <- results[results$FDR<0.05,]
  results <- results[order(results$region, results$direction, results$FDR),]
  results$Trait <- trait
  return(results)
}
results_no_poly <- function_to_run_enrichments_no_poly(smoking_hypo, smoking_hyper, trait="Smoking_no_polycomb")

# res <- readRDS("tissues/Lung/DML_results.rds")
age <- res$AGE[res$AGE$adj.P.Val<0.05,]
smoking_hypo <- rownames(signif[signif$logFC<0,])[rownames(signif[signif$logFC<0,]) %in% rownames(age[age$logFC<0,])]
smoking_hyper <- rownames(signif[signif$logFC>0,])[rownames(signif[signif$logFC>0,]) %in% rownames(age[age$logFC>0,])]
results_age_no_poly <- function_to_run_enrichments_no_poly(smoking_hypo, smoking_hyper, trait="Smoking_age_no_polycomb")


#It is not Generatio and BgRatio! Do it properly
results <- rbind(results_smoking, results_age, results_no_poly, results_age_no_poly)
results$ID <- rownames(results)
results <- results[,c(9, 8, 7, 10, 2, 3, 4, 5, 6)]
results$GeneRatio <- as.numeric(results[,7])/as.numeric(results[,6])
results <- results[,c(1:5,10,8:9)]

colnames(results) <- c("Trait", "Direction", "Region", "ID", "Description", "Ratio", "pvalue", "p.adjust")

library("xlsx")
write.xlsx(results, "output/Supplementary_table_7.xlsx", 
           col.names = TRUE, append = FALSE, row.names = F)


