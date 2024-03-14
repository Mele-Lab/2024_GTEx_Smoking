#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Methylation results analysis
# @software version: R=4.2.2

#Set path 
# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")

tissue_info <- read.csv("data/public/tissues_sorted.csv")
source("../figures/safe_colourblind_pallete.R")

# Number of DML per tissue

tissues <- c("BreastMammaryTissue", "ColonTransverse", "Lung", "MuscleSkeletal", "Ovary", "Prostate", "Testis", "WholeBlood")
metadata <- list()
table <- data.frame("tissue"=1, "Hypermethylated"=1, "Hypomethylated"=1)

for(tissue in tissues){
  print(tissue)
  # metadata_s <- readRDS(paste0("tissues/", tissue, "/metadata.rds"))
  metadata_s <- readRDS(paste0("tissues/", tissue, "/methylation_metadata.rds"))
  metadata[[tissue]] <- metadata_s
  results <- readRDS(paste0("tissues/", tissue, "/DML_results.rds"))
  hypo <- sum(results$Smoking2$adj.P.Val<0.05 & results$Smoking2$logFC<0)
  hyper <- sum(results$Smoking2$adj.P.Val<0.05 & results$Smoking2$logFC>0)
  table <- rbind(table, c(tissue, hyper, hypo))
}
table <- table[-1,]
table <- as.data.frame(table)
table$Hypermethylated <- as.numeric(table$Hypermethylated)
table$Hypomethylated <- as.numeric(table$Hypomethylated)
rownames(table) <- table$tissue
table$tissue <- sapply(table$tissue, function(tissue) tissue_info$TISSUENAMEABREV[tissue_info$tissue==tissue])

table <- table[order(table$Hypomethylated, decreasing = T),]
to_plot <- table


#Lung results:
results <- readRDS("tissues/Lung/DML_results.rds")

smoking2 <- results$Smoking2
signif <- smoking2[smoking2$adj.P.Val<0.05,]
100*(table(signif$logFC>0)[1]/nrow(signif)) #43%
100*(table(signif$logFC>0)[2]/nrow(signif)) #56%
hypo_hyper <- signif[,c(1,2)]
hypo_hyper[,2] <- "hypomethylation"
hypo_hyper[hypo_hyper$logFC>0,2] <- "hypermethylation"
hypo_hyper$logFC <- abs(hypo_hyper$logFC)
colnames(hypo_hyper) <- c("logFC", "direction")


#Save data here:
figure_5 <- list(metadata, to_plot, hypo_hyper)
names(figure_5) <- c("metadata", "to_plot", "hypo_hyper")
saveRDS(figure_5, paste0("../figures/data/figure_5.rds"))


smoking2 <- rownames(smoking2[smoking2$adj.P.Val<0.05,])
smoking1 <- results$Smoking1
smoking1 <- rownames(smoking1[smoking1$adj.P.Val<0.05,])
smoking2_1 <- results$`Smoking1-Smoking2`
smoking2_1 <- rownames(smoking2_1[smoking2_1$adj.P.Val<0.05,])

sum(smoking1 %in% smoking2)
sum(smoking2_1 %in% smoking2)

#Explore direction
non_rev <- smoking1[smoking1 %in% smoking2]
non_rev <- signif[rownames(signif) %in% non_rev,]
table(non_rev$logFC>0)
rev <- smoking2_1[smoking2_1 %in% smoking2]
rev <- signif[rownames(signif) %in% rev,]
table(rev$logFC>0)


#Lung and Colon overlap:
colon <- readRDS("tissues/ColonTransverse/DML_results.rds")

smoking2_colon <- colon$Smoking2
signif_colon <- smoking2_colon[smoking2_colon$adj.P.Val<0.05,]

hyper_colon <- rownames(signif_colon[signif_colon$logFC>0,])
hypo_colon <- rownames(signif_colon[signif_colon$logFC<0,])
hyper_lung <- rownames(signif[signif$logFC>0,])
hypo_lung <- rownames(signif[signif$logFC<0,])

common <- c(hyper_colon[hyper_colon %in% hyper_lung], hypo_colon[hypo_colon %in% hypo_lung])

m <- matrix(c(length(common),
            length(rownames(signif_colon)[!rownames(signif_colon) %in% common]), 
            length(rownames(signif)[!rownames(signif) %in% common]),
            length(rownames(colon$AGE)[!rownames(colon$AGE) %in% signif & !rownames(colon$AGE) %in% signif_colon])),
            nrow=2)
fisher.test(m)


data_path <- "~/Documents/mn4/scratch/bsc83/bsc83535/GTEx/v9/Oliva/"
annotation <- read.csv(paste0(data_path, "GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv"))
annotation <- annotation[annotation$IlmnID %in% common, c("IlmnID", "UCSC_RefGene_Name")]




# 
# #AHRR:
# 
# 
# ho <- subset_annotation[subset_annotation$Name %in% rownames(smoking_hypo),] #11 AHRR
# hr <- subset_annotation[subset_annotation$Name %in% rownames(smoking_hyper),] #59 AHRR
# 
# #What if high logFC?
# smoking_hypo_h <- smoking_hypo[smoking_hypo$logFC<(-0.8),] #I tried 0.5, 0.8 and 1
# smoking_hyper_h <- smoking_hyper[smoking_hyper$logFC>0.8,]
# 
# ho <- subset_annotation[subset_annotation$Name %in% rownames(smoking_hypo_h),] #11 AHRR
# hr <- subset_annotation[subset_annotation$Name %in% rownames(smoking_hyper_h),] #59 AHRR
# length(grep("AHRR", ho$UCSC_RefGene_Name))
# length(grep("AHRR", hr$UCSC_RefGene_Name))
# 
# #What if very signif p_val:
# smoking_hypo_h <- smoking_hypo[smoking_hypo$adj.P.Val<(0.0005),] #I tried 0.5, 0.8 and 1
# smoking_hyper_h <- smoking_hyper[smoking_hyper$adj.P.Val<(0.0005),]
# 
# 
# non_reversible <- DML_annotation_1[grep("AHRR", DML_annotation_1$UCSC_RefGene_Name),"Name"]
# non_reversible %in% hr$Name
# 
# 
# 
# 
# #Colon Analysis
# results_c <- readRDS(paste0(project_path, "ColonTransverse/DML_results.rds"))
# 
# smoking2_c <- results_c$Smoking2
# signif_c <- smoking2_c[smoking2_c$adj.P.Val<0.05,]
# table(signif_c$logFC>0)
# smoking2_c <- rownames(smoking2_c[smoking2_c$adj.P.Val<0.05,])
# smoking1_c <- results_c$Smoking1
# smoking1_c <- rownames(smoking1_c[smoking1_c$adj.P.Val<0.05,])
# smoking2_1_c <- results_c$`Smoking1-Smoking2`
# smoking2_1_c <- rownames(smoking2_1_c[smoking2_1_c$adj.P.Val<0.05,])
# 
# sum(smoking1_c %in% smoking2_c)
# sum(smoking2_1_c %in% smoking2_c)
# sum(smoking2_1_c %in% smoking1_c)
# 
# table(smoking2_c %in% smoking2)
# common <- smoking2_c[smoking2_c %in% smoking2]
# hyper_colon <- signif_c[rownames(signif_c) %in% common,"logFC"]>0
# hyper_lung <- signif[rownames(signif) %in% common,"logFC"]>0
# table(hyper_colon & hyper_lung)
# fisher.test(matrix(c(39,48, 87460, 780379), nrow=2)) #39 common between colon and lung, 48 colon specific, 87460 lung specific, 780379 non DML (all probes minus the other numbers)
# 
# #Are the 39 in the same direction for both? No, one is not, 38 are
# smoking2_c_t <- results_c$res_M$Smoking2
# smoking2_c <- rownames(smoking2_c[smoking2_c$adj.P.Val<0.05,])
# 
# common <- smoking2_c[smoking2_c %in% smoking2]
# hyper_colon_r <- sapply(common, function(x) smoking2_c_t$logFC[rownames(smoking2_c_t)==x])
# hyper_lung_r <- sapply(common, function(x) signif$logFC[rownames(signif)==x])
# identic <- cbind(hyper_colon_r>0, hyper_lung_r>0) #Only one difference
# fisher.test(matrix(c(38,49, 87461, 780379), nrow=2)) #39 common between colon and lung, 48 colon specific, 87460 lung specific, 780379 non DML (all probes minus the other numbers)
# common <- common[-35]
# 
# #What genes are associated to these 38 positions?
# 
# library(missMethyl)
# test <- getMappedEntrezIDs(common, array.type = "EPIC") 
# length(test$sig.eg) #AHRR is here
# AHRR <- subset_annotation$Name[grep("AHRR", subset_annotation$UCSC_RefGene_Name)] #AHRR associated probes
# interesting <- common[common %in% AHRR] #7 positions!
# smoking2_c_t[rownames(smoking2_c_t) %in% interesting,]
# #I loaded betas to compare the logFC output by modeling Ms and the delta beta to check directionality and:
# #The delta beta is in the same direction that the M logFC, so the interpretation I was doing is correct
# 
# #What about cg03636183 from F2RL3 and cg05575921 for AHRR? Known https://clinicalepigeneticsjournal.biomedcentral.com/articles/10.1186/s13148-020-00908-3
# smoking2[rownames(smoking2)=="cg05575921",] #AHRR
# smoking2[rownames(smoking2)=="cg03636183",] #F2RL3
# smoking2[rownames(smoking2)=="cg23576855",]
# smoking2_1[rownames(smoking2_1)=="cg05575921",]
# 
# #AHRR in adipose AND blood
# smoking2[rownames(smoking2)=="cg11554391",]
# smoking2[rownames(smoking2)=="cg25648203",]
# #In adipose only?
# smoking2[rownames(smoking2)=="cg04135110",]
# #Other genes such as CYP1A1:
# smoking2[rownames(smoking2)=="cg23680900",]
# smoking2[rownames(smoking2)=="cg26516004",] #NO
# smoking2[rownames(smoking2)=="cg10009577",]
# smoking2[rownames(smoking2)=="cg00353139",]

