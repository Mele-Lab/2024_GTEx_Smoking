#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to explore reversibility in splicing
# @software version: R=4.2.2

#Set path 
# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")


tissue_info <- read.csv("data/public/tissue_abreviation.txt")
tissues <- tissue_info$tissue

data <- matrix(nrow = length(tissues), ncol=3, dimnames = list(tissues, c("Reversible", "Partially reversible", "Non reversible")))
percentage_closer_to_never <- c()
logFC_never_ex <- c()
logFC_ex_smokers <- c()
tissues_splicing <- c()
for(tissue in tissues){ 
  file <- list.files(paste0("tissues/", tissue), full.names = T)[grepl("fractional_regression",list.files(paste0("tissues/", tissue)))] #Models
  if(length(file)==0){next}
  print(tissue)
  
  metadata <- readRDS(paste0("tissues/", tissue, "/metadata.rds"))
  model <- readRDS(file)
  signif <- model$Smoking2[model$Smoking2$adj.P.Val<0.05,]
  if(nrow(signif)==0){
    data[tissue,] <- rep(NA, 3)
    next
  }
  tissue_id <- tissue_info$Name[tissue_info$tissue==tissue]
  # tissues_splicing <- c(tissues_splicing, paste0(tissue_id, " (", nrow(signif), ")"))
  tissues_splicing <- c(tissues_splicing, tissue_id)
  psi <- readRDS(paste0("tissues/", tissue, "/Alternatively_spliced_events.psi_residuals.rds"))
  psi <- psi[rownames(psi) %in% rownames(signif),]

  # data[tissue, 1] <- nrow(signif)
  
  signif1 <- model$Smoking1[model$Smoking1$adj.P.Val<0.05,] #Never vs ex
  signif2 <- model$`Smoking1-2`[model$`Smoking1-2`$adj.P.Val<0.05,] #Ex vs smokers
  reversed <- signif$Event[signif$Event %in% signif2$Event] #I checked and the few we find are in the same direction
  reversed <- reversed[!reversed %in% signif1$Event]
  
  data[tissue, 1] <- length(reversed) #reversed
  
  non_reversed <- signif$Event[signif$Event %in% signif1$Event]
  non_reversed <- non_reversed[!non_reversed %in% signif2$Event]
  data[tissue, 3] <- length(non_reversed) #non reversed
  
  data[tissue, 2] <- length(signif$Event[!signif$Event %in% reversed & !signif$Event %in% non_reversed])
  
  #What about the PRGs? Are they closer to never smokers?
  psi <- psi[rownames(psi) %in% rownames(signif)[!signif$Event %in% reversed & !signif$Event %in% non_reversed],]
  
  never <- metadata$Sample[metadata$Smoking==0]
  ex <- metadata$Sample[metadata$Smoking==1]
  smoker <- metadata$Sample[metadata$Smoking==2]
  
  #I want two values per gene: logFC never-ex and logFC ex-smokers
  never_m <- rowMeans(as.matrix(psi[,colnames(psi) %in% never]))
  ex_m <- rowMeans(as.matrix(psi[,colnames(psi) %in% ex]))
  smoker_m <- rowMeans(as.matrix(psi[,colnames(psi) %in% smoker]))
  never_ex <- abs(log2(never_m/ex_m))
  ex_smokers <- abs(log2(ex_m/smoker_m))
  
  logFC_never_ex <- c(logFC_never_ex, never_ex)
  logFC_ex_smokers <- c(logFC_ex_smokers, ex_smokers)
  
  booleans <- never_ex < ex_smokers
  percentage_closer_to_never <- c(percentage_closer_to_never, sum(booleans))
}

#Are ex-smokers levels closer to never smokers or smokers?

#get p value of wilcoxon
test <- wilcox.test(logFC_never_ex, logFC_ex_smokers, paired = TRUE, alternative = "less")
# sum(percentage_closer_to_never)/sum(data[,2], na.rm = T)
# table(logFC_never_ex > logFC_ex_smokers)/length(logFC_ex_smokers)
# 
# to_plot <- cbind(logFC_never_ex, "never_ex")
# to_plot <- rbind(to_plot, cbind(logFC_ex_smokers, "ex_smokers"))
# to_plot <- as.data.frame(to_plot)
# to_plot$V2 <- factor(to_plot$V2, levels = c("never_ex", "ex_smokers"))
# to_plot$logFC_never_ex <- as.numeric(to_plot$logFC_never_ex)
# boxplot <- to_plot


#Cummulative plot
logFC_ns_ex <- cbind(logFC=logFC_never_ex, rank=rank(-logFC_never_ex))
logFC_es_s <- cbind(logFC=logFC_ex_smokers, rank=rank(-logFC_ex_smokers))
logFC_ranked <- merge(logFC_ns_ex, logFC_es_s, by = "rank")





#Remove NA from data and sort by DSEs
data <- data[!is.na(rowSums(data)),]
rownames(data) <- tissues_splicing
data <- data[order(data[,2], decreasing = T),]


#Data for boxplots:
gene_annotation <- read.delim("data/public/gencode.v26.GRCh38.genes.bed", header=F)[,c(6,7)]
colnames(gene_annotation) <- c("gene", "symbol")

functiom_to_get_psi <- function(event, tissue){
  #Reading data
  metadata <- readRDS(paste0("tissues/", tissue, "/metadata.rds"))

  #From ensembl id to gene symbol
  gene_id <- strsplit(event, ";")[[1]][1]
  name <- gene_annotation$symbol[gene_annotation$gene==gene_id]
  
  #Residuals
  # psi <- readRDS(paste0(folder, "/", tissue, "/", tissue, ".Alternatively_spliced_events.psi_residuals.rds"))
  # colnames(psi) <- sapply(colnames(psi), function(i) paste(unlist(strsplit(i,split = "-"))[1:3],collapse="-" ) )
  
  #Raw PSI
  psi <- readRDS(paste0("tissues/",tissue,"/psi_splicing_events.rds"))
  psi_subset <- psi[rownames(psi)==event,]
  
  data <- as.data.frame(cbind(unlist(psi_subset), colnames(psi_subset)))
  data$Smoking <- sapply(data$V2, function(sample) metadata$Smoking[metadata$Sample==sample])
  data$V1 <- as.numeric(data$V1)
  
  return(data)
}

#Reversible in heart:
mtdh <-  functiom_to_get_psi("ENSG00000147649.9;SE:chr8:97691188-97699754:97699852-97706626:+", "HeartAtrialAppendage")
# atxn2 <-  functiom_to_get_psi("ENSG00000204842.14;A5:chr12:111470742-111485265:111470742-111485271:-", "HeartAtrialAppendage")

# mturn <- functiom_to_get_psi("ENSG00000180354.15;AF:chr7:30134810:30135298-30146177:30145801:30145959-30146177:+", "HeartAtrialAppendage") #MTURN

#I will choose the one that looks best
rgs12 <- functiom_to_get_psi("ENSG00000159788.18;RI:chr4:3414072:3414241-3414752:3414844:+", "Lung") #RGS12 - RI
# adgrg5_whole_blood <- functiom_to_get_psi("ENSG00000159618.15;RI:chr16:57565034:57565150-57566599:57566751:+", "WholeBlood") #ADGRG5


#Save data here:
figure_S6 <- list(logFC_ranked, data, mtdh)
names(figure_S6) <- c("logFC_ranked", "data", "mtdh")
saveRDS(figure_S6, paste0("../figures/data/figure_S5.rds"))


#Genes we find DS are in Xu et al.? https://www.sciencedirect.com/science/article/pii/S0888754321003864#s0120 
#From ST6:
library(readxl)
genes <- read_excel("data/public/Xu_et_al.xlsx", sheet=1)
genes <- unique(genes[-1,7]) #genes they find have differential isoform usage
genes <- genes$...7

#iterate over tissues and do fisher's
blood <- readRDS("tissues/WholeBlood/fractional_regression_results.rds")
blood <- blood$Smoking2$Ensembl_id[blood$Smoking2$adj.P.Val<0.05]
blood <- sapply(blood, function(gene) gene_annotation$symbol[gene_annotation$gene==gene])

blood %in% genes

blood <- readRDS("tissues/Lung/fractional_regression_results.rds")
blood <- blood$Smoking2$Ensembl_id[blood$Smoking2$adj.P.Val<0.05]
blood <- sapply(blood, function(gene) gene_annotation$symbol[gene_annotation$gene==gene])
table(unique(blood) %in% genes)

blood <- readRDS("tissues/Thyroid/fractional_regression_results.rds")
blood <- blood$Smoking2$Ensembl_id[blood$Smoking2$adj.P.Val<0.05]
blood <- sapply(blood, function(gene) gene_annotation$symbol[gene_annotation$gene==gene])
table(unique(blood) %in% genes)

blood <- readRDS("tissues/HeartAtrialAppendage/fractional_regression_results.rds")
blood <- blood$Smoking2$Ensembl_id[blood$Smoking2$adj.P.Val<0.05]
blood <- sapply(blood, function(gene) gene_annotation$symbol[gene_annotation$gene==gene])
table(unique(blood) %in% genes)

#I could try with transcripts with differential expression
