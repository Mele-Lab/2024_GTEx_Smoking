#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to prepare the data for figure 2
# @software version: R=4.2.2

#Set path 
setwd(system("pwd", intern = T)) #If in linux
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
# setwd("..")

tissue_info <- read.csv("data/public/tissues_sorted.csv")
tissues <- tissue_info$tissue


n_samples <- c()
table <- data.frame("tissue"=1, "DSE"=1, "DSG"=1, "DEG"=1, "DS&EG"=1, "signif"=1, "odds ratio"=1, "p-value"=1)
events <- data.frame("event"="1", "tissue"="1", "gene"="1") #To overlap our DSEs with previous literature
for(i in 1:length(tissues)){ #I am creating a table with the info from all diseases
  tissue <- tissues[i]
  print(tissue)
  metadata <- readRDS(paste0("tissues/", tissue, "/metadata.rds"))

  file <- list.files(paste0("tissues/", tissue, "/"))[grepl("fractional_regression_results.rds",list.files(paste0("tissues/", tissue, "/")))] #We do this in case there is no file for splicing
  
  if(length(file)==0){next}
  n_samples <- c(n_samples, nrow(metadata)) #Computing n
  to_append <- c("Healthy"=sum(metadata[["Smoking"]]==0), "Diseased"=sum(metadata[["Smoking"]]==1))
  dsa_res <- readRDS(paste0("tissues/", tissue, "/", file))
  dsa_res_genes <- dsa_res[["Smoking2"]][dsa_res[["Smoking2"]]$adj.P.Val<0.05,]
  dsa_events <- nrow(dsa_res_genes)
  to_add <- cbind(rownames(dsa_res_genes), rep(tissue, dsa_events), dsa_res_genes$Ensembl_id)
  colnames(to_add) <- c("event", "tissue", "gene")
  events <- rbind(events, to_add)
  dsa_genes <- unique(dsa_res_genes$Ensembl_id)
  #Genes DS and DE:
  dea_res <- readRDS(paste0("tissues/", tissue, "/voom_limma_results.rds"))
  dea_genes <- dea_res[["Smoking2"]][dea_res[["Smoking2"]]$adj.P.Val<0.05,]
  dea_genes <- rownames(dea_genes)
  common_genes <- intersect(dea_genes, dsa_genes)
  #significant overlap?
  potential_spliced_genes <- unique(dsa_res$Age$Ensembl_id)   #Background for fisher's test
  x12 <- dsa_genes[!dsa_genes %in% common_genes]
  x21 <- dea_genes[!dea_genes %in% common_genes] #x21 is hte number of DEG from our background
  x21 <- x21[x21 %in% potential_spliced_genes] #The DEG from the subset of potential spliced genes, in this way we will see if from the potetntailly spliced genes we have more or less overlapping with changes in expression
  x22 <- potential_spliced_genes[!potential_spliced_genes %in% x12]
  x22 <- x22[!x22 %in% x21]
  x22 <- x22[!x22 %in% common_genes]
  
  m <- matrix(ncol = 2, nrow=2, c(length(common_genes), length(x12), length(x21), length(x22)), byrow = T)
  test <- fisher.test(m)
  p_value <- test$p.value
  o_r <- test$estimate
  if(p_value<0.05){
    signif <- 1
  } else{
    signif <- 0 
  }
  table <- rbind(table, c(tissue, dsa_events, length(dsa_genes), length(dea_genes), length(common_genes), signif, o_r, p_value))
}
table <- table[-1,]
events <- events[-1,]

#Barplot
table <- as.data.frame(table)
table$DSE <- as.numeric(table$DSE)
rownames(table) <- table$tissue
table$tissue <- sapply(table$tissue, function(tissue) tissue_info$TISSUENAMEABREV[tissue_info$tissue==tissue])

splicing_table <- table[table$DSE>0,]
splicing_table$tissue <- factor(splicing_table$tissue, levels=splicing_table$tissue[order(splicing_table$DSE, decreasing = T)])

colours <- tissue_info$color
names(colours) <- tissue_info$tissue_abbrv
colours <- colours[names(colours) %in% splicing_table$tissue]
# sort colours to paint in the desidered order
colours <- colours[order(match(names(colours), levels(splicing_table$tissue)))]



#Do expression and splicing overlap? Only in the lung
table[rownames(table)=="Lung",]

# boolean <- disease$Diseased<10

#Piechart for type of event, is smoking different than the other demographic traits?
splicing.events <-  c("SE","MX","AF","AL", "A5","A3","RI")

barplot <- data.frame("tissue"=1, "SE"=1, "MX"=1, "AF"=1, "AL"=1,  "A5"=1, "A3"=1, "RI"=1)
p_values <- data.frame("tissue"=1, "SE"=1, "MX"=1, "AF"=1, "AL"=1,  "A5"=1, "A3"=1, "RI"=1)
for(i in 1:length(tissues)){ #I am creating a table with the info from all diseases
  tissue <- tissues[i]
  print(tissue)
  file <- list.files(paste0("tissues/", tissue, "/"))[grepl("fractional_regression_results.rds",list.files(paste0("tissues/", tissue, "/")))] #We do this in case there is no file for splicing
  if(length(file)==0){next}
  dsa_res <- readRDS(paste0("tissues/", tissue, "/", file))
  dsa_res_types <- dsa_res[["Smoking2"]][dsa_res[["Smoking2"]]$adj.P.Val<0.05,]$Type
  dsa_res_types <- factor(dsa_res_types, levels = splicing.events)
  if(length(dsa_res_types)==0){
    # p_values <- rbind(p_values, c(tissue, 1, 1, 1, 1, 1, 1, 1))
    barplot <- rbind(barplot, c(tissue, summary(as.factor(dsa_res_types))))
    next
  }
  #Does smoking affect different event types?
  age_types <- dsa_res[["Age"]][dsa_res[["Age"]]$adj.P.Val<0.05,]$Type
  EUR_events <- dsa_res[["AncestryEUR"]][dsa_res[["AncestryEUR"]]$adj.P.Val<0.05,]
  AMR_events <- dsa_res[["AncestryAMR"]][dsa_res[["AncestryAMR"]]$adj.P.Val<0.05,]
  ancestry_events <- rbind(EUR_events, AMR_events)
  ancestry_types <- ancestry_events$Type[!duplicated(ancestry_events$Event)]
  BMI_types <- dsa_res[["BMI"]][dsa_res[["BMI"]]$adj.P.Val<0.05,]$Type
  sex_types <- dsa_res[["Sex"]][dsa_res[["Sex"]]$adj.P.Val<0.05,]$Type
  demographic_types <- c(age_types, ancestry_types, BMI_types, sex_types)
  demographic_types <- factor(demographic_types, levels = splicing.events)
  # table(demographic_types)/length(demographic_types)
  # table(dsa_res_types)/length(dsa_res_types)
  p_val <- c()
  for(j in 1:length(splicing.events)){
    demo_p <- table(demographic_types)[j]/length(demographic_types)
    test <- binom.test(table(dsa_res_types)[j], length(dsa_res_types), demo_p)
    p_val <- c(p_val, test$p.value)
  }
  p_values <- rbind(p_values, c(tissue, p_val))
  barplot <- rbind(barplot, c(tissue, summary(as.factor(dsa_res_types))))
}
barplot <- barplot[-1,]
p_values <- p_values[-1,]
rownames(barplot) <- barplot$tissue

#Adjusting p-value for multiple correction:
p_values <- as.data.frame(p_values[,2:8])
p_values <- sapply(p_values, as.numeric)
p_adjusted <- p.adjust(p_values[18,], method="BH")

library(reshape2)
barplot_plot <- barplot
barplot_plot <- melt(barplot_plot, id.vars = c("tissue"))
barplot_plot$tissue <- factor(barplot_plot$tissue, levels=rev(tissues))
barplot_plot$value <- as.numeric(barplot_plot$value)


#pie chart of the colsums
barplot_test <- barplot[,2:8]
for(i in 1:7){
  barplot_test[,i] <- as.numeric(barplot_test[,i])
}
data <- colSums(barplot_test)
data <- as.data.frame(cbind(data, names(data)))
data$data <- as.numeric(data$data)
library(dplyr)

pie_data <- data %>% 
  mutate(csum = rev(cumsum(rev(data))), 
         pos = data/2 + lead(csum, 1),
         pos = if_else(is.na(pos), data/2, pos))




#PSI values:
disease <- "Smoking"

gene_annotation <- read.delim("data/public/gencode.v26.GRCh38.genes.bed", header=F)[,c(6,7)]
colnames(gene_annotation) <- c("gene", "symbol")

functiom_to_get_psi <- function(event, tissue){
  #Reading data
  metadata <- readRDS(paste0("tissues/", tissue, "/metadata.rds"))
  metadata$Smoking <- as.factor(metadata$Smoking)
  
  #From ensembl id to gene symbol
  gene_id <- strsplit(event, ";")[[1]][1]
  name <- gene_annotation$symbol[gene_annotation$gene==gene_id]
  
  #Residuals
  # psi <- readRDS(paste0(folder, "/", tissue, "/", tissue, ".Alternatively_spliced_events.psi_residuals.rds"))
  # colnames(psi) <- sapply(colnames(psi), function(i) paste(unlist(strsplit(i,split = "-"))[1:3],collapse="-" ) )
  
  #Raw PSI
  psi <- readRDS(paste0("tissues/", tissue, "/psi_splicing_events.rds"))
  psi_subset <- psi[rownames(psi)==event,]
  
  data <- as.data.frame(cbind(unlist(psi_subset), colnames(psi_subset)))
  data$Smoking <- sapply(data$V2, function(sample) metadata$Smoking[metadata$Sample==sample])
  data$V1 <- as.numeric(data$V1)
  
  return(data)
}

#Are there some tissue shared events?
shared <- data.frame("event"=1, "n"=1)
shared$n <- as.numeric(shared$n)
for(i in 1:length(tissues)){
  tissue <- tissues[i]
  print(tissue)
  file <- list.files(paste0("tissues/", tissue, "/"))[grepl("fractional_regression_results.rds",list.files(paste0("tissues/", tissue, "/")))] #We do this in case there is no file for splicing
  if(length(file)==0){next}
  dsa_res <- readRDS(paste0("tissues/", tissue, "/", file))
  dsa_res_events <- rownames(dsa_res[["Smoking2"]][dsa_res[["Smoking2"]]$adj.P.Val<0.05,])
  for(event in dsa_res_events){
    if(event %in% shared[,1]){
      shared[shared[,1]==event, 2] <- shared[shared[,1]==event, 2] + 1 
    } else{
      shared[nrow(shared)+1,] <- list(event, 1) #List so that we conserve data types
    }
  }
}
#There is only one shared event in lung and whole blood: ADGRG5
#If we check other splicing events with high logFC, we find RGS12 and CLDN7 in the lung, candidates to plot and report

data_rgs12 <- functiom_to_get_psi("ENSG00000159788.18;RI:chr4:3414072:3414241-3414752:3414844:+", "Lung") #RGS12 - RI
data_cldn7 <- functiom_to_get_psi("ENSG00000181885.18;SE:chr17:7260536-7260642:7260726-7260821:-", "Lung") #CLDN7
data_adgrg5_lung <- functiom_to_get_psi("ENSG00000159618.15;RI:chr16:57565034:57565150-57566599:57566751:+", "Lung") #ADGRG5
data_adgrg5_lung_whole_blood <- functiom_to_get_psi("ENSG00000159618.15;RI:chr16:57565034:57565150-57566599:57566751:+", "WholeBlood") #ADGRG5


metadata <- readRDS("output/metadata.rds")

#Save data here:
figure_2 <- list(splicing_table, colours, pie_data, data_rgs12, data_cldn7, data_adgrg5_lung, data_adgrg5_lung_whole_blood, metadata)
names(figure_2) <- c("splicing_table", "colours", "pie_data", "data_rgs12", "data_cldn7", "data_adgrg5_lung", "data_adgrg5_lung_whole_blood", "metadata")
saveRDS(figure_2, paste0("../figures/data/figure_2_new_splicing.rds"))



#Do we replicate our findings in Xu et al.?
#We check the isoforms that either include or exclude the given event and we check if any of those 
# coordinates <- readRDS("SUPPA/gencode.v26.splicing_events_coordinates.rds")
#No need, just check if our DSGs are genes with DTU in their results
library("readxl")
xu_et_al <- read_excel("data/public/Xu_et_al.xlsx", sheet = 1)
colnames(xu_et_al) <- xu_et_al[1,]
xu_et_al <- xu_et_al[-1,]
xu_et_al <- as.data.frame(xu_et_al)
xu_et_al$FDR <- as.numeric(xu_et_al$FDR)
# xu_et_al <- xu_et_al[xu_et_al$FDR<0.1,]
xu_et_al <- xu_et_al[xu_et_al$FDR<0.05,]

parker_et_al <- c("ERAP1", "LRRN3", "MAN1A2", "UTRN", "AREL1", "SASH1", "GALNT7", "EPS15")
parker_et_al <- sapply(parker_et_al, function(gene) gene_annotation$gene[gene_annotation$symbol==gene])

#Check if the DSEs in table events are found in previous literature
events$gene %in% parker_et_al #None in Parker

events$gene_id <- gsub("\\..*","",events$gene)
events$replicated <- "no"
events$replicated[events$gene_id %in% xu_et_al$gene_id] <- "yes"

replication <- events[events$replicated=="yes",c("tissue", "gene_id")]
replication <- replication[!duplicated(replication),]
table(replication$tissue)
