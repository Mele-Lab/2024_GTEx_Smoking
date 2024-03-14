#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to correlate gene expression and DNA methylation
# @software version: R=4.2.2

library(valr)
library(dplyr)
library(ggplot2)
library(ggpubr)

#Set path 
# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")

cor <- readRDS("tissues/Lung/Correlations_DMP_DEG_new.rds")

table(cor$p.adj<0.05)
signif <- cor[cor$p.adj<0.05,]
# colon <- signif$probe
head(sort(table(signif$gene), decreasing = T))

length(unique(signif$probe)) #1356 DMPs correlate with DEGs
length(unique(signif$probe))/length(unique(cor$probe))*100 #12%

length(unique(signif$gene)) #523 DEGs correlate with DMPs
length(unique(signif$gene))/length(unique(cor$gene))*100 # 27%

#What is the percentage of anticorrelations?
sum(signif$cor<0)/nrow(signif)

#Compare barplots to see anticorrelations
dmp_pos <- length(unique(signif$probe[signif$cor>0]))/length(unique(cor$probe))*100 
dmp_neg <- length(unique(signif$probe[signif$cor<0]))/length(unique(cor$probe))*100 

deg_pos <- length(unique(signif$gene[signif$cor>0]))/length(unique(cor$gene))*100  #Some genes have both positive and negative correlations with different DMPs
deg_neg <- length(unique(signif$gene[signif$cor<0]))/length(unique(cor$gene))*100

to_plot <- data.frame(value=c(dmp_pos, dmp_neg, deg_pos, deg_neg), 
                      Correlation=c("Positive", "Negative", "Positive", "Negative"),
                      variable=c("DMP", "DMP", "DEG", "DEG"))
saveRDS(to_plot, "../figures/data/correlation_percentages.rds")



#Perform enrichments on promoter/enhancer/gene_body. Are most correlations happening in promoters?
pairs <- readRDS("output/gene_probe_pairs.rds") #This file contains the cpg type: promoter/enhancer/gene_body
res <- readRDS("tissues/Lung/DML_results.rds")[["Smoking2"]]
pairs <- merge(cor, pairs, by.x=c("probe", "gene"), by.y=c("Name", "UCSC_RefGene_Name"))

#Here the question is: from the DMP-DEG pairs, how many significant correlations are from hypomethylated promoters for instance
fisher_function <- function(variable, direction){
  print(variable)
  if(variable=="promoter"){
    type_df <- unique(pairs$probe[pairs$category=="promoter"]) #Maybe one DMP is associated to two DEGs
    other_type <- unique(pairs$probe[pairs$category!="promoter"])
  } else if(variable=="enhancer"){
    type_df <- unique(pairs$probe[pairs$category=="enhancer"])
    other_type <- unique(pairs$probe[pairs$category!="enhancer"])
  } else if(variable=="gene_body"){
    type_df <- unique(pairs$probe[pairs$category=="gene_body"])
    other_type <- unique(pairs$probe[pairs$category!="gene_body"])
  } 
  
  if(direction=="hypo"){
    type_diff <- length(type_df[type_df %in% signif$probe & type_df %in% rownames(res)[res$logFC<0]]) #Is the probe significantly correlated and hypomethylated?
    type_notdiff <- length(type_df) - type_diff
    other_type_diff <- length(other_type[other_type %in% signif$probe & other_type %in% rownames(res)[res$logFC<0]])
    other_type_notdiff <- length(other_type) - other_type_diff
  } else if(direction=="hyper"){
    type_diff <- length(type_df[type_df %in% signif$probe & type_df %in% rownames(res)[res$logFC>0]])
    type_notdiff <- length(type_df) - type_diff
    other_type_diff <- length(other_type[other_type %in% signif$probe & other_type %in% rownames(res)[res$logFC>0]])
    other_type_notdiff <- length(other_type) - other_type_diff
  }
  
  ### test significance
  m <- matrix(c(type_diff, type_notdiff, other_type_diff, other_type_notdiff), 2,2, byrow = T)
  f <- fisher.test(m)
  return(list("f" = f, "m" = type_diff))
}

families <- c("promoter", "enhancer", "gene_body")
hypo <- lapply(families, function(region) fisher_function(region, "hypo"))
names(hypo) <- families

hyper <- lapply(families, function(region) fisher_function(region, "hyper"))
names(hyper) <- families

#save these hypo and hyper
saveRDS(hypo, 'output/enrichment_correlated_hypo_lung.rds')
saveRDS(hyper, 'output/enrichment_correlated_hyper_lung.rds')


#We only see enrichments in hypomethylated enhancers, these are positive or negative correlations?
enh <- pairs[pairs$category=="enhancer",]
enh <- enh[enh$p.adj<0.05,]
enh_hypo <- enh[enh$probe %in% rownames(res)[res$logFC<0 & res$adj.P.Val<0.05],]
table(enh_hypo$cor>0)/nrow(enh_hypo)*100



#Plot rho values:
library(ggplot2)
ggplot(pairs, aes(category, abs(cor))) + geom_violin() + #changes in promoters are not stronger in gene expression, if anything enhancers
  geom_boxplot(outlier.shape = NA,
               notch = T,
               width = 0.25) + theme_bw() + xlab("") + ylab("abs(rho)")

ggplot(pairs[pairs$p.adj<0.05,], aes(category, abs(cor))) + geom_violin() + 
  geom_boxplot(outlier.shape = NA,
               notch = T,
               width = 0.25) + theme_bw() + xlab("") + ylab("abs(rho)")







ecpg <- signif$probe
noncpg <- cor$probe[!cor$probe %in% ecpg]

results_DML <- readRDS("tissues/Lung/DML_results.rds")
smoking <- results_DML$Smoking2
smoking <- smoking[rownames(smoking) %in% cor$probe,] #Probes associated to a gene
smoking$ecpg <- "Non-correlated"
smoking$ecpg[rownames(smoking) %in% ecpg] <- "Correlated"

smoking$probe <- rownames(smoking)
smoking <- merge(smoking, cor, by="probe")


get_box_stats <- function(y, upper_limit) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "N =", length(y), "\n"
    )
  ))
}

smoking$correlation <- "positive"
smoking$correlation[smoking$cor<0] <- "negative"
smoking$dummy <- paste0(smoking$ecpg, smoking$correlation)
smoking$dummy <- factor(smoking$dummy, levels=c("Correlatedpositive", "Non-correlatedpositive", "Correlatednegative", "Non-correlatednegative"))
p_value_1_pos <- wilcox.test(abs(smoking$logFC[smoking$dummy=="Correlatedpositive"]), abs(smoking$logFC[smoking$dummy=="Non-correlatedpositive"]))
p_value_2_pos <- wilcox.test(abs(smoking$logFC[smoking$dummy=="Correlatednegative"]), abs(smoking$logFC[smoking$dummy=="Non-correlatednegative"]))
saveRDS(list("smoking"=smoking, "p_value_1_pos"=p_value_1_pos, "p_value_2_pos"=p_value_2_pos), "../figures/data/logFC_pos_neg.rds")



#The more DMPs the more likely to correlate with expression?
table(pairs$gene)

library(dplyr)
result <- pairs %>%
  group_by(gene) %>%
  # group_by(gene, category) %>%
  summarize(mean_abs_cor = mean(abs(cor)),
            num_instances = n())


lm_fit <- lm(mean_abs_cor ~ num_instances, data = result)
rho <- cor(result$mean_abs_cor, result$num_instances)
p_value <- summary(lm_fit)$coefficients[2, 4] #same as cor.test(result$mean_abs_cor, result$num_instances)

ggplot(result, aes(num_instances, mean_abs_cor)) + geom_point(alpha=0.5) +
  geom_smooth(method = "lm", se = FALSE, col = "red") +
  annotate("text", y = max(result$mean_abs_cor), x = min(result$num_instances),
           label = paste("Correlation (rho):", round(rho, 3), ", P-value:", round(p_value, 5)),
           hjust = -0.25, vjust = 0, col = "blue", size = 3) + theme_bw()


name <-"gene_body"
name <-"enhancer"
name <-"promoter"
test <- result[result$category==name,]
lm_fit <- lm(mean_abs_cor ~ num_instances, data = test)
rho <- cor(test$mean_abs_cor, test$num_instances)
p_value <- summary(lm_fit)$coefficients[2, 4] #same as cor.test(result$mean_abs_cor, result$num_instances)

ggplot(data=test, aes(num_instances, mean_abs_cor)) + geom_point() +
  geom_smooth(method = "lm", se = FALSE, col = "red") + ggtitle(name) +
  annotate("text", y = max(test$mean_abs_cor), x = min(test$num_instances),
           label = paste("Correlation (rho):", round(rho, 3), ", P-value:", round(p_value, 5)),
           hjust = -0.25, vjust = 0, col = "blue", size = 3)


#I did more enrichments and nothing really poped up


#Scatterplot of examples:
Sys.time()
beta <- readRDS("tissues/Lung/methylation_residuals.rds") #Reading beta residuals
Sys.time() 

#Reading expression residuals in the lung:
expression <- readRDS("tissues/Lung/expression_residuals_demographic_traits_no_smoking.rds")

#From ensembl id to gene symbol
gene_annotation <- read.delim("data/public/gencode.v26.GRCh38.genes.bed", header=F)[,c(6,7)]
colnames(gene_annotation) <- c("gene", "symbol")

rownames(expression) <- sapply(rownames(expression), function(gene) gene_annotation$symbol[gene_annotation$gene==gene])
#From sample id to donor id
colnames(expression) <- sapply(colnames(expression), function(id) paste0(strsplit(id, "-")[[1]][-3], collapse="-"))
expression <- expression[,colnames(expression) %in% colnames(beta)]
beta <- beta[,colnames(beta) %in% colnames(expression)]


#other examples: AHRR - cg25648203 or BCL11B - cg16452866 -> highest anticorrelation and p.adjusted
#Plot AHRR - cg07943658 -> highest correlation and p.adjusted
scatter_map4k1 <- beta[rownames(beta)=="cg07943658",] #the name map4k1 comes from a previous code
names(scatter_map4k1) <- colnames(beta)
scatter_map4k1 <- rbind(scatter_map4k1, expression[rownames(expression)=="AHRR", colnames(expression) %in% colnames(beta)])
scatter_map4k1 <- as.data.frame(t(scatter_map4k1))
colnames(scatter_map4k1) <- c("Methylation", "Expression")
scatter_map4k1$Methylation <- as.numeric(scatter_map4k1$Methylation)
scatter_map4k1$Expression <- as.numeric(scatter_map4k1$Expression)

correlation <- round(cor$cor[cor$probe=="cg07943658"], 2)
p_value <- signif(cor$p.adj[cor$probe=="cg07943658"], 2)
saveRDS(list("scatter_map4k1"=scatter_map4k1, "correlation"=correlation, "p_value"=p_value), "../figures/data/AHRR_correlation_example.rds")




#Compute functional enrichment of correlated genes:
library(clusterProfiler)
library(org.Hs.eg.db)

smoking <- merge(smoking, pairs, by=c("probe", "gene"))

results_go <- data.frame("ID"=1, "Description"=1, "GeneRatio"=1, "BgRatio"=1, "pvalue"=1, "p.adjust"=1, "region"=1, "direction"=1)

bg <- unique(smoking$gene)
for(region in families){
  print(region)
  #Enrichment of hypermethylated correlations:
  gl <- unique(smoking$gene[smoking$p.adj.x<0.05 & smoking$logFC>0 & smoking$category==region])
  ora.go <- enrichGO(gene     = gl, universe =    bg,
                     keyType = "SYMBOL", OrgDb        = org.Hs.eg.db,
                     ont          = "BP", 
                     pAdjustMethod = "BH", readable = F)
  if(sum(ora.go@result$p.adjust<0.05)>0){
    print(paste("There are enrichments:", region, "hypermethylation"))
    pdf(paste0("Plots/TFBS_", region, "_", "hypermethylation", ".pdf"))
    print(dotplot(ora.go, showCategory = 15, title=paste0(region, "_", "hypermethylation")))
    dev.off()
    results_go <- rbind(results_go, cbind(ora.go[ora.go$p.adjust<0.05,1:6], region, "direction"="hypermethylation"))
  }
  #Enrichment of negative correlations:
  # gl <- unique(pairs$gene[pairs$p.adj<0.05 & pairs$cor<0 & pairs$category==region])
  gl <- unique(smoking$gene[smoking$p.adj.x<0.05 & smoking$logFC<0 & smoking$category==region])
  ora.go <- enrichGO(gene     = gl, universe =    bg,
                     keyType = "SYMBOL", OrgDb        = org.Hs.eg.db,
                     ont          = "BP", 
                     pAdjustMethod = "BH", readable = F)
  if(sum(ora.go@result$p.adjust<0.05)>0){
    print(paste("There are enrichments:", region, "hypomethylation"))
    pdf(paste0("Plots/TFBS_", region, "_", "hypomethylation", ".pdf"))
    print(dotplot(ora.go, showCategory = 15, title=paste0(region, "_", "hypomethylation")))
    dev.off()
    results_go <- rbind(results_go, cbind(ora.go[ora.go$p.adjust<0.05,1:6], region, "direction"="hypomethylation"))
  }  
  
}
results_go <- results_go[-1,]
results_go <- results_go[,c(7,8,1:6)]
colnames(results_go) <- c("Region", "Direction", "ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust")

library("xlsx")
write.xlsx(results_go, "output/Supplementary_table_9_new.xlsx",
           col.names = TRUE, row.names = F, append = FALSE)
