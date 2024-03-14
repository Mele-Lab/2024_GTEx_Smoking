#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez and Rogério Ribeiro
# @E-mail: jose.ramirez1@bsc.es and rogerio.e.ramos.ribeiro@gmail.com
# @Description: Code to run PCA analysis in Lung
# @software version: R=4.2.2

suppressPackageStartupMessages(library(tidyverse))
library(data.table)
library(psych)
library(ComplexHeatmap)


gtex_data <- fread(file = "public/expression/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz") #Expression table in TPM
metadata <- readRDS(file = "output/metadata.rds")
biotype <- read.delim("data/public/gencode.v26.GRCh38.genes.bed")
colnames(biotype) <- c("chr","start","end","strand","feature","ensembl.id","gene.name", "biotype","source") #Renaming variables

lung_metadata <- metadata$Lung

biotype.pcod_lincRNA <- biotype %>% filter(biotype %in%  c("protein_coding","lincRNA"))


row.names(gtex_data) <- gtex_data$Name

gtex_data.filtered <- gtex_data %>% 
  filter(Name %in% biotype.pcod_lincRNA$ensembl.id) %>% 
  select(-Name, Description) %>% 
  t()

row.names(gtex_data.filtered) <- gsub("-SM-.*", "", row.names(gtex_data.filtered))
gtex_data.filtered <- gtex_data.filtered[lung_metadata$Sample, ]


sample_names <- row.names(gtex_data.filtered)

gtex_data.filtered <- gtex_data.filtered %>% 
  as_tibble() %>%
  mutate_at(.vars = 1:ncol(gtex_data.filtered), as.numeric)


## Run PCA
#Remove rows with unit 0 variance
unit0variancegenes <- unname(apply(gtex_data.filtered, 2, sd) == 0)

pca.lung <- prcomp(log(gtex_data.filtered[,!unit0variancegenes] +1), center = TRUE) 
summ = summary(pca.lung)
pc1.var = summ$importance[2,1]
pc2.var = summ$importance[2,2]


smoking_samples <- metadata$Lung %>% 
  select(Sample, Smoking)

all(smoking_samples$Sample == sample_names) #Order of the samples is the same
labels <- ifelse(smoking_samples$Smoking == 2, "Smoker", ifelse(smoking_samples$Smoking  == 0, "Never Smoker", "Ex Smoker"))


fviz_pca_ind(pca.lung, geom = "point", col.ind = labels, pointsize = 2.5, invisible="quali") + 
  scale_shape_manual(values = c(16,16,16)) + 
  scale_color_manual(values = c("green", "royalblue", "tomato3")) + 
  coord_fixed() + 
  theme(axis.title = element_text(size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.text = element_text(size = 16), 
        legend.position = "top", 
        legend.title = element_blank()) +
  xlab(paste0("PC1 (", round(pc1.var*100,2), "%)")) + ylab(paste0("PC2 (", round(pc2.var*100,2), "%)")) +
  ggtitle("")



## Correlation of the first 10 PC with the covars
covars <- colnames(metadata$Lung)[c(2,4:17)]

pca.lung.pcs <- cbind(Sample = sample_names, pca.lung$x %>% as.data.frame())

PCS_to_correlated <- paste0("PC", seq(1,20, 1))

lung_gtex_pcs <- pca.lung.pcs %>%
  dplyr::select(Sample,all_of(PCS_to_correlated))

pca.lung.pcs <- pca.lung.pcs %>%
  as_tibble() %>%
  gather(key = "PC", value = "value", -Sample) %>%
  filter(PC == "PC1") %>%
  inner_join(metadata$Lung[, covars], by = "Sample")


pca.lung.pcs.covars <- pca.lung.pcs %>%
  dplyr::select(-starts_with("PC")) %>%
  mutate(Sex = factor(Sex, labels = c("1", "2"))) %>%
  mutate(Smoking = factor(Smoking, labels = c("Smoking Never Smoker", "Smoking Ex Smoker", "Smoking Smoker"))) %>%
  mutate(Ancestry = factor(Ancestry, labels = c("AncestryAFR", "AncestryAMR", "AncestryEUR"))) %>%
  mutate_if(is.character, as.factor)

one_hot_encoded <- model.matrix(~ Smoking + Ancestry - 1, data = pca.lung.pcs.covars)

# bind the one-hot encoded columns with the original numeric columns
pca.lung.pcs.covars <- cbind(one_hot_encoded, pca.lung.pcs.covars[,covars])
pca.lung.pcs.covars <- pca.lung.pcs.covars %>% 
  select(-c(Sample,Ancestry ,Smoking )) %>% 
  mutate_all(as.numeric)


lung_gtex_cors <- cor(lung_gtex_pcs[, -c(1)], pca.lung.pcs.covars, use = "pairwise.complete.obs")
lung_gtex_cors.pvalue <- corr.test(lung_gtex_pcs[, -c(1)], pca.lung.pcs.covars, adjust = "BH",  use = "pairwise.complete.obs") #uses the psych library
lung_gtex_cors.pvalue.adjust <- lung_gtex_cors.pvalue$p.adj ##Inspect this df to check if cor is statistical significant

importance_percent <- round(summ$importance[2,] * 100, 2)

row.names(lung_gtex_cors) <- paste0(row.names(lung_gtex_cors), "(", importance_percent[1:20], "%)")

Heatmap(
  matrix = lung_gtex_cors,
  name = "Pearson's r",
  cluster_rows = F,
  cluster_columns = F,
  row_names_gp = gpar(fontsize = 12),
  column_names_gp = gpar(fontsize = 12),
  rect_gp = gpar(col = "white", lwd = 2))



# Save the results
data <- list()
data$PCA <- list()
data$PCA$PCA <- pca.lung
data$PCA$PCA_sum <- summ
data$PCA$PCA_labels <- labels

data$heatmap <- lung_gtex_cors

saveRDS(data, "../figures/data/pca.rds")

