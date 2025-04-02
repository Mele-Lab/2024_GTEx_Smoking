#!/usr/bin/env Rscript
# @Author: Rogério Eduardo Ramos Ribeiro
# @E-mail: rogerio.e.ramos.ribeiro
# @Description:Generate Principal Factors using MOFA. Divide methylation probes into several types of data.
# @software version: R=4.2.2

setwd("D:/PhD/smoking_analysis.v4")


# Load libraries ----
library(tidyverse)
library(data.table)
library(MOFA2)
library(psych)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(tidyverse)
library(MOFAdata)
library(PCGSE)
library(ggpubr)
library(DOSE)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DESeq2)


normalize_count_data_vst <- function(count_matrix){
  col_data <- data.frame(
    row.names = colnames(count_matrix)
  )
  
  dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = col_data, design = ~ 1)
  dds <- DESeq(dds)
  vsd <- vst(dds, blind = TRUE)
  
  return (assay(vsd))
}


# Function to remove everything after the second "-"
remove_after_second_dash <- function(x) {
  parts <- unlist(strsplit(x, "-"))
  if (length(parts) > 2) {
    return(paste(parts[1:2], collapse = "-"))
  } else {
    return(x)
  }
}


## Enrichment analysis functions
enrichement_analysis_MOFA <- function(model, ontology_data, factor, statistical_test, verbose, view = "gexp"){
  enrichement.positive <- run_enrichment(model,
                                    view = view, 
                                    factors = factor,
                                    feature.sets = ontology_data,
                                    sign = "positive",
                                    statistical.test = statistical_test,
                                    min.size = 5,
                                    alpha = 0.05,
                                    verbose = verbose
  )
  
  enrichement.negative <- run_enrichment(model,
                                         view = view, 
                                         factors = factor,
                                         feature.sets = ontology_data,
                                         sign = "negative",
                                         statistical.test = statistical_test,
                                         min.size = 5,
                                         alpha = 0.05,
                                         verbose = verbose
  )
  
  return(list("positive"=enrichement.positive, "negative"= enrichement.negative))
}

enrichement_analysis_MOFA_methylation <- function(model, probe_list, view_name, ontology_data, factor, statistical_test, verbose){
  
  ## Filter the cpg in the data
  cpg <- probe_list %>% filter(Name %in% row.names(model@expectations$W[[view_name]]))
  ## Merge with ensembl
  cpg <- merge(cpg, biotype, by.x = "UCSC_RefGene_Name", by.y = "Name")
  cpg$ensembl <- gsub("\\.\\d+", "", cpg$ensembl)
  # Make a cpg vector
  cpg_vector <- cpg$Name
  names(cpg_vector) <-  cpg$ensembl
  
  # Replace the colnames in ontology_data
  ontology_data_cpg <- ontology_data[,colnames(ontology_data) %in% cpg$ensembl]
  colnames(ontology_data_cpg) <- unname(cpg_vector[colnames(ontology_data_cpg)])
  
  enrichement.positive <- run_enrichment(model,
                                         view = view_name, 
                                         factors = factor,
                                         feature.sets = ontology_data_cpg,
                                         sign = "positive",
                                         statistical.test = statistical_test,
                                         min.size = 5,
                                         alpha = 0.05,
                                         verbose = verbose
  )
  
  enrichement.negative <- run_enrichment(model,
                                         view = view_name, 
                                         factors = factor,
                                         feature.sets = ontology_data_cpg,
                                         sign = "negative",
                                         statistical.test = statistical_test,
                                         min.size = 5,
                                         alpha = 0.05,
                                         verbose = verbose
  )
  return(list("positive"=enrichement.positive, "negative"= enrichement.negative))
}



### Set directories ###
analysis <- "MOFA"
base_image_folder <- analysis

if (!dir.exists(base_image_folder)){
  dir.create(base_image_folder)
}

#Set variables----
sample_size <- list()
corrs <- list()
res_enrichement <- list()
res_corr <- list()
sex_tissues <- c("testis", "prostate", "ovary")
biotype <- readRDS("metadata/biotype.rds")
data_folder <- "D:/PhD/age_molecular_characterization/data/"

annot <- read.csv(paste0(data_folder, "metadata/methylation_epic_v1.0b5.csv"))
load_resids_from_memory <- TRUE


# Classify probes ----
#We first classify as promoters the CpG annotated as "Promoter_Associated" and TSS. If a CpG is associated to two genes, 
#One as promoter and another as gene body, we keep the CpG only as associated to the gene as promoter and exclude the other gene (very few cases) 
test <- annot %>% separate_rows(UCSC_RefGene_Name, UCSC_RefGene_Group, sep = ';')
promoter <- distinct(test[test$Regulatory_Feature_Group=="Promoter_Associated" | grepl("TSS200|TSS1500", test$UCSC_RefGene_Group), c("Name", "UCSC_RefGene_Name")])
promoter_cpg <- unique(promoter$Name)

#We then classify the remaining cpgs as enhancers:
test <- test[!test$Name %in% promoter_cpg,]
enhancer <- distinct(test[test$Phantom5_Enhancers!="", c("Name", "UCSC_RefGene_Name")])
enhancer_cpg <- unique(enhancer$Name)

#We then assign the other cpg as either gene body or intergenic:
test <- test[!test$Name %in% enhancer_cpg,]
body <- test[grepl("Body|1stExon|ExonBnd|5'UTR|3'UTR", test$UCSC_RefGene_Group), c("Name", "UCSC_RefGene_Name")]
body_cpg <- unique(body$Name)

test <- test[!test$Name %in% body_cpg,]
intergenic <- test[,c("Name", "UCSC_RefGene_Name")]
intergenic_cpg <- test$Name

list_of_cpg <- list("Promoter" = promoter_cpg, 
                    "Enhancer" = enhancer_cpg,
                    "Body" = body_cpg, 
                    "Intergenic" = intergenic_cpg)


n_factors_per_tissue <- list( "lung" = 12)

### Apply MOFA ####

tissue <- "lung"

#for (tissue in names(n_factors_per_tissue)){

cat("Running MOFA analysis for", tissue, "\n")
tissue2 <- tissue
if (tissue == "breast_mammarytissue"){
  tissue2 <- "breast"
}
if (tissue == "colon_transverse"){
  tissue2 <- "colon"
}
if (tissue == "muscle_skeletal"){
  tissue2 <- "muscle"
}

# Load metadata ----
gexp_metadata <- readRDS(file.path(data_folder, "metadata", "processed_metadata", "gexp_metadata.rds"))[[tissue]]
meth_metadata <- readRDS(paste0("metdata/", tissue2, "_methylation_metadata.rds"))


gexp <- fread(paste0(data_folder, "/processed_gexp/", tissue, ".counts.csv")) %>% dplyr::select(-Description) %>% column_to_rownames("Name")
meth <- fread(paste0(data_folder, "/methylation/methylation_", tissue, ".csv")) %>% column_to_rownames("probe")

colnames(gexp)[1:ncol(gexp)] <- unname(sapply(colnames(gexp)[1:ncol(gexp)], remove_after_second_dash))
colnames(meth)[1:ncol(meth)] <- unname(sapply(colnames(meth)[1:ncol(meth)], remove_after_second_dash))

common <- intersect(gexp_metadata$subject_id, meth_metadata$SUBJID)
gexp_metadata <- gexp_metadata %>% filter(subject_id %in% common) %>% 
  filter(SmokerStatus == "Smoker" | SmokerStatus == "Non Smoker")

print(paste("Sample size for ", tissue, " ", nrow(gexp_metadata)))

gexp <- gexp %>% dplyr::select(gexp_metadata$subject_id)
meth <- meth %>% dplyr::select(gexp_metadata$subject_id)


#1 - Process dataset 
## Normalize gene expression data by size factors and apply variance stabilizing transformation
normalized_counts <- normalize_count_data_vst(gexp)
M <- convert_to_M_values(meth)


#2 - Keep Highly variable features

## Divide methylation into enhancers, promoters, gene-body and intergenic 
#list_of_cpg

promoters_cps <- M %>%
  filter(row.names(.) %in% list_of_cpg$Promoter)

enhancers_cpg <- M %>%
  filter(row.names(.) %in% list_of_cpg$Enhancer)

body_cpg <- M %>%
  filter(row.names(.) %in% list_of_cpg$Body)

intergenic_cpg <- M %>%
  filter(row.names(.) %in% list_of_cpg$Intergenic)

## Keep Highly variable features
top_n_variable_features <- function(data, n) {
  feature_variances <- apply(data, 1, var, na.rm = TRUE)
  top_features_indices <- order(feature_variances, decreasing = TRUE)[1:n]
  top_features <- data[top_features_indices, , drop = FALSE]
  
  return(top_features)
}

gexp.top <- top_n_variable_features(normalized_counts, n = 5000)

meth.promoters.top <- top_n_variable_features(promoters_cps, n = 5000)
meth.enhancer.top <- top_n_variable_features(enhancers_cpg, n = 5000)
meth.body.top <- top_n_variable_features(body_cpg, n = 5000)
meth.intergenic.top <- top_n_variable_features(intergenic_cpg, n = 5000)


#3 Prepare MOFA object 
data <- list("gexp" = as.matrix(gexp.top),
             "promoters" = as.matrix(meth.promoters.top),
             "enhancer" = as.matrix(meth.enhancer.top),
             "body" = as.matrix(meth.body.top),
             "intergenic" = as.matrix(meth.intergenic.top)
)

MOFAobject <- create_mofa(data)

## Define data options 
data_opts <- get_default_data_options(MOFAobject)
data_opts$scale_views <- TRUE
## Define models options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- n_factors_per_tissue[[tissue]]
## Define train options
train_opts <- get_default_training_options(MOFAobject) 
train_opts$drop_factor_threshold <- 0.01
train_opts$convergence_mode <- "slow"

MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)


#4 Train the model ----
outfile <- paste0(base_image_folder, "/", tissue, "_model.hdf5")
MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = TRUE)
#MOFAobject.trained <- load_model(outfile)

#5 Add metadata ----

metadata_meth_tissue <- meth_metadata %>% 
  mutate(meth.PEER1 = PEER1, meth.PEER2 = PEER2) %>%
  dplyr::select(SUBJID, meth.PEER1, meth.PEER2)

metadata <- gexp_metadata %>% 
  mutate(gexp.PEER1 = PEER1, gexp.PEER2 = PEER2) %>%
  dplyr::select(-c(sample_id, starts_with("PEER"))) %>% 
  mutate(sample = subject_id) %>% 
  merge(metadata_meth_tissue, by.y= "SUBJID", by.x= "sample") %>% 
  dplyr::select(sample, SEX, AGE, BMI, TRISCHD, DTHHRDY, SmokerStatus, meth.PEER1, meth.PEER2, gexp.PEER1, gexp.PEER2)

metadata$SmokerStatus <- ifelse(metadata$SmokerStatus == "Smoker", 1, 0)

cat("Samples in the correct order:", all(metadata$sample_id == colnames(gexp.top)), "\n")

samples_metadata(MOFAobject.trained) <- metadata


##6 Variance decomposition and other plots ----

#pdf(paste0(base_image_folder, "/", tissue, "_variance_per_modality_per_view.pdf"), h = 6)
#print(plot_variance_explained(MOFAobject.trained, x="view", y="factor"))
#dev.off()


#pdf(paste0(base_image_folder, "/", tissue, "_variance_per_modality.pdf"), h = 4)
#print(plot_variance_explained(MOFAobject.trained, x="group", y="factor", plot_total = T)[[2]])
#dev.off()

#pdf(paste0(base_image_folder, "/", tissue, "_correlation_among_factors.pdf"), h = 4)
#print(plot_factor_cor(MOFAobject.trained))
#dev.off()


##7 Heatmap ----

factors <- get_factors(MOFAobject.trained, factors = "all")$group1

metadata_to_correlate <- metadata %>% 
  dplyr::select(-c(sample, DTHHRDY, TRISCHD)) %>% 
  mutate(SEX = as.numeric(SEX))


# Function to remove constant columns
remove_constant_columns <- function(data) {
  # Calculate the standard deviation for each column
  sds <- apply(data, 2, sd, na.rm = TRUE)
  sds <- na.omit(sds)
  # Keep only columns with non-zero standard deviation
  data <- data[, sds != 0, drop = FALSE]
  
  data
}

# Remove constant columns from factors and metadata_to_correlate
factors <- remove_constant_columns(factors)
metadata_to_correlate <- remove_constant_columns(metadata_to_correlate)

if (tissue %in% c("ovary", "testis", "prostate")){
  metadata_to_correlate <- metadata_to_correlate %>% dplyr::select(-SEX)
}

corr <- psych::corr.test(factors, metadata_to_correlate, method = "spearman", adjust = "fdr")

corr_cor <- corr$r
corr_fdr <- corr$p.adj

res_corr[[tissue]] <- list("cor_spearman" = corr_cor, 
                           "padjust" = corr_fdr)



#8 Enrichement analysis (GSEA).
#statistical_test <- "cor.adj.parametric" # One of "parametric", "permutation" or cor.adj.parametric
statistical_test <- "permutation"
list_of_cpg_dataframe <- list("Promoter" = promoter, 
                    "Enhancer" = enhancer,
                    "Body" = body, 
                    "Intergenic" = intergenic)



row.names(MOFAobject.trained@expectations$W$gexp) <- gsub("\\.\\d+", "", row.names(MOFAobject.trained@expectations$W$gexp))
row.names(MOFAobject.trained@data$gexp$group1) <- gsub("\\.\\d+", "", row.names(MOFAobject.trained@data$gexp$group1))
MOFAobject.trained@features_metadata$feature <- gsub("\\.\\d+", "", MOFAobject.trained@features_metadata$feature)

names(MOFAobject.trained@intercepts$gexp$group1) <- gsub("\\.\\d+", "", names(MOFAobject.trained@intercepts$gexp$group1) )

cols <- c("SYMBOL", "GO")

#Build association tables between GO terms and genes
# For gene expression
gene_expression_to_go <- select(org.Hs.eg.db, keys=row.names(MOFAobject.trained@data$gexp$group1), columns=cols, keytype="ENSEMBL")
gene_expression_to_go <- na.omit(gene_expression_to_go) %>% filter(ONTOLOGY == "BP")

go_matrix <- matrix(0, nrow = length(unique(gene_expression_to_go$ENSEMBL)), ncol = length(unique(gene_expression_to_go$GO)),
                    dimnames = list(unique(gene_expression_to_go$ENSEMBL), unique(gene_expression_to_go$GO)))

for (i in 1:nrow(gene_expression_to_go)) {
  gene <- gene_expression_to_go$ENSEMBL[i]
  go_term <- gene_expression_to_go$GO[i]
  go_matrix[gene, go_term] <- 1
}

go_matrix <- t(go_matrix)

# For Methylation probes 

build_go_table_for_methy <- function(model, probe_list, view_name){
  ## Filter the cpg in the data
  cpg <- probe_list %>% filter(Name %in% row.names(model@expectations$W[[view_name]]))
  ## Merge with ensembl -> Get gene nam
  cpg <- merge(cpg, biotype, by.x = "UCSC_RefGene_Name", by.y = "Name")
  cpg$ensembl <- gsub("\\.\\d+", "", cpg$ensembl)
  
  # For each gene associated with a probe, get the annotation
  methylation_to_go <- select(org.Hs.eg.db, keys=cpg$ensembl, columns=cols, keytype="ENSEMBL")
  methylation_to_go <- na.omit(methylation_to_go) %>% filter(ONTOLOGY == "BP")
  methylation_to_go <- merge(methylation_to_go, cpg, by.x = "ENSEMBL", by.y = "ensembl") %>% 
    mutate(gene_name_combination = paste0(ENSEMBL, Name))
  
  # Build probe -> GO term association
  go_matrix <- matrix(0, nrow = length(unique(methylation_to_go$Name)), ncol = length(unique(methylation_to_go$GO)),
                      dimnames = list(unique(methylation_to_go$Name), unique(methylation_to_go$GO)))
  
  for (i in 1:nrow(methylation_to_go)) {
    gene <- methylation_to_go$Name[i]
    go_term <- methylation_to_go$GO[i]
    go_matrix[gene, go_term] <- 1
  }
  go_matrix <- t(go_matrix)
  return(go_matrix)
}

promoter_to_go <- build_go_table_for_methy(MOFAobject.trained, list_of_cpg_dataframe$Promoter, "promoters")
enhancer_to_go <- build_go_table_for_methy(MOFAobject.trained, list_of_cpg_dataframe$Enhancer, "enhancer")
body_to_go <- build_go_table_for_methy(MOFAobject.trained, list_of_cpg_dataframe$Body, "body")


n_factors <- ncol(get_weights(MOFAobject.trained)[[names(get_weights(MOFAobject.trained)[1])]])
res_enrichement[[tissue]][["Gene Expression"]] <- enrichement_analysis_MOFA(MOFAobject.trained, go_matrix, 1:n_factors, statistical_test, verbose = FALSE)
saveRDS(res_enrichement, paste0(base_image_folder, "/enrichement_MOFA.rds"))
res_enrichement[[tissue]][["Promoters"]] <- enrichement_analysis_MOFA(MOFAobject.trained, promoter_to_go, 1:n_factors, statistical_test, verbose = FALSE, view = "promoters")
saveRDS(res_enrichement, paste0(base_image_folder, "/enrichement_MOFA.rds"))
res_enrichement[[tissue]][["Enhancers"]] <- enrichement_analysis_MOFA(MOFAobject.trained, enhancer_to_go, 1:n_factors, statistical_test, verbose = FALSE, view = "enhancer")
saveRDS(res_enrichement, paste0(base_image_folder, "/enrichement_MOFA.rds"))
res_enrichement[[tissue]][["Gene Body"]] <- enrichement_analysis_MOFA(MOFAobject.trained, body_to_go, 1:n_factors, statistical_test, verbose = FALSE, view = "body")

saveRDS(res_enrichement, paste0(base_image_folder, "../figure/data/enrichement_MOFA.rds"))
saveRDS(res_corr, paste0("../figure/data/corr_mofa_factors.rds"))