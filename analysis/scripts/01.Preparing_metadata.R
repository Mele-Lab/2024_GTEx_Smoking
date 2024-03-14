#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to prepare the necessary metadata per tissue
# @software version: R=4.2.2

#Set path
setwd(system("pwd", intern = T)) #If in linux
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")

# Reading tissue information: names, abbreviations and colors
tissue_info <- read.csv("data/public/tissue_abreviation.txt")
tissues <- tissue_info$SMTSD #4 tissues that we do not use
tissues <- tissues[!tissues %in% c("Cells - EBV-transformed lymphocytes", "Cells - Cultured fibroblasts")] #Non-tissues are excluded

#Reading protected metadata
donor_metadata <- read.delim("data/protected/GTEx_Subject_Phenotypes.GRU.txt.gz")
donor_metadata <- donor_metadata[, colnames(donor_metadata) %in% c("SUBJID", "DTHHRDY", "AGE", "SEX", "BMI")] #Variables of interest

sample_metadata <- read.delim("data/protected/GTEx_Sample_Attributes.GRU.txt.gz")
ancestry_metadata <- read.delim("data/protected/inferred_ancestry_838donors.txt")
  
#Reading and filtering gene information
gene_annotation <- read.delim("data/public/gencode.v26.GRCh38.genes.bed")
colnames(gene_annotation) <- c("chr","start","end","strand","feature","ensembl.id","gene.name", "biotype","source") #Renaming variables
#Keeping only PC genes and lincRNAs
gene_annotation <- gene_annotation[gene_annotation$biotype %in% c("protein_coding","lincRNA"),] 
#Excluding PAR genes
PAR_genes <- sapply(grep(".Y", gene_annotation$ensembl.id, value = T), function(gene)
  unlist(strsplit(gene, split = "_"))[[1]])
gene_annotation <- gene_annotation[-unlist(lapply(PAR_genes, function(gene)
  grep(gene, gene_annotation$ensembl.id))),]
write.csv(gene_annotation, "data/public/gene_annotation.csv", row.names = F)


#Reading expression data (these files are not shared by us. They can be downloaded directly from the GTEx portal: https://gtexportal.org/home/datasets)
# data_dir <- "~/Documents/mn4/scratch/bsc83/bsc83535/GTEx/v8/" #My path to these heavy files
print("About to read expression data")
data_dir <- "/gpfs/scratch/bsc83/bsc83535/GTEx/v8/" #My path to these heavy files
tpm <- read.delim(paste0(data_dir, "expression_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"), skip=2)
colnames(tpm) <- gsub("\\.","-", colnames(tpm))
rownames(tpm) <- tpm[,1]
tpm <- tpm[,-c(1:2)]
print("tpm read")

counts <- read.delim(paste0(data_dir, "expression_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"), skip = 2)
colnames(counts) <- gsub("\\.","-", colnames(counts))
rownames(counts) <- counts[,1]
counts <- counts[,-c(1:2)]
print("counts read")

# Subset TPM and count matrix for the genes we are interested in
tpm <- tpm[rownames(tpm) %in% gene_annotation$ensembl.id,]
counts <- counts[rownames(counts) %in% gene_annotation$ensembl.id,]

#Reading the PEER factors. The files are publicly available in the GTEx portal under "GTEx_Analysis_v8_eQTL_covariates.tar.gz"
peer_path <- paste0(data_dir, "cis_QTLs/cis_eQTLs/GTEx_Analysis_v8_eQTL_covariates/") # This is the path where we have located the extracted information

#Reading smoking annotation:
library("readxl")
smoking <- read_excel("data/protected/exSmoker_annotation.xlsx") #814 Manual curation of the protected smoking annotation
smoking <- smoking[!smoking$MHSMKTP %in% c("Cigar", "Pipe", "Other"),] #792  Excluding non-cigarette smokers
smoking <- smoking[!smoking$SmokerStatus=="unknown",] #718
smoking$Smoking <- 0  #never-smokers
smoking$Smoking[smoking$SmokerStatus=="smoker"] <- 2
smoking$Smoking[smoking$SmokerStatus=="ex-smoker"] <- 1
smoking <- smoking[,c("SUBJID", "SmokerStatus", "Smoking")]
smoking <- smoking[!is.na(smoking$SUBJID),]
names(smoking)[1] <- c("Donor")
write.table(smoking, "data/protected/Donor_IDs_with_smoking_status.txt", quote = F, row.names = F) #File that will be used in downstream analysis
smoking <- smoking[,-2]
smoking$Smoking <- as.factor(smoking$Smoking) #We need smoking as a factor for downstream analysis

# Parsing the clinical data into the format we want: From one variable with all the diseases (comma separated) to one variable per disease (levels 0=healthy and 1=disease)
clinical_data <- read.csv("data/public/histology_annotation.csv") # File downloaded directly from the GTEx portal: https://gtexportal.org/home/histologyPage# 
df1 <- na.omit(stack(setNames(strsplit(clinical_data$Pathology.Categories, ","), seq_len(nrow(clinical_data))))[, 2:1])
df1$values <- gsub(" ","",df1$values)
tab <- as.data.frame.matrix(table(df1))
clinical_data <- cbind(clinical_data, tab)

clinical_data[clinical_data$monckeberg==1,"sclerotic"] <- 1  #Monckeberg samples should be annotated as sclerotic, as Monckeberg is a type of sclerosis
#Atherosclerosis, atherosis, calcification and sclerotic (including Monckeberg) will be combined into a single variable that I will call atherosclerosis
clinical_data$atherosclerosis <- as.numeric(clinical_data$atherosclerosis | clinical_data$atherosis | clinical_data$sclerotic | clinical_data$calcification)
clinical_data <- clinical_data[,-which(colnames(clinical_data) %in% c("atherosis", "calcification", "sclerotic", "monckeberg"))]
#Remove non relevant phenotypes (e.g., non diseases)
# clinical_data <- clinical_data[!colnames(clinical_data) %in% c("clean_specimens", "macrophages", "pigment", "no_abnormalities", "post_menopausal", "monckeberg", "spermatogenesis", "congestion")] #spermatogenesis is not correctly annotated, a 1 in the variable sometimes refers to "active spermatogenesis" in the pathology notes and some others to "reduced spermatogenesis", as these categories are automatically extracted from the pathology annotations
clinical_data <- clinical_data[!colnames(clinical_data) %in% c("clean_specimens", "macrophages", "pigment", "no_abnormalities", "post_menopausal", "spermatogenesis", "congestion")] #spermatogenesis is not correctly annotated, a 1 in the variable sometimes refers to "active spermatogenesis" in the pathology notes and some others to "reduced spermatogenesis", as these categories are automatically extracted from the pathology annotations

#Save this histology annotation
# write.csv(clinical_data, "data/public/histological_data.csv", na="NA")

#Function to call later on 
keep_diseases <- function(table){
  if(length(table)==2 & sum(table>20)==2){
    return(TRUE)
  } else{return(FALSE)}
}

create_metadata <- function(tissue){ #Function to get the metadata per tissue
  metadata_subset <- sample_metadata[sample_metadata$SMTSD==tissue,] #Sample metadata for the tissue of interest
  metadata_subset <- metadata_subset[,colnames(metadata_subset) %in% c("SAMPID", "SMTSISCH", "SMRIN", "SMEXNCRT")] #Variables of interest
  metadata_subset$SUBJID <- sapply(metadata_subset$SAMPID, function(i) paste(unlist(strsplit(i,split = "-"))[1:2],collapse="." ) ) #Getting donor ID based on the sample ID
  
  #Adding PEER information in the metadata
  tissue_id <- tissue
  tissue_id <- gsub("\\(|\\)|-", "", tissue_id) #Replace "(", ")" and "-"
  tissue_id <- gsub("[[:space:]]+", "_", tissue_id) #Replace any amount of blank for an "_"
  if(tissue_id=="Brain_Spinal_cord_cervical_c1"){ 
    tissue_id <- "Brain_Spinal_cord_cervical_c-1" #Adding back the "-"
  }
  if(length(grep(tissue_id, list.files(peer_path)))==0){return(NA)}
  print(tissue)
  peer_metadata <- read.delim(list.files(peer_path, full.names = T)[grep(tissue_id, list.files(peer_path))])
  peer_metadata <- peer_metadata[6:7,2:ncol(peer_metadata)] #PEER1 and PEER2
  peer_metadata <- t(peer_metadata)
  colnames(peer_metadata) <- c("PEER1", "PEER2")
  metadata_subset <- merge(metadata_subset, peer_metadata, by.x="SUBJID", by.y="row.names")
  
  #Adding donor metadata:
  metadata_subset$SUBJID <- gsub("\\.", "-", metadata_subset$SUBJID)
  metadata_subset <- merge(metadata_subset, donor_metadata, by="SUBJID")
  
  #Adding genetically inferred ancestry:
  metadata_subset <- merge(metadata_subset, ancestry_metadata, by.x="SUBJID", by.y="ID")
  metadata_subset <- metadata_subset[metadata_subset$inferred_ancestry!="ASN",]
  
  #Renaming and reordering variables
  colnames(metadata_subset) <- c("Donor", "Sample", "RIN", "IschemicTime", "ExonicRate", "PEER1", "PEER2", "Sex", "Age", "BMI", "HardyScale", "Ancestry")
  metadata_subset <- metadata_subset[,c("Donor", "Sample", "HardyScale", "IschemicTime", "RIN", "ExonicRate", "PEER1", "PEER2", "Age", "Ancestry", "Sex", "BMI")]
  metadata_subset$HardyScale <- as.factor(metadata_subset$HardyScale)
  metadata_subset$Ancestry <- as.factor(metadata_subset$Ancestry)
  metadata_subset$Sex <- as.factor(metadata_subset$Sex)
  
  metadata_subset <- na.omit(metadata_subset) #Removing samples with missing data
  
  #Adding smoking annotation
  metadata_subset <- merge(metadata_subset, smoking, by= "Donor")
  
  if(nrow(metadata_subset)<80){ #We keep tissues with at least 80 samples
    print(paste("Skipping", tissue))
    return(NA)
  }
  
  tissue_id <- tissue_info$tissue[tissue_info$SMTSD==tissue] #Getting tissue name without space or _
  
  #Create a directory where all analysis will be performed
  dir.create(paste0("tissues/", tissue_id), showWarnings = FALSE) #If /tissues/ exists
  
  #Subsetting the tpm and counts per tissue
  count_tissue <- counts[,colnames(counts) %in% metadata_subset$Sample]
  tpm_tissue <- tpm[,colnames(tpm) %in% metadata_subset$Sample]
  
  #Subset of metadata for which we have expression data
  metadata_subset <- metadata_subset[metadata_subset$Sample %in% colnames(count_tissue),]
  #Renaming sample IDs
  metadata_subset$Sample <- sapply(metadata_subset$Sample, function(i) paste(unlist(strsplit(i,split = "-"))[1:3],collapse="-" ) ) #This shorter id will be to match the clinical annotation
  colnames(count_tissue) <- sapply(colnames(count_tissue), function(i) paste(unlist(strsplit(i,split = "-"))[1:3],collapse="-" ) ) #This shorter id will be to match the clinical annotation
  colnames(tpm_tissue) <- sapply(colnames(tpm_tissue), function(i) paste(unlist(strsplit(i,split = "-"))[1:3],collapse="-" ) ) #This shorter id will be to match the clinical annotation
  
  #Removing genes that are lowly expressed:
  #TPM>=0.1 in at least 20% of the tissue samples. 
  exprs_genes.tpm <- rownames(tpm)[apply(tpm, 1, function(x) sum(x>=0.1) ) >= 0.2*ncol(tpm)  ]
  #Count >=6 in at least 20% of the tissue samples. 
  exprs_genes.counts <- rownames(counts)[ apply(counts, 1, function(x) sum(x>=6) ) >= 0.2*ncol(counts)  ]
  exprs_genes <- intersect(exprs_genes.tpm,
                           exprs_genes.counts)
  
  # Excluding chrY genes in tissues associated to females
  if(tissue %in% c("Vagina", "Uterus", "Ovary")){
    Y_genes <- gene_annotation[gene_annotation$chr=="chrY",]$ensembl.id 
    exprs_genes <- exprs_genes[!exprs_genes %in% Y_genes]
  }
  
  #saving expression count data
  count_tissue <- count_tissue[exprs_genes,]
  tpm_tissue <- tpm_tissue[exprs_genes,]
  saveRDS(count_tissue, paste0("tissues/", tissue_id, "/counts.rds"))
  saveRDS(tpm_tissue, paste0("tissues/", tissue_id, "/tpm.rds"))
  saveRDS(exprs_genes, paste0("tissues/", tissue_id, "/expressed_genes.rds"))
  
  #Adding clinical traits to metadata
  histology_annotation_tissue <- clinical_data[clinical_data$Tissue == tissue, ] #Annotation for our tissue of itnerest
  histology_annotation_tissue <- histology_annotation_tissue[histology_annotation_tissue$Tissue.Sample.ID %in% metadata_subset$Sample, ] #Annotation for our samples of interest
  histology_annotation_tissue <- rbind(histology_annotation_tissue, clinical_data[clinical_data$Subject.ID %in% metadata_subset$Donor[!metadata_subset$Donor %in% histology_annotation_tissue$Subject.ID] & clinical_data$Tissue==tissue, ])

  count_tables <- mapply(table, histology_annotation_tissue[10:ncol(histology_annotation_tissue)])
  count_tables <- count_tables[sapply(count_tables, function(table) keep_diseases(table))] #Clinical traits with more than 20 healthy and diseased individuals
  if(length(count_tables)>0){ #If there are clinical traits that fit our criteria
    histology_annotation_tissue_to_save <- histology_annotation_tissue[,colnames(histology_annotation_tissue) %in% c("Subject.ID", names(count_tables))]
    #Clinical traits need to be factors for downstream analysis
    for(disease in colnames(histology_annotation_tissue_to_save)){ histology_annotation_tissue_to_save[[disease]] <- as.factor(histology_annotation_tissue_to_save[[disease]]) }
    
    #Merging the previous annotation to the disease annotation
    metadata_subset <- merge(metadata_subset, histology_annotation_tissue_to_save, by.x="Donor", by.y="Subject.ID")
  }
  saveRDS(metadata_subset, paste0("tissues/", tissue_id, "/metadata.rds"))
  return(metadata_subset)
}

metadata <- lapply(tissues, create_metadata)
