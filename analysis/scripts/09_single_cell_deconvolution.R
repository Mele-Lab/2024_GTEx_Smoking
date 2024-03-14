#!/usr/bin/env Rscript
# @Author: Rogério Eduardo Ramos Ribeiro and Jose Miguel Ramirez
# @E-mail: zroger499@gmail.com and jose.ramirez1@bsc.es
# @Description: Code to run BayesPrism on the lung single cell analysis
# @software version: R=4.2.2

# Load libraries
library(Seurat)
library(BayesPrism)
suppressPackageStartupMessages(library(tidyverse))

#Set work dir 
setwd(system("pwd", intern = T)) #If in linux
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio


## Create out folder
out.folder <- "output/scRNA_seq/"
if (!dir.exists(out.folder)){
  dir.create(out.folder, recursive = T)
}

# Load single cell data
processed_data <- readRDS("output/scRNA_seq_processed.rds")
featureList <- readRDS(file = "output/featureList_GTEx.rds")

# Load Bulk data 
expression_data <- read.delim("data/expression/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz", skip = 2)
metadata <- readRDS(file = "data/metadata/metadata.rds")
degs <- readRDS(file = "data/DEGS/differential_expression_results.rds")


#Subset bulk expression dataset for expressed genes
expression_data <- expression_data %>% 
  filter(Name %in% row.names(degs$Lung$`SmokingSMOKER-NEVER`))


#Subset tissue samples to filter
tissue_samples <- metadata[["Lung"]]

#Subset dataset
colnames(expression_data) <- gsub("\\.SM.*", "", colnames(expression_data))
colnames(expression_data) <- gsub("\\.", "-", colnames(expression_data))


exp.lung <- expression_data %>% 
  select(Name, all_of(tissue_samples$Sample))

cat("The expression dataset has ", nrow(exp.lung)," genes\n")

# Assign celltype/subtype labels  

## Remove rare cell types (min 50 cell per celltype)
min_cells <- 50
filtered_celltypes <- table(Idents(processed_data))[table(Idents(processed_data)) > min_cells]

processed_data <- processed_data[,unname(Idents(processed_data)) %in% names(filtered_celltypes)] # Select only celltypes that pass this threshold

# Get count assay and filter it for pcoding and lincRNA in common with GTEx
featureList.filtered <- featureList %>% 
  filter(biotype == "protein_coding" | biotype == "lincRNA")

scLung.data <- t(GetAssayData(object = processed_data, slot = "counts"))
scLung.data <- scLung.data[, featureList.filtered$V1]

## Filter outlier genes in single-cell and in Bulk
#Gene expressed at high magnitude, such as ribosomal protein genes and mitochondrial genes, may dominate the distribution and bias the inference. 
#These genes are often not informative in distinguishing cell types and can be a source of large spurious variance. 
#As a result, they can be detrimental to deconvolution. We recommend the removal of these genes.

sc.stat <- plot.scRNA.outlier(
  input=scLung.data %>% as.matrix(), 
  cell.type.labels=Idents(processed_data),
  species="hs",
  return.raw=TRUE, #return the data used for plotting. 
  pdf.prefix="sc.stat" #specify pdf.prefix if need to output to pdf
)
#Genes like ribossomal proteins show high expression and low expression specificity (check the pdf)

# Remove the genes from selected groups. Also lowly expressed genes (I set min cells to 50, since this is the min required number of cell we need to be considered)
min_cells <- 50
sc.dat.filtered <- cleanup.genes (input=scLung.data %>% as.matrix(),
                                  input.type="count.matrix",
                                  species="hs", 
                                  gene.group=c("Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                  exp.cells=50)

# Check plot after filtering
sc.stat.filtered <- plot.scRNA.outlier(
  input=sc.dat.filtered, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels=Idents(processed_data),
  species="hs",
  return.raw=TRUE, #return the data used for plotting. 
  pdf.prefix="sc.stat.filtered" #specify pdf.prefix if need to output to pdf
)




#Next, we check the concordance of gene expression for different types of genes. As bulk and single cell data are usually collected by different experimental 
#protocols, they may have different sensitivity to different types of genes.
exp.lung.t <- exp.lung %>% 
  mutate(Name = gsub("\\.\\d+", "",  Name)) %>%
  column_to_rownames("Name") %>% 
  t()

# Check concordance
#plot.bulk.vs.sc(sc.input = sc.dat.filtered,
#                 bulk.input = exp.lung.t,
#                 pdf.prefix="bk.vs.sc") 


# Set subcelltypes 
cellTypeMetadata <- data.frame(as.character(Idents(processed_data))) #subtypes
colnames(cellTypeMetadata) <- "subCelltype"

## Merge celltype
processed_data <- RenameIdents(processed_data , `CD8 M/E Tcell-1` = "CD8 M/E Tcell")
processed_data <- RenameIdents(processed_data, `CD8 M/E Tcell-2` = "CD8 M/E Tcell")
processed_data <- RenameIdents(processed_data, `Macrophage-1` = "Macrophage")
processed_data <- RenameIdents(processed_data, `Macrophage-2` = "Macrophage")
processed_data <- RenameIdents(processed_data, `AT2-1` = "AT2")
processed_data <- RenameIdents(processed_data, `AT2-2` = "AT2")


cellTypeMetadata$Celltype <- as.character(unname(Idents(processed_data)))



diff.exp.stat <- get.exp.stat(sc.dat=sc.dat.filtered[,colSums(sc.dat.filtered>0)>3],# filter genes to reduce memory use
                              cell.type.labels=cellTypeMetadata$Celltype,
                              cell.state.labels=cellTypeMetadata$subCelltype,
                              psuedo.count=0.1, #a numeric value used for log2 transformation. =0.1 for 10x data
                              cell.count.cutoff=50, # a numeric value to exclude cell state with number of cells fewer than this value for t test. Default=50.
                              n.cores=6 #number of threads
)

sc.dat.filtered.pc.sig <- select.marker (sc.dat=sc.dat.filtered,
                                         stat=diff.exp.stat,
                                         pval.max=0.01,
                                         lfc.min=0.1)

## Run BayesPrism

Prism <- new.prism(
  reference=sc.dat.filtered.pc.sig, 
  mixture=exp.lung.t,
  input.type="count.matrix", 
  key = NULL,
  cell.type.labels = cellTypeMetadata$Celltype , 
  cell.state.labels = cellTypeMetadata$subCelltype,
  outlier.cut=0.01,
  outlier.fraction=0.1,
)

gtex_lung.ted <- run.prism(prism = Prism, n.cores=6)

#Save data
saveRDS(gtex_lung.ted, file  = paste0(out.folder, "/bayesPrism.lung.signature_genes.rds"))

#Write lines
writeLines(capture.output(sessionInfo()), "GtexTedCountDeconvolutionSessionInfo.txt")