## Preprocess the single cell data for analysis
library(Seurat)
suppressPackageStartupMessages(library(tidyverse))


#Load the dataset
processed_data <- readRDS("data/public/scRNASEQ/GSE173896_COPD.rds") #Downloaded from the papper

#Row.names are gene symbol and I prefere ensembl genes id
#Load ensembl gene set
featureList <- read.table(file = "data/public/scRNASEQ/featureList.txt", sep = "\t")

#Rename the genes in Seurat 
RenameGenesSeurat <- function(obj, newnames) { 
  RNA <- obj@assays$RNA
  
  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
    if (length(RNA@scale.data)) row.names(RNA@scale.data)    <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}

processed_data <- RenameGenesSeurat(processed_data, featureList$V1)


#Load GTEx gene biotype
biotype <- read.delim("data/public/gencode.v26.GRCh38.genes.bed")
biotype <- biotype %>% mutate(ensembl = gsub("\\.\\d+", "", ensembl))

#Keep only genes that are also in the biotype table (since all the analysis involve comparisons with the lung bulk data)
featureList <- featureList[featureList$V1 %in% biotype$ensembl, ]
featureList <- merge(featureList, biotype, by.x = "V1", by.y = "ensembl")

## Preprocess metadata
metadata <- processed_data@meta.data
#Add celltype information
celltype_ident <- Idents(processed_data)
metadata$Celltype <- celltype_ident[row.names(metadata)]

##Add a new colummn for COPD and protocol, filter ex-smokers from the dataset (information from Supplemental Table1：Patient demographics and clinical data and personal correspondence)
ex_smokers <- c("JK01", "JK08", "JK03")
protocol1 <- c("JK01", "JK02", "JK03", "JK04")

metadata <- metadata %>%
  mutate(COPD = ifelse(Class == "COPD", "Yes", "No")) %>%
  mutate(Protocol = ifelse(orig.ident %in% protocol1, "A", "B"))  %>%
  mutate(Smoker = ifelse(Class == "COPD" | Class == "Smoker", "Smoker", "Non-Smoker")) %>%
  mutate(Smoker = ifelse(orig.ident %in%  ex_smokers, "Ex-Smoker", Smoker))

metadata$COPD <- as.factor(metadata$COPD)
metadata$Smoker <- as.factor(metadata$Smoker)

processed_data <- AddMetaData(processed_data, metadata)

## Save data 
saveRDS(featureList, file = "output/featureList_GTEx.rds")
saveRDS(processed_data, file = "output/scRNA_seq_processed.rds")


##Data was log normalized with NormalizeData(COPD, normalization.method = "LogNormalize", scale.factor = 10000) and scaled (by regression the percent.mt)
##Assuming lognorm data is in the "data slot"