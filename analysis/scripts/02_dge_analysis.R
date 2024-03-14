suppressPackageStartupMessages(library(tidyverse))

degs_analysis <- readRDS(file = "output/differential_expression_results.rds")
metadata <- readRDS(file = "output/metadata/metadata.rds")

get_degs_tables_per_tissue <- function(dge_analysis, analysis_name, fdr = 0.05){
  degs_per_tissue <- c()
  for (tissue in names(dge_analysis)){
    degs <- dge_analysis[[tissue]][[analysis_name]] %>% 
      filter(adj.P.Val < fdr) %>% 
      mutate(tissue = tissue)
    degs_per_tissue <- rbind(degs_per_tissue, degs)
  }
  degs_per_tissue %>% mutate(direction = ifelse(degs_per_tissue$logFC > 0, "Upregulated", "Downregulated"))
}


DEGS <- get_degs_tables_per_tissue(dge_analysis,  analysis_name = "SmokingSMOKER-NEVER", fdr = 0.05)
samples_per_tissue <- sapply(metadata, function(x) nrow(x))
genes_per_tissue <- DEGS %>% group_by(tissue) %>% summarise(n = n())
genes_per_tissue$Sample <- unname(samples_per_tissue[genes_per_tissue$tissue])
cor.test(genes_per_tissue$n, genes_per_tissue$Sample, method = "spearman")

degs.per.tissue <- DEGS %>% 
  group_by(tissue, direction) %>% 
  summarise(n_degs = n())

fdr <- 0.05

degs.per.tissue.wider <- t(sapply(names(dge_analysis), function(tissue)
  sapply("SmokingSMOKER-NEVER", function(trait)
         c(
           sum(dge_analysis[[tissue]][[trait]]$adj.P.Val < fdr & dge_analysis[[tissue]][[trait]]$logFC > 0), 
           sum(dge_analysis[[tissue]][[trait]]$adj.P.Val < fdr & dge_analysis[[tissue]][[trait]]$logFC < 0), 
           sum(dge_analysis[[tissue]][[trait]]$adj.P.Val < fdr)
    )
  )
))


colnames(degs.per.tissue.wider) <- c("Up", "Down", "Total")


tissue_per_gene <- DEGS %>% 
  mutate(ensembl = gsub("\\.\\d+", "", row.names(.))) %>%
  group_by(ensembl, direction) %>% 
  summarise(n_tissue = n()) %>%
  arrange(-n_tissue)


sum(tissue_per_gene$n_tissue == 1) / nrow(tissue_per_gene)


binom_test_up_vs_down <- function(degs.per.tissue.wider){
  #DEGs per tissue in a dataframe with at least 2 colums (Up and Down) and the tissue name set as the row.names
  res <- c()
  for (i in 1:nrow(degs.per.tissue.wider)){
    tissue <- row.names(degs.per.tissue.wider)[i]
    up_tissue <- degs.per.tissue.wider[i, "Up"]
    down_tissue <- degs.per.tissue.wider[i, "Down"]
    total <- up_tissue+ down_tissue
    if (total == 0){
      next
    }
    test <- binom.test(up_tissue, total, p = 0.5)$p.value
    res <- rbind(res, c(tissue, up_tissue, down_tissue, test))
  }
  colnames(res) <- c("tissue","Up", "Down", "p.value")
  res <- res %>% as.data.frame() %>% mutate_at(2:4, as.numeric) %>% mutate(FDR = p.adjust(p.value, method = "fdr"))
}

binom_test_per_tissue <- binom_test_up_vs_down(degs.per.tissue.wider)

DEGS.ExvsNever <- get_degs_tables_per_tissue(dge_analysis,  analysis_name = "SmokingEX-NEVER", fdr = 0.05)


## Summarise degs (per tissue basis)

degs.ex_never.per.tissue <- DEGS.ExvsNever %>% 
  group_by(tissue, direction) %>% 
  summarise(n_degs = n())


### Summarise DEGS (4 colum tables)

fdr <- 0.05

degs.ex_never.per.tissue.wider <- t(sapply(names(dge_analysis), function(tissue)
  sapply("SmokingEX-NEVER", function(trait)
         c(
           sum(dge_analysis[[tissue]][[trait]]$adj.P.Val < fdr & dge_analysis[[tissue]][[trait]]$logFC > 0), 
           sum(dge_analysis[[tissue]][[trait]]$adj.P.Val < fdr & dge_analysis[[tissue]][[trait]]$logFC < 0), 
           sum(dge_analysis[[tissue]][[trait]]$adj.P.Val < fdr)
    )
  )
))


colnames(degs.per.tissue.wider) <- c("Up", "Down", "Total")


### Recurrent genes
tissue_per_gene.ex_vs_never <- DEGS.ExvsNever %>% 
  group_by(gene_name) %>% 
  summarise(n_tissue = n()) %>%
  arrange(-n_tissue)

DEGS.ExvsS <- get_degs_tables_per_tissue(dge_analysis,  analysis_name = "SmokingEX-SMOKER", fdr = 0.05)


## Summarise degs (per tissue basis)

degs.ex_smoker.per.tissue <- DEGS.ExvsS %>% 
  group_by(tissue, direction) %>% 
  summarise(n_degs = n())


### Summarise DEGS (4 colum tables)
fdr <- 0.05

degs.ex_smoker.per.tissue.wider <- t(sapply(names(dge_analysis), function(tissue)
  sapply("SmokingEX-SMOKER", function(trait)
         c(
           sum(dge_analysis[[tissue]][[trait]]$adj.P.Val < fdr & dge_analysis[[tissue]][[trait]]$logFC > 0), 
           sum(dge_analysis[[tissue]][[trait]]$adj.P.Val < fdr & dge_analysis[[tissue]][[trait]]$logFC < 0), 
           sum(dge_analysis[[tissue]][[trait]]$adj.P.Val < fdr)
    )
  )
))


colnames(degs.ex_smoker.per.tissue.wider) <- c("Up", "Down", "Total")


### Recurrent genes
tissue_per_gene.ex_vs_smoker <- DEGS.ExvsS %>% 
  group_by(gene_name) %>% 
  summarise(n_tissue = n()) %>%
  arrange(-n_tissue)



degs_summary <- list()

degs_summary$DEGS <- DEGS
degs_summary$DEGS.per.tissue <- degs.per.tissue.wider
degs_summary$tissue.per.gene <- tissue_per_gene

degs_summary$DEGS.per.tissue_ex_ns <- degs.ex_never.per.tissue.wider
degs_summary$tissue.per.gene_ex_ns <- tissue_per_gene.ex_vs_never

degs_summary$DEGS.per.tissue_ex_s <- degs.ex_smoker.per.tissue.wider
degs_summary$tissue.per.gene_ex_s <- tissue_per_gene.ex_vs_smoker


degs_summary$metadata <- metadata

saveRDS(degs_summary, "../figures/data/degs_summary.rds")