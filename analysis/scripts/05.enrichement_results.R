suppressPackageStartupMessages(library(tidyverse))
library(UpSetR)
library(DOSE)


enrichement_NS_S <- readRDS(file = "output/DEGS_enrichement/never_vs_smoker_enrichement.rds") # "Raw enrichment" analysis
enrichement_EX_S <- readRDS(file = "output/DEGS_enrichement/ex_vs_smoker_enrichement.rds") # "Raw enrichment" analysis

abrev <- read.csv(file = "../../metadata/tissue_abreviation.txt", sep = "\t")



get_all_terms <- function(data, tissues_to_plot, slots, FDR, minCount, abreviatures){
  data_to_plot <- c()
  for (slot in slots){ 
    for (tissue in tissues_to_plot){ 
      if (!is.null(data[[tissue]]) & !is.null(data[[tissue]][[slot]])){
        tissue.v2 <- abreviatures %>% filter(SMTSDNOSEP == tissue) %>% pull(SMTSCUSTOM)
        temp_data <- data[[tissue]][[slot]]@result %>% filter(p.adjust < 0.05 & Count > minCount)
        if(nrow(temp_data) > 0){
          if (is.null(data_to_plot[[tissue]])){
            data_to_plot[[tissue.v2]] <- temp_data$ID
          }else{
            data_to_plot[[tissue.v2]] <- c(data_to_plot[[tissue.v2]], temp_data$ID)
          }
          
        }
      }
    }  
  }
  data_to_plot
}

# Get top terms
go_bp_enriched_terms.up <- get_all_terms(enrichement_NS_S, 
                                         names(enrichement_NS_S), 
                                         c("BP_up"), 
                                         0.05, 5, 
                                         abrev)

go_bp_enriched_terms.up.list = fromList(go_bp_enriched_terms.up)


UpSetR::upset(go_bp_enriched_terms.up.list, nsets = ncol(go_bp_enriched_terms.up.list), 
              nintersects = NA, 
              mainbar.y.label = "Common Enriched GO Terms (UP)",
              sets.x.label = "Number of Enriched GO Terms (UP)", 
              text.scale = 2.5, 
              main.bar.color = "black", sets.bar.color = "black")


go_bp_enriched_terms.down <- get_all_terms(enrichement_NS_S, 
                                         names(enrichement_NS_S), 
                                         c("BP_down"), 
                                         0.05, 5, 
                                         abrev)

go_bp_enriched_terms.down.list = fromList(go_bp_enriched_terms.down)


UpSetR::upset(go_bp_enriched_terms.down.list, nsets = ncol(go_bp_enriched_terms.down.list), 
              nintersects = NA, 
              mainbar.y.label = "Common Enriched GO Terms (Down)",
              sets.x.label = "Number of Enriched GO Terms (Down)", 
              text.scale = 2.5, 
              main.bar.color = "black", sets.bar.color = "black")


terms.per.tissue.up <- data.frame(tissue = rep(names(go_bp_enriched_terms.up), sapply(go_bp_enriched_terms.up, length)),
                 term = unlist(go_bp_enriched_terms.up))

tissue.per.terms.up <- terms.per.tissue.up %>% 
  group_by(term) %>%
  summarise(n_tissues = n()) %>% 
  arrange(-n_tissues)

terms.per.tissue.down <- data.frame(tissue = rep(names(go_bp_enriched_terms.down), sapply(go_bp_enriched_terms.down, length)),
                 term = unlist(go_bp_enriched_terms.down))

tissue.per.terms.down <- terms.per.tissue.down %>% 
  group_by(term) %>%
  summarise(n_tissues = n()) %>% 
  arrange(-n_tissues)

#Extract terms from the enrichment analysis (KEGG)
kegg_enriched_terms <- get_all_terms(enrichement_NS_S, 
                                      names(enrichement_NS_S), 
                                      c("kegg_up", "kegg_down"), 
                                      0.05, 5, 
                                      abrev)

kegg_enriched_terms.list = fromList(kegg_enriched_terms)


UpSetR::upset(kegg_enriched_terms.list, nsets = ncol(kegg_enriched_terms.list), 
              nintersects = NA, 
              mainbar.y.label = "Common Enriched KEGG ",
              sets.x.label = "Number of Enriched KEGG Pathways", 
              text.scale = 2.5, 
              main.bar.color = "black", sets.bar.color = "black")

##Functions
get_top_enrichement_results_up_down <- function(data, 
                                             tissues_to_plot, 
                                             slot_up, 
                                             slot_down, 
                                             terms_per_tissue, 
                                             FDR = 0.05, 
                                             minCount = 5){
  
  data_to_plot <- c()
  terms_to_plot <- c()
  
  #Get list of terms in slot up
  for (tissue in tissues_to_plot){ 
    if (!is.null(data[[tissue]]) & !is.null(data[[tissue]][[slot_up]])){
      temp_terms <- data[[tissue]][[slot_up]]@result %>% 
        filter(p.adjust < FDR & Count > minCount) %>% 
        slice_head(n = terms_per_tissue)
      if  (nrow(temp_terms) > 0){
        terms_to_plot <- c(terms_to_plot, temp_terms$Description)
      }
    }
  }
  
  #Get list of terms in slot down
  for (tissue in tissues_to_plot){ 
    if (!is.null(data[[tissue]]) & !is.null(data[[tissue]][[slot_down]])){
      temp_terms <- data[[tissue]][[slot_down]]@result %>% 
        mutate(qvalue = ifelse(is.na(qvalue), 1, qvalue)) %>% 
        filter(p.adjust < FDR & Count > minCount) %>% 
        slice_head(n = terms_per_tissue)
      if  (nrow(temp_terms) > 0){
        terms_to_plot <- c(terms_to_plot, temp_terms$Description)
      }
    }
  }
  
  #Get all selected terms in all tissues in slot up
  for (tissue in tissues_to_plot){ 
    if (!is.null(data[[tissue]]) & !is.null(data[[tissue]][[slot_up]])){
      temp_data <- data[[tissue]][[slot_up]]@result %>% 
        mutate(qvalue = ifelse(is.na(qvalue), 1, qvalue)) %>% 
        filter(p.adjust < FDR & Count > minCount & Description %in% terms_to_plot)
      if(nrow(temp_data) > 0){
        temp_data <- temp_data %>% 
          mutate(sample= tissue) %>%
          mutate(Direction = "Upregulated")
        data_to_plot <- rbind(data_to_plot, temp_data)
      }
    }
  }
  
  #Get all selected terms in all tissue in slot down
  for (tissue in tissues_to_plot){ 
    if (!is.null(data[[tissue]]) & !is.null(data[[tissue]][[slot_down]])){
      temp_data <- data[[tissue]][[slot_down]]@result %>% 
        mutate(qvalue = ifelse(is.na(qvalue), 1, qvalue)) %>% 
        filter(p.adjust < FDR & qvalue < 0.2 & Count > minCount & Description %in% terms_to_plot)
      if(nrow(temp_data) > 0){
        temp_data <- temp_data %>% 
          mutate(sample= tissue) %>%
          mutate(Direction = "Downregulated")
        data_to_plot <- rbind(data_to_plot, temp_data)
      }
    }
  }
  
  data_to_plot <- data.frame(data_to_plot) 
  return (data_to_plot)
}

## Plot Orsum results
orsum_bp.all.up <- read_tsv("output/DEGS_enrichement/BP.up/filteredResult-Summary.tsv")
orsum_bp.all.down <- read_tsv("output/DEGS_enrichement/BP.down/filteredResult-Summary.tsv")

# Parse ORSUM results to a wider table
## Note: Not all terms represented here are enriched in the specific tissues if I take a look at the individual analysis due to the way the terms are summarised (the GO terms from all tissues are pooled together)

## Get top terms from the summarise term set
#I defined the top term as terms with "Representing term rank of 6 or lower". Term rank is defined as the lower rank a term has in any tissue

max_rank <- 6

orsum_bp.all.up.top_terms <- orsum_bp.all.up %>% 
  filter(`Representing term rank` <= max_rank) %>% 
  mutate(direction = "Upregulated")

orsum_bp.all.down.top_terms <- orsum_bp.all.down %>% 
  filter(`Representing term rank` <= max_rank)  %>%
  mutate(direction = "Downregulated")

# Add missing collums (the tissues are not the same)
collums <- union(colnames(orsum_bp.all.up.top_terms), colnames(orsum_bp.all.down.top_terms))

for(colname in collums){
  if (!colname %in% colnames(orsum_bp.all.up.top_terms)){
    orsum_bp.all.up.top_terms[[colname]] <- NA
  }
  
  if (!colname %in% colnames(orsum_bp.all.down.top_terms)){
    orsum_bp.all.down.top_terms[[colname]] <- NA
  }
  
}

orsum_bp.all.top_terms <- rbind(orsum_bp.all.up.top_terms, orsum_bp.all.down.top_terms)
tissues <- colnames(orsum_bp.all.top_terms)[grep("^(?!.*Representing).*\\bterm\ rank\\b", colnames(orsum_bp.all.top_terms), perl = T)]

orsum_bp.all.top_terms <- orsum_bp.all.top_terms %>%
  pivot_longer(cols = all_of(tissues)) %>% 
  filter(value != "None") %>% 
  mutate(tissue = gsub(" term rank", "", name)) %>%
  mutate(direction = factor(direction, levels = c("Upregulated", "Downregulated"))) %>% 
  mutate(`Representing term name` = factor(`Representing term name`, levels = unique(`Representing term name`[nrow(.):1]))) %>% 
  merge(abrev, by.x = "tissue", by.y = "SMTSDNOSEP")



ggplot(orsum_bp.all.top_terms, aes(x = SMTSCUSTOM, y = `Representing term name`, color = direction)) +
  geom_point(size = 4) +
  scale_color_manual(name = "Gene Set", values =  c("tomato3", "royalblue")) +
  theme_dose(10) +
  xlab("") + ylab("") + 
  theme(axis.text.x = element_text(size = 16, angle = 90, hjust = 1, vjust = 0.1), 
        axis.text.y = element_text(size = 18), 
        legend.title = element_text(size = 25), 
        legend.text = element_text(size = 22), 
        legend.position = "top", 
        plot.margin = margin(l = 100))



## ----fig.height=14, fig.width=22-------------------------------------------------------------------------------------------------------------------------
kegg_top_terms <- get_top_enrichement_results_up_down(enrichement_NS_S, names(enrichement_NS_S), 
                                    "kegg_up", "kegg_down", 10, 0.05, 5)


#kegg_top_terms$ID[duplicated(kegg_top_terms$ID)] ##Check for duplicated terms  (and if the direction matches)


kegg_top_terms %>% 
  mutate(Direction = factor(Direction, levels = c("Upregulated", "Downregulated"))) %>%
  merge(abrev, by.x = "sample", by.y = "SMTSDNOSEP") %>%
  arrange(SMTSCUSTOM) %>% 
  mutate(Description = factor(Description, levels = unique(Description)[length(unique(Description)):1])) %>% 
  ggplot(aes(x=SMTSCUSTOM, y=Description, color = Direction)) + 
  geom_point(size = 10) +
  scale_color_manual(name = "Gene Set", values = c("tomato3", "royalblue")) +
  theme_dose(10) +
  scale_size(name = "Gene Number", range=c(2, 8)) +
  xlab("") + ylab("") + 
  theme(axis.text.x = element_text(size = 16, angle = 90, hjust = 1, vjust = 0.1), 
        axis.text.y = element_text(size = 18), 
        legend.title = element_text(size = 22), 
        legend.text = element_text(size = 20), 
        legend.position = "top", 
        plot.margin = margin(l = 200), 
        element_line(linewidth =1.5))


lung_DO_terms_up <- enrichement_NS_S$Lung$do_up@result %>% filter(p.adjust < 0.05, Count > 5) 
thyroid_DO_top_terms <- enrichement_NS_S$Thyroid$do_down@result %>% filter( p.adjust < 0.05, Count > 5) %>% 
  slice_head(n = 15)

ggplot(thyroid_DO_top_terms, aes(y = Description, x = Count)) + 
  geom_bar(stat = "identity") + 
  xlab("") + ylab("") + 
  theme_classic()



## Run enrichment analysis on common genes between Lung and Thyroid, pacnreas and Skin 
library(clusterProfiler)
library(org.Hs.eg.db)


degs_analysis <- readRDS(file = "output/differential_expression_results.rds")

run_enrichement_on_common_genes <- function(degs_analysis, tissues, n_tissue_in_common, direction){
  #@tissues tissues to consider
  #@n_tissue_in_common Number of tissues a gene has to appear to be considered in common. Used to run the enrichement
  #@direction - Direction of the change (Up, Down or both)
  
  gene_list <- c()
  background_common <- c()
  for (tissue in tissues){
    if (direction == "Up"){
      gene_tissue <- row.names(degs_analysis[[tissue]][["SmokingSMOKER-NEVER"]] %>% filter(adj.P.Val < 0.05 & logFC > 0))
    }
    #else if (direction == "Down")
    #else if (direction == "both")
    gene_tissue <- data.frame("gene" = gene_tissue, "tissue" = tissue)
    gene_list <- rbind(gene_list, gene_tissue)
    background_union <- unique(c(background_common, row.names(degs_analysis[[tissue]][["SmokingSMOKER-NEVER"]])))
  }
  
  #Filter the gene list 
  gene_list_to_run <- gene_list %>% 
    group_by(gene) %>% 
    summarise(n = n()) %>% 
    filter(n >= n_tissue_in_common)
  
  cat("Running enrichement with ",  nrow(gene_list_to_run))
  
  # Run enrichment
  genes <- gsub("\\.\\d+", "", gene_list_to_run$gene)
  background <- gsub("\\.\\d+", "", background_union)
  
  
  enrichGO(gene = genes, keyType = "ENSEMBL", OrgDb='org.Hs.eg.db', pvalueCutoff = 0.05, ont = "BP", pAdjustMethod = "BH", universe = background)
}

# Call function
enrichement_common_genes_3_tissues <- run_enrichement_on_common_genes(degs_analysis, tissues = c("Lung", "Thyroid","Pancreas", "SkinSunExposedLowerleg"), 3, direction = "Up")

enrichement_common_genes_4_tissues <- run_enrichement_on_common_genes(degs_analysis, tissues = c("Lung", "Thyroid","Pancreas", "SkinSunExposedLowerleg"), 4, direction = "Up")

enrichement_common_genes_3_tissues_no_lung <- run_enrichement_on_common_genes(degs_analysis, tissues = c("Thyroid","Pancreas", "SkinSunExposedLowerleg"), 3, direction = "Up")

enrichement_common_genes_lung_thyroid <- run_enrichement_on_common_genes(degs_analysis, tissues = c("Lung", "Thyroid"), 2, direction = "Up")

# Common terms among Lung, Thyroid, Pacnreas and Skin Sun exposed
common_terms <- intersect(go_bp_enriched_terms.up$Lung, intersect(go_bp_enriched_terms.up$Pancreas, intersect(go_bp_enriched_terms.up$Thyroid, go_bp_enriched_terms.up$`Skin - Sun Exposed`)))

enrichement_NS_S$Lung$BP_up@result %>% filter(ID %in% common_terms)


# In enrichement_common_genes_4_tissues we see terms like positive thymic T cell selection enriched (not the exact same terms but they are related)

data <- list()
data$enrichement_ns_s_go_bp.up <- go_bp_enriched_terms.up.list
data$enrichement_ns_s_go_bp.down <- go_bp_enriched_terms.down.list
data$enrichement_ns_s_kegg <- kegg_enriched_terms.list
data$tissue_per_terms.up <- tissue.per.terms.up
data$tissue_per_terms.down <- tissue.per.terms.down
data$orsum_top_terms <- orsum_bp.all.top_terms
data$KEGG_top_terms <- kegg_top_terms
data$thyroid_do_down <- thyroid_DO_top_terms
data$lung_do_up <- lung_DO_terms_up


saveRDS(data, file = "../figures/data/enrichement_results.summary.rds")