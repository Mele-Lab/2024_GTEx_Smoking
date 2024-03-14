#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez and Rogério Ribeiro
# @E-mail: jose.ramirez1@bsc.es and rogerio.e.ramos.ribeiro@gmail.com
# @Description: Code to compare with literature (transcriptome)
# @software version: R=4.2.2


suppressPackageStartupMessages(library(tidyverse))
library(Seurat)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggvenn)

smoking_degs <- readRDS(file = "output/differential_expression_results.rds")


data <- list()
data$reproduced_data <- c() #("Reproduced", "Total", "Study name", "tissue", "OR")

#### Blood

gtex_whole_blood <- smoking_degs$WholeBlood$`SmokingSMOKER-NEVER` %>% filter(adj.P.Val < 0.05)
gtex_whole_blood_expressed <- smoking_degs$WholeBlood$`SmokingSMOKER-NEVER`


row.names(gtex_whole_blood) <- gsub("\\.\\d+","",row.names(gtex_whole_blood))
row.names(gtex_whole_blood_expressed) <- gsub("\\.\\d+","",row.names(gtex_whole_blood_expressed))

gtex_whole_blood$gene_dir <- paste0(row.names(gtex_whole_blood), "_",ifelse(gtex_whole_blood$logFC > 0, "+", "-"))

## Huan 2016
huan_2016 <- read_tsv(file = "public/related_report/huan2016/huan_2016_never_vs_smoker_degs.tsv")

#Give esembl_id to the table
huan_2016$esembl_id <-  mapIds(org.Hs.eg.db,
                                  keys=as.character(huan_2016$Entrez.Gene.ID), 
                                  keytype ="ENTREZID",
                                  column ="ENSEMBL",
                                  multiVals="first")


#Note: Some Id did not map
sum(is.na(huan_2016$esembl_id)) #113 are missing.

huan_2016 <- huan_2016 %>% 
  mutate(esembl_id_dir = paste0(esembl_id, "_", ifelse(Meta.Beta > 0, "+", "-")))

# From the paper:
# We at first identified differentially expressed genes for smoking in the discovery set at FDR<0.1, and

gtex_whole_blood.fdr0_1 <- gtex_whole_blood %>% filter(adj.P.Val < 0.1)

#How many genes DEG in common with GTEx study? 
gtex_whole_blood_intersect_huan <- merge(gtex_whole_blood.fdr0_1, huan_2016, by.x = 0, by.y = "esembl_id")

#any(duplicated(row.names(gtex_whole_blood_intersect_huan)))
#any(duplicated(gtex_whole_blood_intersect_huan$GeneSymbol))
#any(duplicated(gtex_whole_blood_intersect_huan$esembl_id_dir))

# Subset genes also expressed in the GTEx data
huan_2016_no_na <- huan_2016 %>% 
  filter(!is.na(esembl_id)) %>% 
  filter(esembl_id %in% row.names(gtex_whole_blood_expressed))


huan_2016_venn_data <- list(`GTEx\nWholeBlood` = gtex_whole_blood.fdr0_1$gene_dir, 
                            `Huan 2016\nFDR 0.1` = huan_2016_no_na$esembl_id_dir)

huan_2016_venn_01 <- ggvenn::ggvenn(huan_2016_venn_data)


#Plot results
print(huan_2016_venn_01)


# Fisher test
sig_huan_sig_GTEx <- nrow(gtex_whole_blood_intersect_huan)
sig_huan_non_sig_GTEx <- nrow(huan_2016_no_na) - nrow(gtex_whole_blood_intersect_huan)
non_sig_huan_sig_GTEx <- nrow(gtex_whole_blood) - nrow(gtex_whole_blood_intersect_huan)
non_sig_huan_gtex <- nrow(gtex_whole_blood_expressed) - sig_huan_sig_GTEx - sig_huan_non_sig_GTEx - non_sig_huan_sig_GTEx


fisher_matrix <- matrix(c(sig_huan_sig_GTEx, sig_huan_non_sig_GTEx,
                   non_sig_huan_sig_GTEx, non_sig_huan_gtex),
                   nrow = 2, ncol = 2, byrow = T)

huan_gtex_01_fisher <- fisher.test(fisher_matrix)

huan_gtex_01_fisher


data$reproduced_data <- rbind(data$reproduced_data, 
                              c(sig_huan_sig_GTEx, nrow(huan_2016_no_na), "Huan et al.", "Whole Blood", huan_gtex_01_fisher$estimate, huan_gtex_01_fisher$conf.int[[1]], huan_gtex_01_fisher$conf.int[[2]]))



## Vink et al
vink_2017 <- read_tsv(file = "../../related_report/vink2017/degs.tsv") ## All these genes are DEG between never and smokers 
#write.csv(vink_2017$GeneSymbol, file = "../../related_report/vink2017/temp_gene_symbols.tsv", quote = F, row.names = F, col.names = F)

#Submit gene name here to update https://www.genenames.org/tools/multi-symbol-checker/
updated_symbols_vink2017 <- read_tsv(file = "data/public/related_report/vink2017/update_smbols.txt")
updated_symbols_vink2017_ <- updated_symbols_vink2017$`Approved symbol`
names(updated_symbols_vink2017_) <- updated_symbols_vink2017$Input

vink_2017$GeneSymbolUpdated <- ifelse(vink_2017$GeneSymbol %in% updated_symbols_vink2017$Input, 
                                      updated_symbols_vink2017_[vink_2017$GeneSymbol], vink_2017$GeneSymbol)

vink_2017$esembl_id <-  mapIds(org.Hs.eg.db,
                                        keys=as.character(vink_2017$GeneSymbolUpdated), 
                                        keytype ="SYMBOL",
                                        column ="ENSEMBL",
                                        multiVals="first")

#Add a missing esembl 
vink_2017[7, "esembl_id"] <- "ENSG00000211829"

#Add direction information
vink_2017$esembl_id_dir <- paste0(vink_2017$esembl_id, "_", ifelse(vink_2017$`Beta_C vs N` > 0, "-", "+")) # Beta is reversed, "-" are upregulated in smoker

#Subet col and get only genes expressed in GTEx blood samples 
vink_2017_never_vs_smoker <- vink_2017[c(32,2,3, 15,31,33)] %>% 
  filter(esembl_id %in% row.names(gtex_whole_blood_expressed))


#Compare with never vs smokers in GTEx Whole blood 

# From the Paper
#To correct for multiple testing, a Bonferroni correction was applied (P<0.05/44 241 =1206).//Bonferroni 0.05
gtex_whole_blood_bon <- smoking_degs$WholeBlood$`SmokingSMOKER-NEVER` %>% 
  mutate(adj.P.Val = p.adjust(`P.Value`, method = "bonferroni")) %>% 
  filter(adj.P.Val < 0.05)


row.names(gtex_whole_blood_bon) <- gsub("\\.\\d+","",row.names(gtex_whole_blood_bon))
gtex_whole_blood_bon$gene_dir <- paste0(row.names(gtex_whole_blood_bon), "_",ifelse(gtex_whole_blood_bon$logFC > 0, "+", "-"))


vink_2017_gtex_never_vs_smk_data <- list("GTEx\nWhole Blood" = gtex_whole_blood_bon$gene_dir, 
                                         "Vink 2017" = vink_2017_never_vs_smoker$esembl_id_dir)


ggvenn(vink_2017_gtex_never_vs_smk_data)

## Fisher test for overlap between gene sets 
# No gene background was available Assume the background are similar

fisher_matrix <- matrix(c(length(intersect(vink_2017_never_vs_smoker$esembl_id_dir, gtex_whole_blood_bon$gene_dir)),
                          length(setdiff(vink_2017_never_vs_smoker$esembl_id_dir,  gtex_whole_blood_bon$gene_dir)),
                          length(setdiff(gtex_whole_blood_bon$gene_dir, vink_2017_never_vs_smoker$esembl_id_dir)),
                          nrow(gtex_whole_blood_expressed) - length(union(gtex_whole_blood_bon$gene_dir, vink_2017_never_vs_smoker$esembl_id_dir))),
                        nrow = 2, ncol = 2, byrow = T)

vink_gtex_fisher <- fisher.test(fisher_matrix)

vink_gtex_fisher


data$reproduced_data <- rbind(data$reproduced_data, 
                              c(length(intersect(vink_2017_never_vs_smoker$esembl_id_dir, gtex_whole_blood$gene_dir)), nrow(vink_2017_never_vs_smoker), "Vink et al.", "Whole Blood", vink_gtex_fisher$estimate, vink_gtex_fisher$conf.int[[1]], vink_gtex_fisher$conf.int[[2]]))


## Lung tissue

gtex_lung_degs <- smoking_degs$Lung$`SmokingSMOKER-NEVER`%>% filter(adj.P.Val < 0.05)
gtex_lung_expression <- smoking_degs$Lung$`SmokingSMOKER-NEVER`

row.names(gtex_lung_degs) <- gsub("\\.\\d+", "", row.names(gtex_lung_degs))
gtex_lung_degs$gene_dir <- paste0(row.names(gtex_lung_degs), "_", ifelse(gtex_lung_degs$logFC < 0, "-", "+"))

row.names(gtex_lung_expression) <- gsub("\\.\\d+", "", row.names(gtex_lung_expression)) 


### Bose 2012
bose_degs <- read.csv(file = "data/public/related_report/bosse2012/data_bose_et_all_2012.csv", sep = ";") %>%
  filter(!GeneSymbol == "") %>%
  filter(!duplicated(GeneSymbol))  %>%
  mutate(adjust_bh = p.adjust(gsub(",", "\\.", wilcox.pvalue)))


bose_degs[320, "GeneSymbol"] <- "MARCHF1"
bose_degs[873, "GeneSymbol"] <- "SEPTIN3"
bose_degs[1822, "GeneSymbol"] <- "SEPTIN10"
# 
# 
# write.csv(bose_degs$GeneSymbol, file = "../../related_report/bosse2012/temp_gene_symbols.tsv", quote = F, row.names = F)

#Submit gene name here to update https://www.genenames.org/tools/multi-symbol-checker/
updated_symbol_bosse_2012 <- read_csv(file = "data/public/related_report/bosse2012/updated_gene_symbols.csv")

updated_symbol_bosse_2012_ <- updated_symbol_bosse_2012$`Approved symbol`
names(updated_symbol_bosse_2012_) <-updated_symbol_bosse_2012$Input

bose_degs$GeneSymbolUpdated <- ifelse(bose_degs$GeneSymbol %in% updated_symbol_bosse_2012$Input, 
                                      updated_symbol_bosse_2012_[bose_degs$GeneSymbol], bose_degs$GeneSymbol)



bose_degs$esembl_id <-  mapIds(org.Hs.eg.db,
                               keys=as.character(bose_degs$GeneSymbolUpdated), 
                               keytype ="SYMBOL",
                               column ="ENSEMBL",
                               multiVals="first")

#ST20-AS1 <- ENSG00000259642 #908
#LINC01560-> ENSG00000196741 #395
#AFDN-DT -> ENSG00000198221 #689

bose_degs$esembl_id[908] <- "ENSG00000259642"
bose_degs$esembl_id[395] <- "ENSG00000196741"
bose_degs$esembl_id[689] <- "AFDN-DT"

bose_degs$esembl_id_gene_dir <- paste0(bose_degs$esembl_id, "_", ifelse(bose_degs$FoldChange > 0,"+", "-"))

bose_degs <- bose_degs %>% 
  filter(esembl_id %in% row.names(gtex_whole_blood_expressed))

#Note on the logFC: fold change in Bosse 2012 was compute as 2 ^ smkExp - nonSmkExp. 
#Hence, fc < 1 the gene is downregualted and fc > 1 the gene is upregulated

sum(is.na(bose_degs$esembl_id)) 


#Intersect with the entire dataset
# From the paper
#We applied a Bonferroni correction to correct for multiple testing (0.05/38,820 probe sets, P value < 1.29  106).

gtex_lung_degs_bon <- smoking_degs$Lung$`SmokingSMOKER-NEVER`%>% mutate(adj.P.Val = p.adjust(`P.Value`, method = "bonferroni")) %>% filter(adj.P.Val < 0.05)

row.names(gtex_lung_degs_bon) <- gsub("\\.\\d+", "", row.names(gtex_lung_degs_bon))
gtex_lung_degs_bon$gene_dir <- paste0(row.names(gtex_lung_degs_bon), "_", ifelse(gtex_lung_degs_bon$logFC < 0, "-", "+"))

# Intersect with bose
bose_gtex_lung_data2 <- list("Bosse 2012" = bose_degs$esembl_id_gene_dir, "GTEx Lung\nFDR" = gtex_lung_degs_bon$gene_dir)
bose_venn_2 <- ggvenn(bose_gtex_lung_data2, text_size = 7)

bose_venn_2


#Fisher test
fisher_matrix <- matrix(c(length(intersect(bose_degs$esembl_id_gene_dir, gtex_lung_degs_bon$gene_dir)),
                          length(setdiff(bose_degs$esembl_id_gene_dir,  gtex_lung_degs_bon$gene_dir)),
                          length(setdiff( gtex_lung_degs_bon$gene_dir, bose_degs$esembl_id_gene_dir)),
                          nrow(smoking_degs$Lung$`SmokingSMOKER-NEVER`) - length(union(gtex_lung_degs_bon$gene_dir, bose_degs$esembl_id_gene_dir))),
                          nrow = 2, ncol = 2, byrow = T)

bose_gtex_bon_fisher <- fisher.test(fisher_matrix)

bose_gtex_bon_fisher


data$reproduced_data <- rbind(data$reproduced_data, 
                              c(length(intersect(bose_degs$esembl_id_gene_dir, gtex_lung_degs$gene_dir)), 
                                nrow(bose_degs), "Bossé et al.", "Lung", bose_gtex_bon_fisher$estimate, bose_gtex_bon_fisher$conf.int[[1]], bose_gtex_bon_fisher$conf.int[[2]]))


### Landi at al 2008
landi_degs <- read_tsv(file = "data/public/related_report/landi2008/non_tumour_samples_landi2008.tsv") 
landi_degs$dir <- c(rep("+", 28), rep("-", 75)) #First 28 genes are up, the other are down

#write.csv(landi_degs$`Gene Symbol`, file = "../../related_report/landi2008/temp_gene_symbols.tsv", row.names = F, quote = F)
#Submit gene name here to update https://www.genenames.org/tools/multi-symbol-checker/

updated_symbol_landi_2008 <- read_csv(file = "../../related_report/landi2008/hgnc-symbol-check.csv")

updated_symbol_landi_2008_ <- updated_symbol_landi_2008$`Approved symbol`
names(updated_symbol_landi_2008_) <-updated_symbol_landi_2008$Input



landi_degs$GeneSymbolUpdated <- ifelse(landi_degs$`Gene Symbol` %in% updated_symbol_landi_2008$Input, 
                                      updated_symbol_landi_2008_[landi_degs$`Gene Symbol`], landi_degs$`Gene Symbol`)

landi_degs$esembl_id <-  mapIds(org.Hs.eg.db,
                               keys=as.character(landi_degs$GeneSymbolUpdated), 
                               keytype ="SYMBOL",
                               column ="ENSEMBL",
                               multiVals="first")

landi_degs$esembl_id_gene_dir <- paste0(landi_degs$esembl_id, "_", landi_degs$dir)


landi_degs <- landi_degs %>% 
  filter(esembl_id %in% row.names(gtex_lung_expression))

## Compare with GTEx
# From Landi et al 2008: 
#To limit false positive findings, genes were considered statistically significant if their p-values were less than the stringent threshold of 0.001. Under the null hypothesis of no difference in expression profiles, and considering the analysis of 22,283 probes, we expect that by chance the average number of false positive findings will be < 23. We used the Benjamini-Hochberg[2] procedure to calculate the False Discovery Rate (FDR). Gene selection based on p,0.001 (two-sided) and foldchange 1.5 are referred to as ‘‘stringent criteria’’.

# Its unclear if the p.values are BH corrected, but I think they are
# Further filter GTEx DEGs by logFC

treshold_fc <- log(1.5, 2)

gtex_lung_degs_0.001_foldChange1_5 <- gtex_lung_degs %>% 
  filter(adj.P.Val < 0.001) %>% 
  filter(abs(logFC) > treshold_fc)

landi_gtex_data <- list(`GTEx Lung` = gtex_lung_degs_0.001_foldChange1_5$gene_dir, "Landi et al. 2008" = landi_degs$esembl_id_gene_dir)
landi_gtex_venn_1 <- ggvenn(landi_gtex_data, text_size =  8)
landi_gtex_venn_1

## Fisher test

fisher_matrix <- matrix(c(length(intersect(landi_degs$esembl_id_gene_dir, gtex_lung_degs_0.001_foldChange1_5$gene_dir)),
                          length(setdiff(landi_degs$esembl_id_gene_dir,  gtex_lung_degs_0.001_foldChange1_5$gene_dir)),
                          length(setdiff( gtex_lung_degs_0.001_foldChange1_5$gene_dir, landi_degs$esembl_id_gene_dir)),
                          nrow(smoking_degs$Lung$`SmokingSMOKER-NEVER`) - length(union(gtex_lung_degs_0.001_foldChange1_5$gene_dir, landi_degs$esembl_id_gene_dir))),
                        nrow = 2, ncol = 2, byrow = T)

landi_gtex_fisher <- fisher.test(fisher_matrix)

data$reproduced_data <- rbind(data$reproduced_data, 
                              c(length(intersect(landi_degs$esembl_id_gene_dir, gtex_lung_degs_0.001_foldChange1_5$gene_dir)), 
                                nrow(landi_degs), "Landi et al.", "Lung", landi_gtex_fisher$estimate, landi_gtex_fisher$conf.int[[1]], landi_gtex_fisher$conf.int[[2]]))


### Pintarelli 2019
pintarelli <- read.csv(file ="data/public/related_report/pintarelli2019/pintarelli_2019_degs.csv")[,-1] %>% 
  filter(Gene.symbol != "")

pintarelli$esembl_id <-  mapIds(org.Hs.eg.db,
                               keys=as.character(pintarelli$Gene.symbol), 
                               keytype ="SYMBOL",
                               column ="ENSEMBL",
                               multiVals="first")

#sum(is.na(pintarelli$esembl_id))

pintarelli[111, "esembl_id"] <- "ENSG00000134698"
pintarelli[118, "esembl_id"] <- "ENSG00000204257"
pintarelli[251, "esembl_id"] <- "ENSG00000173421"

pintarelli[82, "esembl_id"] <- "ENSG00000153485"
pintarelli[174, "esembl_id"] <- "ENSG00000177628"


pintarelli <- pintarelli %>% 
  filter(esembl_id %in% row.names(gtex_lung_expression))
  

#Since we do not have information about the logFC I can not consider it for this test
gtex_lung_degs <- gtex_lung_degs %>% 
  mutate(gene = gsub("\\.\\d+", "", row.names(.)))


pintarelli_gtex_lung_data <- list("Pintarelli 2019" = pintarelli$esembl_id, "GTEx Lung" = gtex_lung_degs$gene)


ggvenn(pintarelli_gtex_lung_data)


# Fisher test
fisher_matrix <- matrix(c(length(intersect(pintarelli$esembl_id, gtex_lung_degs$gene)),
                          length(setdiff(pintarelli$esembl_id,  gtex_lung_degs$gene)),
                          length(setdiff(gtex_lung_degs$gene, pintarelli$esembl_id)),
                          nrow(smoking_degs$Lung$`SmokingSMOKER-NEVER`) - length(union(gtex_lung_degs$gene, pintarelli$esembl_id))),
                        nrow = 2, ncol = 2, byrow = T)

pintarelli_gtex_fisher <- fisher.test(fisher_matrix)
pintarelli_gtex_fisher

data$reproduced_data <- rbind(data$reproduced_data, 
                              c(length(intersect(pintarelli$esembl_id, gtex_lung_degs$gene)), 
                                nrow(pintarelli), "Pintarelli et al.", "Lung", pintarelli_gtex_fisher$estimate, pintarelli_gtex_fisher$conf.int[[1]], pintarelli_gtex_fisher$conf.int[[2]]))



## Adipose tissue (Tsai 2018)
adipose_tissue_subcutenous <- smoking_degs$AdiposeSubcutaneous$`SmokingSMOKER-NEVER` %>% 
  dplyr::filter(adj.P.Val < 0.1) %>% 
  mutate(gene_dir = paste0(gsub("\\.\\d+", "", row.names(.)), "_", ifelse(logFC > 0, "+", "-")))

adipose_tissue_subcutenous_expression <- gsub("\\.\\d+", "", row.names(smoking_degs$AdiposeSubcutaneous$`SmokingSMOKER-NEVER`))


tsai_adipose_degs <- read.csv(file  = "data/public/related_report/tsai2018/degs_adipose_tissue.tsv") %>%
  mutate(name = gsub("\\.\\d+", "", ID)) %>%
  mutate(gene_dir = paste0(gsub("\\.\\d+", "", ID), "_", ifelse(Coef. > 0, "+", "-"))) %>%
  filter(name %in% adipose_tissue_subcutenous_expression)

# From the paper: 
#For each CpG site or gene, a full model, that regressed all of the covariates was compared to a null model that excluded smoking status. The models were compared using the ANOVA F statistic. A genome-wide significance level was set at 1% false discovery rate for all analyses.

tsai_2018_gtex_adipose_venn_data <- list("GTEx Adipose" = adipose_tissue_subcutenous$gene_dir, 
                                         "Tsai et al 2018" = tsai_adipose_degs$gene_dir)

ggvenn(tsai_2018_gtex_adipose_venn_data, text_size =  8)

fisher_matrix <- matrix(c(length(intersect(tsai_adipose_degs$gene_dir, adipose_tissue_subcutenous$gene_dir)),
                          length(setdiff(tsai_adipose_degs$gene_dir,  adipose_tissue_subcutenous$gene_dir)),
                          length(setdiff(adipose_tissue_subcutenous$gene_dir,  tsai_adipose_degs$gene_dir)),
                          length(adipose_tissue_subcutenous_expression) - length(union(adipose_tissue_subcutenous$gene_dir,  tsai_adipose_degs$gene_dir))),
                        nrow = 2, ncol = 2, byrow = T)

tsai_gtex_fisher <- fisher.test(fisher_matrix)
tsai_gtex_fisher

data$reproduced_data <- rbind(data$reproduced_data, 
                              c(length(intersect(tsai_adipose_degs$gene_dir, adipose_tissue_subcutenous$gene_dir)), 
                                nrow(adipose_tissue_subcutenous), "Tsai et al.", "Adipose Tissue", tsai_gtex_fisher$estimate, tsai_gtex_fisher$conf.int[[1]], tsai_gtex_fisher$conf.int[[2]]))


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
saveRDS(data, "../figures/data/pca.rds")

