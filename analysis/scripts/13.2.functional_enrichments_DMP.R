#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Functional enrichment differentially methylated positions/loci
# @software version: R=4.2.2

library(dplyr)
library(tidyr)
library(tidyverse)
library(missMethyl)
library(ggplot2)

#Set path 
# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")


#Functional enrichments

#Reading annotation
Sys.time()
data_path <- "~/Documents/mn4/scratch/bsc83/bsc83535/GTEx/v9/Oliva/"
annotation <- read.csv(paste0(data_path, "GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv"))
# annotation <- read.csv(paste0(first_dir, "/scratch/bsc83/bsc83535/GTEx/v9/Oliva/GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv"))
Sys.time()

#Reading data
results <- readRDS("tissues/Lung/DML_results.rds")
smoking2 <- results$Smoking2
signif <- smoking2[smoking2$adj.P.Val<0.05,]
# results <- readRDS("tissues/ColonTransverse/DML_results.rds")
# smoking2 <- results$Smoking2
# signif <- smoking2[smoking2$adj.P.Val<0.05,]

#Functions to plot

plot_go <- function(go, direction){
  go <- go[go$ONTOLOGY=="BP",] #Only BP
  sig <- sum(go$FDR<0.05)
  print(sig)
  if(sig>20){
    max<-20
  } else{max<-sig}
  topgo <- topGSA(go, max)
  ggplot(data = topgo, aes(x=DE/N, y = factor(TERM, levels=rev(TERM)), 
                           color = FDR, size = DE)) + 
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    theme_bw() + ylab("") + xlab("Gene Ratio") + 
    ggtitle(paste0("GO ", direction, " BP"))
}


#Enrichment separated by genomic location:
#We first classify as promoters the CpG annotated as "Promoter_Associated" and TSS. If a CpG is associated to two genes, one as promoter and another as gene body, we keep the CpG only as associated to the gene as promoter and exclude the other gene (very few cases) 
test <- annotation %>% separate_rows(UCSC_RefGene_Name, UCSC_RefGene_Group, sep = ';')
promoter <- distinct(test[test$Regulatory_Feature_Group=="Promoter_Associated" | grepl("TSS200|TSS1500", test$UCSC_RefGene_Group), c("Name", "UCSC_RefGene_Name")])
promoter_cpg <- unique(promoter$Name)

#We then classify the remaining cpgs as enhancers:
test <- test[!test$Name %in% promoter_cpg,]
enhancer <- distinct(test[test$Phantom5_Enhancers!="", c("Name", "UCSC_RefGene_Name")])
enhancer_cpg <- unique(enhancer$Name)

#We then assign the other cpg as either gene body or intergenic:
test <- test[!test$Name %in% enhancer_cpg,]
body <- test[grepl("Body|1stExon|5URT|3UTR|ExonBnd", test$UCSC_RefGene_Group), c("Name", "UCSC_RefGene_Name")]
body_cpg <- unique(body$Name)

test <- test[!test$Name %in% body_cpg,]
intergenic_cpg <- test$Name


smoking_hypo <- signif[signif$logFC<0,]
prom_hypo <- rownames(smoking_hypo)[rownames(smoking_hypo) %in% promoter_cpg]

# go_hypo <- gometh(sig.cpg=prom_hypo, all.cpg = rownames(signif), collection="GO", array.type="EPIC") #Output with genes takes a lot
go_hypo <- gometh(sig.cpg=prom_hypo, all.cpg = annotation$IlmnID, collection="GO", array.type="EPIC") #Output with genes takes a lot
# go_hypo <- gometh(sig.cpg=prom_hypo, all.cpg = promoter_cpg, collection="GO", array.type="EPIC") #Output with genes takes a lot
go_hypo <- go_hypo[go_hypo$ONTOLOGY=="BP",]
sum(go_hypo$FDR<0.05)
plot_go(go_hypo, "hypo")
results <- cbind(go_hypo[go_hypo$FDR<0.05,], "region"="promoters", "direction"="hypo")

smoking_hyper <- signif[signif$logFC>0,]
prom_hyper <- rownames(smoking_hyper)[rownames(smoking_hyper) %in% promoter_cpg]

# go_hyper <- gometh(sig.cpg=prom_hyper, all.cpg = rownames(signif), collection="GO", array.type="EPIC") #Output with genes takes a lot
# go_hyper <- gometh(sig.cpg=prom_hyper, all.cpg = promoter_cpg, collection="GO", array.type="EPIC") #Output with genes takes a lot
go_hyper <- gometh(sig.cpg=prom_hyper, all.cpg = annotation$IlmnID, collection="GO", array.type="EPIC") #Output with genes takes a lot
go_hyper <- go_hyper[go_hyper$ONTOLOGY=="BP",]
sum(go_hyper$FDR<0.05)
plot_go(go_hyper, "hyper")
if(sum(go_hyper$FDR<0.05)==0){
  results <- rbind(results, cbind(ONTOLOGY=NA, TERM=NA, N=0, DE=0, P.DE=1, FDR=1, "region"="promoters", "direction"="hyper"))
} else{
  results <- rbind(results, cbind(go_hyper[go_hyper$FDR<0.05,], "region"="promoters", "direction"="hyper"))
}
#If promoter hypo was 0:
# results <- cbind(go_hyper[go_hyper$FDR<0.05,], "region"="promoters", "direction"="hyper")

#Enhancer
enh_hypo <- rownames(smoking_hypo)[rownames(smoking_hypo) %in% enhancer_cpg]
enh_hyper <- rownames(smoking_hyper)[rownames(smoking_hyper) %in% enhancer_cpg]

# go_hypo <- gometh(sig.cpg=enh_hypo, all.cpg = union(enh_hypo, enh_hyper), collection="GO", array.type="EPIC") #Output with genes takes a lot
# go_hypo <- gometh(sig.cpg=enh_hypo, all.cpg = rownames(signif), collection="GO", array.type="EPIC") #Output with genes takes a lot
go_hypo <- gometh(sig.cpg=enh_hypo, all.cpg = annotation$IlmnID, collection="GO", array.type="EPIC") #Output with genes takes a lot
# go_hypo <- gometh(sig.cpg=enh_hypo, all.cpg = enhancer_cpg, collection="GO", array.type="EPIC") #Output with genes takes a lot
# go_hypo <- gometh(sig.cpg=enh_hypo, collection="GO", array.type="EPIC") #Output with genes takes a lot
go_hypo <- go_hypo[go_hypo$ONTOLOGY=="BP",]
sum(go_hypo$FDR<0.05)
plot_go(go_hypo, "hypo")
if(sum(go_hypo$FDR<0.05)==0){
  results <- rbind(results, cbind(ONTOLOGY=NA, TERM=NA, N=0, DE=0, P.DE=1, FDR=1, "region"="enhancers", "direction"="hypo"))
} else{
  results <- rbind(results, cbind(go_hypo[go_hypo$FDR<0.05,], "region"="enhancers", "direction"="hypo"))
}
# #Orsum needs the file to be sorted
# go_hypo <- go_hypo[go_hypo$FDR<0.05,]
# go_hypo <- go_hypo[go_hypo$ONTOLOGY=="BP",] #Orsum only showd BP anyways
# to_save <- rownames(go_hypo)
# write.table(to_save, "output/enrichments_hypo_enhancers.txt", quote = F, col.names = F, row.names = F)


# go_hyper <- gometh(sig.cpg=enh_hyper, all.cpg = rownames(signif), collection="GO", array.type="EPIC") #Output with genes takes a lot
go_hyper <- gometh(sig.cpg=enh_hyper, all.cpg = annotation$IlmnID, collection="GO", array.type="EPIC") #Output with genes takes a lot
# go_hyper <- gometh(sig.cpg=enh_hyper, all.cpg = enhancer_cpg, collection="GO", array.type="EPIC") #Output with genes takes a lot
go_hyper <- go_hyper[go_hyper$ONTOLOGY=="BP",]
sum(go_hyper$FDR<0.05)
plot_go(go_hyper, "hyper")
if(sum(go_hyper$FDR<0.05)==0){
  results <- rbind(results, cbind(ONTOLOGY=NA, TERM=NA, N=0, DE=0, P.DE=1, FDR=1, "region"="enhancers", "direction"="hyper"))
} else{
  results <- rbind(results, cbind(go_hyper[go_hyper$FDR<0.05,], "region"="enhancers", "direction"="hyper"))
}

# go_hyper <- go_hyper[go_hyper$FDR<0.05,]
# go_hyper <- go_hyper[go_hyper$ONTOLOGY=="BP",] #Orsum only showd BP anyways
# to_save <- rownames(go_hyper)
# write.table(to_save, "output/enrichments_hyper_enhancers.txt", quote = F, col.names = F, row.names = F)

go_hypo$direction <- "Hypomethylated"
go_hyper$direction <- "Hypermethylated"
to_plot_f <- rbind(go_hypo, go_hyper)
to_plot_f <- to_plot_f[order(to_plot_f$FDR),]

saveRDS(to_plot_f, "../figures/data/figure_6_enrichment_enhancers.rds")



#Gene body
body_hypo <- rownames(smoking_hypo)[rownames(smoking_hypo) %in% body_cpg]
body_hyper <- rownames(smoking_hyper)[rownames(smoking_hyper) %in% body_cpg]

# go_hypo <- gometh(sig.cpg=body_hypo, all.cpg = rownames(signif), collection="GO", array.type="EPIC") #Output with genes takes a lot
go_hypo <- gometh(sig.cpg=body_hypo, all.cpg = annotation$IlmnID, collection="GO", array.type="EPIC") #Output with genes takes a lot
# go_hypo <- gometh(sig.cpg=body_hypo, all.cpg = body_cpg, collection="GO", array.type="EPIC") #Output with genes takes a lot
go_hypo <- go_hypo[go_hypo$ONTOLOGY=="BP",]
plot_go(go_hypo, "hypo")
if(sum(go_hypo$FDR<0.05)==0){
  results <- rbind(results, cbind(ONTOLOGY=NA, TERM=NA, N=0, DE=0, P.DE=1, FDR=1, "region"="gene_body", "direction"="hypo"))
} else{
  results <- rbind(results, cbind(go_hypo[go_hypo$FDR<0.05,], "region"="gene_body", "direction"="hypo"))
}

# go_hyper <- gometh(sig.cpg=body_hyper, all.cpg = rownames(signif), collection="GO", array.type="EPIC") #Output with genes takes a lot
go_hyper <- gometh(sig.cpg=body_hyper, all.cpg = annotation$IlmnID, collection="GO", array.type="EPIC") #Output with genes takes a lot
# go_hyper <- gometh(sig.cpg=body_hyper, all.cpg = body_cpg, collection="GO", array.type="EPIC") #Output with genes takes a lot
go_hyper <- go_hyper[go_hyper$ONTOLOGY=="BP",]
sum(go_hyper$FDR<0.05)
plot_go(go_hyper, "hyper")
if(sum(go_hyper$FDR<0.05)==0){
  results <- rbind(results, cbind(ONTOLOGY=NA, TERM=NA, N=0, DE=0, P.DE=1, FDR=1, "region"="gene_body", "direction"="hyper"))
} else{
  results <- rbind(results, cbind(go_hyper[go_hyper$FDR<0.05,], "region"="gene_body", "direction"="hyper"))
}

#Share table of results
results <- results[!is.na(results$ONTOLOGY),]
results$FDR <- as.numeric(results$FDR)
results$direction <- factor(results$direction, levels=c("hypo", "hyper"))
results <- results[order(results$region, results$direction, results$FDR), ]
results <- results[,-1]
results$ID <- rownames(results)
rownames(results) <- NULL
colnames(results) <- c("Description", "Number of genes in the category", "Number of DEGs in the category",
                       "p-value", "FDR", "region", "direction", "ID")
results <- results[,c("ID", "Description", "Number of genes in the category", "Number of DEGs in the category",
                       "p-value", "FDR", "region", "direction")]

library("xlsx")
write.xlsx(results, "output/Supplementary_table_14.xlsx", 
           col.names = TRUE, append = FALSE)


# Intergenic has obviously 0 terms
# inter_hyper <- rownames(smoking_hyper)[rownames(smoking_hyper) %in% intergenic_cpg]
# inter_hypo <- rownames(smoking_hypo)[rownames(smoking_hypo) %in% intergenic_cpg]
# 
# go_hypo <- gometh(sig.cpg=inter_hypo, all.cpg = annotation$IlmnID, collection="GO", array.type="EPIC") #Output with genes takes a lot
# go_hypo <- go_hypo[go_hypo$ONTOLOGY=="BP",]
# sum(go_hypo$FDR<0.05)
# plot_go(go_hypo, "hypo")
# 
# 
# go_hyper <- gometh(sig.cpg=inter_hyper, all.cpg = annotation$IlmnID, collection="GO", array.type="EPIC") #Output with genes takes a lot
# go_hyper <- go_hyper[go_hyper$ONTOLOGY=="BP",]
# sum(go_hyper$FDR<0.05)
# plot_go(go_hyper, "hyper")


