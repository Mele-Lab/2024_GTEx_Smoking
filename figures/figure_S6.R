#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Generate publication Figure S6


#Set path 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Loading libraries
suppressPackageStartupMessages(library(ComplexHeatmap))
library(RColorBrewer)
suppressPackageStartupMessages(library(circlize))
library(ggplot2)
library(ggpubr)


#General methylation
library(readr)
chromhmm_cpgs <- read_csv("../analysis/data/public/lung_chromhmm.csv") #Preprocessed chromHMM data, we preprocess it in a preivous step of the pipeline

beta <- readRDS("../analysis/tissues/Lung/methylation_data.rds") #We need to run the first part of the pipeline to get this object, it is not in the github
probes <- rownames(beta)
beta <- sapply(beta, as.numeric)
rownames(beta) <- probes

to_plot <- rowMeans(beta)
to_plot <- as.data.frame(to_plot)
to_plot$cpg <- rownames(to_plot)
to_plot <- merge(to_plot, chromhmm_cpgs, by.x="cpg", by.y="name_ann")

library(ggplot2)

library(dplyr)
to_plot %>%  group_by(region_chromhmm) %>% summarise(n=n()) ->Summary.data

to_plot$region_chromhmm <- factor(to_plot$region_chromhmm, levels=c("Active enhancer",
                                                                    "Genic enhancer",
                                                                    "Weak enhancer",
                                                                    "Active TSS", "Flanking TSS",
                                                                    "Strong transcription",
                                                                    "Weak transcription",
                                                                    "Heterochromatin",
                                                                    "Quiescent",
                                                                    "Bivalent enhancer",
                                                                    "Bivalent TSS",
                                                                    "Repressed polycomb",
                                                                    "Weak repressed polycomb",
                                                                    "ZNF genes & repeats"))
saveRDS(to_plot, 'data/general_methylation.rds')
to_plot <- readRDS('data/general_methylation.rds')

g <- ggplot(to_plot) + geom_violin(aes(region_chromhmm, to_plot), scale="width", fill="gray90") + xlab("") + ylab("Beta values") +
  geom_text(data=Summary.data ,aes(x = region_chromhmm, y = 1.1, label=n), fontface =1, size = 3.2) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,
                                                size = 14, colour="black"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title = element_text(size = 15),
                     axis.text.y = element_text(size = 12.5, colour = "black"))

ggsave("figures/figure_s6/methylation.pdf", g, device = "pdf", width = 6.8, height = 4.5)



#Replication in blood:
read_data <- function(variables, data){ #Function to prepare data to plot and compute adjusted p value
  
  odds_ratio <- lapply(variables, function(trait) data[[trait]][['f']]$estimate)
  adj.P.Val <- p.adjust(sapply(variables, function(trait) data[[trait]][['f']]$p.value), method = "BH")
  CI_down <- lapply(variables, function(trait) data[[trait]][['f']]$conf.int[1])
  CI_up <- lapply(variables, function(trait) data[[trait]][['f']]$conf.int[2])
  sample_size <- lapply(variables, function(trait) data[[trait]][['m']])
  
  
  names(odds_ratio) <- variables
  names(adj.P.Val) <- variables
  names(CI_down) <- variables
  names(CI_up) <- variables
  names(sample_size) <- variables
  
  odds_ratio_df <- as.data.frame(unlist(odds_ratio))
  odds_ratio_df$label <- variables
  odds_ratio_df$type <- deparse(substitute(data)) #Either hypo or hyper
  colnames(odds_ratio_df) <- c('oddsRatio', 'region','type')
  
  adj.P.Val_df <- as.data.frame(unlist(adj.P.Val))
  adj.P.Val_df$label <- variables
  adj.P.Val_df$type <- deparse(substitute(data))
  colnames(adj.P.Val_df) <- c('adjPvalue','region','type')
  
  CI_down_df <- as.data.frame(unlist(CI_down))
  CI_down_df$label <- variables
  CI_down_df$type <- deparse(substitute(data))
  colnames(CI_down_df) <- c('CI_down','region','type')
  
  CI_up_df <- as.data.frame(unlist(CI_up))
  CI_up_df$label <- variables
  CI_up_df$type <- deparse(substitute(data))
  colnames(CI_up_df) <- c('CI_up','region','type')
  
  sample_size_df <- as.data.frame(unlist(sample_size))
  sample_size_df$label <- variables
  sample_size_df$type <- deparse(substitute(data))
  colnames(sample_size_df) <- c('sample_size','region','type')
  
  all <- Reduce(function(x, y) merge(x, y, all=TRUE), list(odds_ratio_df, adj.P.Val_df, CI_down_df, CI_up_df, sample_size_df))
  head(all)
  all$sig <- 'not Sig'
  all$sig[all$adjPvalue<0.05] <- 'Sig'
  all <- all[,c("region","oddsRatio","adjPvalue","CI_down","CI_up","sig","type", "sample_size")]
  return(all)
}

hypo <- readRDS(paste0('data/enrichment_chromhmm_hypo_blood.rds'))
hyper <- readRDS(paste0('data/enrichment_chromhmm_hyper_blood.rds'))

hypo_d <- read_data(names(hypo), hypo)
hyper_d <- read_data(names(hyper), hyper)
hyper_hypo <- rbind(hypo_d, hyper_d)
#Change order to plot
hyper_hypo$region <- factor(hyper_hypo$region, levels=rev(c("Active TSS", "Flanking TSS", "Bivalent TSS",
                                                            "ZNF genes & repeats", "Heterochromatin", "Quiescent",
                                                            "Weak repressed polycomb", "Repressed polycomb",
                                                            "Weak transcription", "Strong transcription",
                                                            "Weak enhancer", "Active enhancer", "Genic enhancer", "Bivalent enhancer")))
hyper_hypo$type[hyper_hypo$type =="hypo"] <- "Hypomethylation"
hyper_hypo$type[hyper_hypo$type =="hyper"] <- "Hypermethylation"

hyper_hypo$sig[hyper_hypo$sig =="Sig"] <- "FDR < 0.05"
hyper_hypo$sig[hyper_hypo$sig =="not Sig"] <- "FDR >= 0.05"
hyper_hypo$sig <- factor(hyper_hypo$sig, levels=c("FDR >= 0.05", "FDR < 0.05"))

size <- 12
g <- ggplot(hyper_hypo, aes(x=log(oddsRatio), y=region, colour=type, alpha=sig)) + 
  geom_errorbar(aes(xmin=log(CI_down), xmax=log(CI_up)), width=.3) + 
  geom_vline(xintercept = 0) +
  geom_point(size=3) + ylab('') + theme_bw() +
  scale_colour_manual(values=c("#CC6677", "#88CCEE")) +
  xlab("Log(Odds ratio)") +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour="black", size=size),
        axis.text.y = element_text(colour="black", size=size+1),
        axis.title.x = element_text(size=size+1),
        legend.text = element_text(size=size+1),
        legend.key.size = unit(0.9,"cm")) +
  scale_alpha_discrete(range = c(0.3, 1)) +
  guides(color = guide_legend(override.aes = list(size = 3)),
         alpha = guide_legend(override.aes = list(size = 3))) 

g2 <- ggplot(hyper_hypo) + geom_col(aes(sample_size, region, fill=type), width = 0.7) +
  theme_bw() + xlab("Number of DMPs") + ylab("") +
  scale_fill_manual(values=c("#CC6677", "#88CCEE")) +
  theme(legend.position = "none",
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour="black", size=size),
        axis.title.x = element_text(size=size+1)) +
  scale_x_continuous(breaks=c(0, 1000, 2000))

pdf(file = paste0("figures/figure_s6/chromHMM_blood.pdf"), w = 7.2, h = 4.2)
ggarrange(g, g2, ncol=2, common.legend = TRUE, legend="right", widths = c(1, 0.35))
dev.off()

#TFBS:

library(ComplexHeatmap)
library(RColorBrewer)


table <- readRDS("data/TFBS_conservative.rds")
table <- table[table$TF!="Epitope",] #This should be removed from previous code
table$region <- factor(table$region, levels=c("Active TSS", "Flanking TSS", "Bivalent TSS",
                                              "ZNF genes & repeats", "Heterochromatin", "Quiescent",
                                              "Weak repressed polycomb", "Repressed polycomb",
                                              "Weak transcription", "Strong transcription",
                                              "Weak enhancer", "Active enhancer", "Genic enhancer", "Bivalent enhancer"))
hypo <- table[table$direction=="hypo" & table$adj_p_val<0.05,]
hyper <- table[table$direction=="hyper" & table$adj_p_val<0.05,]

heatmap_data <- cbind(table(hyper$region), table(hypo$region))
colnames(heatmap_data) <- c("Hyper", "Hypo")

size <- 13
ht_DEGs <- Heatmap(heatmap_data,
                   col = brewer.pal(9,"BuPu")[1:7],
                   cluster_rows = F,
                   cluster_columns = F,
                   name = "Number of TFBS",
                   row_names_side = "left",
                   column_names_side = "bottom",
                   row_names_gp = gpar(fontsize = size),
                   column_names_rot = 0,
                   column_names_gp = gpar(fontsize = size),
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(prettyNum(heatmap_data[i, j], big.mark = ","), x, y, gp = gpar(fontsize = size))},
                   heatmap_legend_param = list(direction = "vertical", title_position = "topcenter",
                                               labels_gp = gpar(fontsize = size+1),
                                               title_gp = gpar(fontsize = size+1, fontface = "bold"))
)
pdf("figures/figure_s6/TFBS_heatmap.pdf", width = 4, height = 4)
ht_DEGs
dev.off()



hypo <- table[table$direction=="hypo",]
hyper <- table[table$direction=="hyper",]
regions <- c("Active TSS", "Flanking TSS", "Bivalent TSS",
             "ZNF genes & repeats", "Heterochromatin", "Quiescent",
             "Weak repressed polycomb", "Repressed polycomb",
             "Weak transcription", "Strong transcription",
             "Weak enhancer", "Active enhancer", "Genic enhancer", "Bivalent enhancer")

#Keep TFBS that are enriched in 4 or more regions:
hyper_subset <- hyper[hyper$adj_p_val<0.05 & hyper$odds_ratio>1,]
hyper_subset <- table(hyper_subset$TF)>=4
tbfs_hyper <- rownames(hyper_subset)[hyper_subset]

hypo_subset <- hypo[hypo$adj_p_val<0.05 & hypo$odds_ratio>1,]
hypo_subset <- table(hypo_subset$TF)>=9
tbfs_hypo <- rownames(hypo_subset)[hypo_subset]


#Plot:
hypo_plot <- matrix(nrow = 14, ncol=length(tbfs_hypo), dimnames = list(regions, tbfs_hypo))
hyper_plot <- matrix(nrow = 14, ncol=length(tbfs_hyper), dimnames = list(regions, tbfs_hyper))

for(region in regions){
  for(tf in tbfs_hyper){
    entry <- hyper$TF==tf & hyper$region==region
    if(hyper$adj_p_val[entry]<0.05 & hyper$odds_ratio[entry]>1){
      hyper_plot[region, tf] <- hyper$odds_ratio[hyper$TF==tf & hyper$region==region]
      # hyper_plot[region, tf] <- hyper$adj_p_val[hyper$TF==tf & hyper$region==region]
    }
  }
  for(tf in tbfs_hypo){
    entry <- hypo$TF==tf & hypo$region==region
    if(hypo$adj_p_val[entry]<0.05 & hypo$odds_ratio[entry]>1){
      hypo_plot[region, tf] <- hypo$odds_ratio[hypo$TF==tf & hypo$region==region]
      # hypo_plot[region, tf] <- hypo$adj_p_val[hypo$TF==tf & hypo$region==region]
      
    }
  }
}

hyper_plot <- apply(hyper_plot, c(1, 2), as.numeric)
hypo_plot <- apply(hypo_plot, c(1, 2), as.numeric)

pdf("figures/figure_5/TFBS_hyper_Polycomb.pdf", width = 6.3, height = 4.5) #if png -> units="in", res=200
Heatmap(log(hyper_plot),
        col = brewer.pal(9,"BuPu")[3:7],
        na_col = "white",
        cluster_rows = F,
        cluster_columns = F,
        name = "Log(Odds ratio)",
        row_names_gp = gpar(fontsize = size+1.5),
        column_names_gp = gpar(fontsize = size),
        row_names_side = "left",
        heatmap_legend_param = list(title_position = "topcenter",
                                    labels_gp = gpar(fontsize = size+1),
                                    title_gp = gpar(fontsize = size+2, fontface = "bold")))
dev.off()

size <- 13
pdf("figures/figure_s6/TFBS_hypo_Polycomb.pdf", width = 8, height = 4.5) #if png -> units="in", res=200
Heatmap(log(hypo_plot),
        col = brewer.pal(9,"BuPu")[3:7],
        na_col = "white",
        cluster_rows = F,
        cluster_columns = F,
        name = "Log(Odds ratio)",
        row_names_gp = gpar(fontsize = size),
        column_names_gp = gpar(fontsize = size),
        row_names_side = "left",
        heatmap_legend_param = list(title_position = "topcenter",
                            labels_gp = gpar(fontsize = size),
                            title_gp = gpar(fontsize = size, fontface = "bold")))
dev.off()



#Check functional enrichments, iterate over hypo and hyper and over regions
library(clusterProfiler)
library(org.Hs.eg.db)

results_go <- data.frame("ID"=1, "Description"=1, "GeneRatio"=1, "BgRatio"=1, "pvalue"=1, "p.adjust"=1, "region"=1, "direction"=1)

for(region in unique(table$region)){
  bg <- unique(table$TF[table$region==region]) #221 TFs
  
  for(direction in c("hypo", "hyper")){
    gl <- table$TF[table$region==region & table$adj_p_val<0.05 & table$direction==direction & table$odds_ratio>1]
    if(rlang::is_empty(gl)){next}
    ora.go <- enrichGO(gene     = gl, universe =    bg,
                       keyType = "SYMBOL", 
                       OrgDb        = org.Hs.eg.db,
                       ont          = "BP", 
                       minGSSize    = 10, maxGSSize    = 500,
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                       readable = F)
    if(sum(ora.go@result$p.adjust<0.05)>0){
      print(paste("There are enrichments:", region, direction))
      print(dotplot(ora.go, showCategory = 15, title=paste0(region, "_", direction)))

      results_go <- rbind(results_go, cbind(ora.go[ora.go$p.adjust<0.05,1:6], region, direction))
    }
  }
  
}
results_go <- results_go[-1,]






#Age-smoking-DMPs
table_age <- readRDS("data/TFBS_Age_conservative.rds")
table_age <- table_age[table_age$TF!="Epitope",] #This should be removed from previous code

table_age$region <- factor(table_age$region, levels=c("Active TSS", "Flanking TSS", "Bivalent TSS",
                                                      "ZNF genes & repeats", "Heterochromatin", "Quiescent",
                                                      "Weak repressed polycomb", "Repressed polycomb",
                                                      "Weak transcription", "Strong transcription",
                                                      "Weak enhancer", "Active enhancer", "Genic enhancer", "Bivalent enhancer"))

hypo <- table_age[table_age$direction=="hypo" & table_age$adj_p_val<0.05,]
hyper <- table_age[table_age$direction=="hyper" & table_age$adj_p_val<0.05,]

heatmap_data <- cbind(table(hyper$region), table(hypo$region))
colnames(heatmap_data) <- c("Hyper", "Hypo")

ht_DEGs <- Heatmap(heatmap_data,
                   col = brewer.pal(9,"BuPu")[1:7],
                   cluster_rows = F,
                   cluster_columns = F,
                   name = "Number of TFBS",
                   row_names_side = "left",
                   column_names_side = "bottom",
                   row_names_gp = gpar(fontsize = 12),
                   column_names_rot = 0,
                   column_names_gp = gpar(fontsize = 13),
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(prettyNum(heatmap_data[i, j], big.mark = ","), x, y, gp = gpar(fontsize = 12))},
                   heatmap_legend_param = list(direction = "vertical", title_position = "topcenter",
                                               labels_gp = gpar(fontsize = 13),
                                               title_gp = gpar(fontsize = 13, fontface = "bold"))
)
pdf("figures/figure_s7/TFBS_heatmap_age.pdf", width = 3.8, height = 3.3)
ht_DEGs
dev.off()


hypo <- table_age[table_age$direction=="hypo",]
hyper <- table_age[table_age$direction=="hyper",]

hyper_subset <- hyper[hyper$adj_p_val<0.05 & hyper$odds_ratio>1,]
hyper_subset <- table(hyper_subset$TF)>=8
tbfs_hyper <- rownames(hyper_subset)[hyper_subset]

hypo_subset <- hypo[hypo$adj_p_val<0.05 & hypo$odds_ratio>1,]
hypo_subset <- table(hypo_subset$TF)>=7
tbfs_hypo <- rownames(hypo_subset)[hypo_subset]


#Plot:
hypo_plot <- matrix(nrow = 14, ncol=length(tbfs_hypo), dimnames = list(regions, tbfs_hypo))
hyper_plot <- matrix(nrow = 14, ncol=length(tbfs_hyper), dimnames = list(regions, tbfs_hyper))

for(region in regions){
  for(tf in tbfs_hyper){
    entry <- hyper$TF==tf & hyper$region==region
    if(hyper$adj_p_val[entry]<0.05 & hyper$odds_ratio[entry]>1){
      hyper_plot[region, tf] <- hyper$odds_ratio[hyper$TF==tf & hyper$region==region]
    }
  }
  for(tf in tbfs_hypo){
    entry <- hypo$TF==tf & hypo$region==region
    if(hypo$adj_p_val[entry]<0.05 & hypo$odds_ratio[entry]>1){
      hypo_plot[region, tf] <- hypo$odds_ratio[hypo$TF==tf & hypo$region==region]
    }
  }
}

hyper_plot <- apply(hyper_plot, c(1, 2), as.numeric)
hypo_plot <- apply(hypo_plot, c(1, 2), as.numeric)

pdf("figures/figure_5/Age_TFBS_hyper_Polycomb_final.pdf", width = 5, height = 4) #if png -> units="in", res=200
Heatmap(log(hyper_plot),
        col = brewer.pal(9,"BuPu")[3:7],
        # col = colorRamp2(c(1, 5), c("gray", "#88CCEE")),
        # col = c("Enriched"="#88CCEE", "Depleted"="#CC6677"),
        na_col = "white",
        cluster_rows = F,
        cluster_columns = F,
        name = "Log(Odds ratio)",
        # name = "Odds ratio",
        row_names_side = "left", 
        heatmap_legend_param = list(direction = "vertical", title_position = "topcenter"))
dev.off()

# pdf("figures/figure_s7/Age_TFBS_hypo_Polycomb_final.pdf", width = 6.5, height = 3.7) #if png -> units="in", res=200
pdf("figures/figure_s7/Age_TFBS_hypo_Polycomb_final.pdf", width = 8.5, height = 3.7) #if png -> units="in", res=200
Heatmap(log(hypo_plot),
        col = brewer.pal(9,"BuPu")[3:7],
        # col = colorRamp2(c(1, 5), c("gray", "#88CCEE")),
        # col = c("Enriched"="#88CCEE", "Depleted"="#CC6677"),
        na_col = "white",
        cluster_rows = F,
        cluster_columns = F,
        name = "Log(Odds ratio)",
        # name = "Odds ratio",
        row_names_side = "left",
        heatmap_legend_param = list(direction = "vertical", title_position = "topcenter",
                                    labels_gp = gpar(fontsize = 13),
                                    title_gp = gpar(fontsize = 13, fontface = "bold")))
dev.off()



#Functional enrichment of aging
results_go_age <- data.frame("ID"=1, "Description"=1, "GeneRatio"=1, "BgRatio"=1, "pvalue"=1, "p.adjust"=1, "region"=1, "direction"=1)

for(region in unique(table_age$region)){
  bg <- unique(table_age$TF[table_age$region==region]) #221 TFs
  
  for(direction in c("hypo", "hyper")){
    gl <- table_age$TF[table_age$region==region & table_age$adj_p_val<0.05 & table_age$direction==direction & table_age$odds_ratio>1]
    if(rlang::is_empty(gl)){next}
    ora.go <- enrichGO(gene     = gl, universe =    bg,
                       keyType = "SYMBOL", 
                       OrgDb        = org.Hs.eg.db,
                       ont          = "BP", 
                       minGSSize    = 10, maxGSSize    = 500,
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                       readable = F)
    if(sum(ora.go@result$p.adjust<0.05)>0){
      print(paste("There are enrichments:", region, direction))
      # pdf(paste0("Plots/Age_", region, "_", direction, ".pdf"))
      print(dotplot(ora.go, showCategory = 15, title=paste0(region, "_", direction)))
      # dev.off()
      
      results_go_age <- rbind(results_go_age, cbind(ora.go[ora.go$p.adjust<0.05,1:6], region, direction))
    }
  }
  
}
results_go_age <- results_go_age[-1,]


#Saving parsed tables:
#Get table with the information we want:
function_to_parse <- function(data, variable){
  colnames(data) <- c("Transcrition factor", "p-value", "Odds ratio", "CI lower bound", "CI upper bound", "Region", "Methylation", "Adjusted p-value" )
  data$Direction <- "Enriched"
  data$Direction[data$`Odds ratio`<1] <- "Depleted"
  data$Trait <- variable
  data <- data[,c(10, 1, 2, 6, 8, 7, 3, 4, 5, 9)]
  data <- data[order(data$Direction, data$`Adjusted p-value`, method= "radix", decreasing = c(T, F)),]
  # data <- data[data$`Adjusted p-value`<0.05,]
  return(data)
}

# table <- table[,c(1:7,9)]
# table_age <- table_age[,c(1:7,9)]
test <- rbind(function_to_parse(table, "Smoking"), 
              function_to_parse(table_age, "Smoking-age"))

library("xlsx")
write.xlsx(test, "output/Supplementary_table_8.xlsx",
           col.names = TRUE, row.names = F, append = FALSE)


