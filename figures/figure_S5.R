#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Generate publication Figure S9


#Set path 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Loading libraries
suppressPackageStartupMessages(library(ComplexHeatmap))
library(RColorBrewer)
suppressPackageStartupMessages(library(circlize))
library(ggplot2)
library(ggpubr)


#PEER correlation:
peer_data <- readRDS("data/PEER_correlations.rds")
for(i in 1:length(peer_data)) assign(names(peer_data)[i], peer_data[[i]])

col_fun = colorRamp2(c(-0.5, 0, 0.5), c("#CC6677", "white", "#88CCEE"))

pdf("../figures/figures/figure_s4/PEER_correlations.pdf", width = 1.7, height = 4)
Heatmap(correlations, cluster_rows = F,
        row_names_gp = grid::gpar(fontsize = 9),
        name = "Correlation\nwith\nsmoking",
        col = col_fun)
decorate_heatmap_body("Correlation\nwith\nsmoking", { for (i in 1:length(columns_to_highlight)) {
  grid.lines(c(0, 1), c(top[i],top[i]), gp = gpar(lty = 1, lwd = 1, col = 1))
  grid.lines(c(0, 1), c(btm[i],btm[i]), gp = gpar(lty = 1, lwd = 1, col = 1))
  grid.lines(c(0, 0), c(btm[i],top[i]), gp = gpar(lty = 1, lwd = 1, col = 1))
  grid.lines(c(1, 1), c(btm[i],top[i]), gp = gpar(lty = 1, lwd = 1, col = 1))
}})
dev.off()

pdf("../figures/figures/figure_s4/PEER_correlations_legend.pdf", width = 1, height = 1) #I edited in inkscape to be gray and have a box
lg <- Legend(labels = c("FDR >= 0.05", "FDR < 0.05"), 
             legend_gp = gpar(fill="#88CCEE"), 
             ncol = 1, gap = unit(1, "cm"))
draw(lg)
dev.off()


# #Replication to previous literature
figure_data <- readRDS(file = "data/methylation_literature.rds")
for(i in 1:length(figure_data)) assign(names(figure_data)[i], figure_data[[i]])

# cols <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#999933", "#882255", "#661100", "#888888")
cols <- c("#88CCEE", "#CC6677", "#DDCC77", "#8a51a4")

names(cols) <- c("Lung", "WholeBlood", "AdiposeSubcutaneous", "ColonTransverse")
dmps$cols[9] <- "ColonTransverse"
a <- ggplot(to_plot, aes(log(ors), names, col=cols)) + ylab("") +
  xlab("Log(Odds ratio)") + theme_bw() + scale_colour_manual(values=cols, guide = "none") +
  geom_vline(aes(xintercept = 0), linewidth = .25, linetype = "dashed") +
  # xlim(c(0, NA)) +
  geom_errorbarh(aes(xmax = log(CI_up), xmin = log(CI_down)),
                 position=position_dodge(width=0.9),
                 linewidth = .5, height =.2) +
  geom_point(aes(size = -log10(p_vals))) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 9, colour = "black"),
        axis.title.y = element_text(margin = margin(r = 5), size = 9)) +
  guides(size=guide_legend(title="-log10(p-value)"))

b <- ggplot(dmps, aes(dmps, name, fill=cols)) + geom_col() + theme_bw() +
  scale_fill_manual(values=cols, guide = "none") + xlab("Number of DMPs") +
  ylab("") + geom_text(aes(label=dmps), size=3, hjust = -0.1) +
  scale_x_continuous(breaks=c(0, 12000, 24000), limits = c(0, 28000)) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 9, colour = "black"),
        axis.title.y = element_text(margin = margin(r = 5), size = 9))



pdf(file =  "figures/figure_s5/replication.pdf", w = 6, h = 2.5)
ggpubr::ggarrange(a, b, widths = c(0.7,0.3), common.legend = TRUE, legend = "right")
dev.off()


#CpG island context plot
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

tissue <- "Lung"
hypo <- readRDS(paste0('../analysis/output/final_enrichment_hypo_', tissue,'.rds'))
hyper <- readRDS(paste0('../analysis/output/final_enrichment_hyper_', tissue,'.rds'))

hypo_d <- read_data(names(hypo)[!names(hypo) %in% c("promoter", "enhancer", "gene_body", "intergenic")], hypo)
hyper_d <- read_data(names(hypo)[!names(hypo) %in% c("promoter", "enhancer", "gene_body", "intergenic")], hyper)
hyper_hypo <- rbind(hypo_d, hyper_d)

hyper_hypo$sig[hyper_hypo$sig =="Sig"] <- "FDR < 0.05"
hyper_hypo$sig[hyper_hypo$sig =="not Sig"] <- "FDR >= 0.05"
hyper_hypo$sig <- factor(hyper_hypo$sig, levels=c("FDR >= 0.05", "FDR < 0.05"))
hyper_hypo$region <- factor(hyper_hypo$region, levels=rev(c("island", "shore", "shelf", "open_sea")))
hyper_hypo$type[hyper_hypo$type =="hypo"] <- "Hypomethylation"
hyper_hypo$type[hyper_hypo$type =="hyper"] <- "Hypermethylation"

g1 <- ggplot(hyper_hypo, aes(x=log(oddsRatio), y=region, colour=type, alpha=sig)) + 
  geom_errorbar(aes(xmin=log(CI_down), xmax=log(CI_up)), width=.4) + xlim(c(-2, 1))+
  geom_vline(xintercept = 0) +
  geom_point(size=2.5) + ylab('') + theme_bw() +
  scale_colour_manual(values=c("#CC6677", "#88CCEE")) +
  xlab("Log(Odds ratio)") +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(colour="black", size=11),
        axis.text.y = element_text(colour="black", size=13),
        legend.text = element_text(colour="black", size=12),
        axis.title.x = element_text(size=13),
        legend.spacing.y = unit(-0.05, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth=1)) +
  scale_alpha_discrete(range = c(0.3, 1), drop=FALSE) +
  scale_y_discrete(breaks=c("island", "shore", "shelf", "open_sea"),
                   labels=c("CpG island", "CpG shore", "CpG shelf", "Open sea")) 



#Plot sample sizes:
g2 <- ggplot(hyper_hypo) + geom_col(aes(sample_size, region, fill=type), width = 0.7) +
  theme_bw() + xlab("Number of probes") + ylab("") +
  scale_fill_manual(values=c("#CC6677", "#88CCEE")) +
  scale_y_discrete(breaks=c("island", "shore", "shelf", "open_sea"),
                   labels=c("CpG island", "CpG shore", "CpG shelf", "Open sea")) +
  theme(legend.position = "none",
        axis.text.x = element_text(colour="black", size=11),
        axis.title.x = element_text(size=13),
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth=1)) +
  scale_x_continuous(breaks=c(0, 30000, 60000)) #Only for lung


pdf(file = paste0("figures/figure_s5/genomic_location_", tissue,"_supp.pdf"), w = 6, h = 2)
ggarrange(g1, g2, ncol=2, common.legend = TRUE, legend="right", widths = c(1, 0.6))
dev.off()



#Colon:
tissue <- "ColonTransverse"
hypo <- readRDS(paste0('../analysis/output/final_enrichment_hypo_', tissue,'.rds'))
hyper <- readRDS(paste0('../analysis/output/final_enrichment_hyper_', tissue,'.rds'))

hypo_d <- read_data(c("promoter", "enhancer", "gene_body", "intergenic"), hypo)
hyper_d <- read_data(c("promoter", "enhancer", "gene_body", "intergenic"), hyper)
hyper_hypo <- rbind(hypo_d, hyper_d)

hyper_hypo$sig[hyper_hypo$sig =="Sig"] <- "FDR < 0.05"
hyper_hypo$sig[hyper_hypo$sig =="not Sig"] <- "FDR >= 0.05"
hyper_hypo$sig <- factor(hyper_hypo$sig, levels=c("FDR >= 0.05", "FDR < 0.05"))
hyper_hypo$region <- factor(hyper_hypo$region, levels=rev(c("promoter", "enhancer", "gene_body", "intergenic")))
hyper_hypo$type[hyper_hypo$type =="hypo"] <- "Hypomethylation"
hyper_hypo$type[hyper_hypo$type =="hyper"] <- "Hypermethylation"

g1 <- ggplot(hyper_hypo, aes(x=log(oddsRatio), y=region, colour=type, alpha=sig)) + 
  geom_errorbar(aes(xmin=log(CI_down), xmax=log(CI_up)), width=.4) + 
  geom_vline(xintercept = 0) + 
  geom_point(size=2.5) + ylab('') + theme_bw() +
  scale_colour_manual(values=c("#CC6677", "#88CCEE")) +
  xlab("Log(Odds ratio)") +
  scale_alpha_discrete(range = c(0.4, 1), drop = FALSE) +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(colour="black", size=11),
        axis.text.y = element_text(colour="black", size=13),
        legend.text = element_text(colour="black", size=12),
        axis.title.x = element_text(size=13),
        legend.spacing.y = unit(-0.05, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth=1)) +
  scale_y_discrete(breaks=c("enhancer", "promoter", "gene_body", "intergenic"),
                   labels=c("Enhancer", "Promoter", "Gene body", "Intergenic")) #+ xlim(-1, 2)


#Plot sample sizes:
g2 <- ggplot(hyper_hypo) + geom_col(aes(sample_size, region, fill=type), width = 0.6) +
  theme_bw() + xlab("Number of DMPs") + ylab("") +
  scale_fill_manual(values=c("#CC6677", "#88CCEE")) +
  theme(legend.position = "none",
        axis.text.x = element_text(colour="black", size=11),
        axis.title.x = element_text(size=13),
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth=1)) +
  scale_y_discrete(breaks=c("enhancer", "promoter", "gene_body", "intergenic"),
                   labels=c("Enhancer", "Promoter", "Gene body", "Intergenic"))


pdf(file = paste0("figures/figure_s5/genomic_location_", tissue,".pdf"), w = 6, h = 2)
ggarrange(g1, g2, ncol=2, common.legend = TRUE, legend="right", widths = c(1, 0.6))
dev.off()


