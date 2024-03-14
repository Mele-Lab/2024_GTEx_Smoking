# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Generate publication Figure 4 (DNA methylation)


#Set path 
# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio

#Load libraries 
library(ggplot2)
library(ggpubr)
library(RColorBrewer) 
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))


#Load data
figure_data <- readRDS(file = "data/methylation_results.rds")
for(i in 1:length(figure_data)) assign(names(figure_data)[i], figure_data[[i]])

source("safe_colourblind_pallete.R")


row_ha_left <- HeatmapAnnotation("Samples" = anno_barplot(t(sapply(rownames(to_plot), function(tissue) table(metadata[[tissue]]$Smoking))),
                                                          gp = gpar(fill = safe_colorblind_palette[c(1,3,2)]),
                                                          border=F),
                                 gap = unit(0.25,"cm"),
                                 show_legend = T, 
                                 show_annotation_name = T,
                                 annotation_name_rot = 90,
                                 annotation_name_gp = gpar(fontsize = 10),
                                 which = "row")


# ht_DSGs <- Heatmap(apply(to_plot[2], 2,function(x) x/max(x,na.rm=T)),
to_plot <- as.matrix(to_plot)
to_plot[,2] <- as.numeric(to_plot[,2])
to_plot[,3] <- as.numeric(to_plot[,3])
ht_DSGs <- Heatmap(apply(to_plot[,2:3], 2, as.numeric),
                   col = colorRamp2(c(0,1, 50000), brewer.pal(9,"BuPu")[c(1,2,7)]),
                   na_col = "white",
                   cluster_rows = F,
                   cluster_columns = F,
                   name = "Number of DMPs",
                   row_names_side = "left",
                   row_labels = to_plot[,1],
                   column_names_side = "bottom",
                   left_annotation = row_ha_left,
                   row_names_gp = gpar(fontsize = 10),
                   column_names_gp = gpar(fontsize = 10),
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(prettyNum(to_plot[,2:3][i, j], big.mark = ","), x, y, gp = gpar(fontsize = 10))}#,
                   # heatmap_legend_param = list(direction = "horizontal")
)

pdf(file = "figures/figure_4/heatmap.pdf", w = 3.9, h = 2.8)
draw(ht_DSGs)#, heatmap_legend_side = "bottom")
dev.off()


#Gene location information:
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

#Main
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
                   labels=c("Enhancer", "Promoter", "Gene body", "Intergenic")) + xlim(-1, 2)


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
                   labels=c("Enhancer", "Promoter", "Gene body", "Intergenic")) +
  scale_x_continuous(breaks=c(0, 20000, 40000)) #Only for lung


pdf(file = paste0("figures/figure_4/genomic_location_", tissue,".pdf"), w = 6, h = 2)
ggarrange(g1, g2, ncol=2, common.legend = TRUE, legend="right", widths = c(1, 0.6))
dev.off()



#ChromHMM plot
tissue <- "Lung"
# tissue <- "ColonTransverse"
hypo <- readRDS(paste0('../analysis/output/enrichment_chromhmm_hypo_', tissue, '.rds'))
hyper <- readRDS(paste0('../analysis/output/enrichment_chromhmm_hyper_', tissue, '.rds'))

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

g <- ggplot(hyper_hypo, aes(x=log(oddsRatio), y=region, colour=type, alpha=sig)) + 
  geom_errorbar(aes(xmin=log(CI_down), xmax=log(CI_up)), width=.3) + 
  geom_vline(xintercept = 0) +
  geom_point() + ylab('') + theme_bw() +
  scale_colour_manual(values=c("#CC6677", "#88CCEE")) +
  xlab("Log(Odds ratio)") +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour="black", size=9),
        axis.text.y = element_text(colour="black", size=10),
        axis.title.x = element_text(size=10)) +
  scale_alpha_discrete(range = c(0.3, 1))

g2 <- ggplot(hyper_hypo) + geom_col(aes(sample_size, region, fill=type), width = 0.7) +
  theme_bw() + xlab("Number of DMPs") + ylab("") +
  scale_fill_manual(values=c("#CC6677", "#88CCEE")) +
  theme(legend.position = "none",
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour="black", size=9),
        axis.title.x = element_text(size=10)) +
scale_x_continuous(breaks=c(0, 10000, 20000))


pdf(file = paste0("figures/figure_5/chromHMM_", tissue,".pdf"), w = 5.5, h = 3)
ggarrange(g, g2, ncol=2, common.legend = TRUE, legend="right", widths = c(1, 0.4))
dev.off()



#OR plot is in figure_S4.R

#Correlations in main
to_plot <- matrix(c(14.6, 11.9, "DMPs\nassociated\nwith\nDEGs", "DEGs\nassociated\nwith\nDMPs"), nrow=2)
to_plot <- as.data.frame(to_plot)
colnames(to_plot) <- c("value", "variable")
to_plot$value <- as.numeric(to_plot$value)
to_plot$variable <- as.factor(to_plot$variable)
g <- ggplot(to_plot) + geom_col(aes(variable, value, fill=variable), width = 0.9) + xlab("") + 
  ylab("Percentage of counts") +
  scale_fill_manual(values=c("#88CCEE", "#CC6677")) + theme_classic() +
  theme(legend.position="none",
        axis.title.y = element_text(margin = margin(r = 2), size = 11),
        axis.text.y = element_text(size = 9, colour = "black"), 
        axis.text.x = element_text(size = 9, colour = "black"))
g
pdf(file = paste0("figures/figure_4/percentages_associations.pdf"), w = 2, h = 2.5)
g
dev.off()

#Correlations in supplements
to_plot <- readRDS("data/correlation_percentages.rds")

g <- ggplot(to_plot) + geom_col(aes(variable, value, fill=Correlation), width = 0.9) + xlab("") + 
  ylab("Percentage of correlations") +
  scale_fill_manual(values=c("#88CCEE", "#CC6677")) + theme_classic() +
  theme(axis.title.y = element_text(margin = margin(r = 2), size = 11),
        axis.text.y = element_text(size = 9, colour = "black"), 
        axis.text.x = element_text(size = 9, colour = "black"))
#This plot to supplements?
pdf(file = paste0("figures/figure_4/percentages.pdf"), w = 2.5, h = 2.2)
g
dev.off()



#Correlation enrichments:
hypo <- readRDS(paste0('../analysis/output/enrichment_correlated_hypo_lung.rds'))
hyper <- readRDS(paste0('../analysis/output/enrichment_correlated_hyper_lung.rds'))

hypo_d <- read_data(c("promoter", "enhancer", "gene_body"), hypo)
hyper_d <- read_data(c("promoter", "enhancer", "gene_body"), hyper)
hyper_hypo <- rbind(hypo_d, hyper_d)
hyper_hypo$sig[hyper_hypo$sig =="Sig"] <- "FDR < 0.05"
hyper_hypo$sig[hyper_hypo$sig =="not Sig"] <- "FDR >= 0.05"
hyper_hypo$sig <- factor(hyper_hypo$sig, levels=c("FDR >= 0.05", "FDR < 0.05"))
hyper_hypo$region <- factor(hyper_hypo$region, levels=rev(c("promoter", "enhancer", "gene_body")))
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
  scale_y_discrete(breaks=c("enhancer", "promoter", "gene_body"),
                   labels=c("Enhancer", "Promoter", "Gene body")) + xlim(-1, 2) #xlim(-5, 5)

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

pdf(file = paste0("figures/figure_4/genomic_location_correlations.pdf"), w = 6, h = 2)
ggarrange(g1, g2, ncol=2, common.legend = TRUE, legend="right", widths = c(1, 0.6))
dev.off()



#AHRR example correlation
figure_data <- readRDS(file = "data/AHRR_correlation_example.rds")
for(i in 1:length(figure_data)) assign(names(figure_data)[i], figure_data[[i]])

g <- ggplot(scatter_map4k1, aes(Methylation, Expression)) + geom_point() + theme_bw() + 
  ggtitle("AHRR - cg07943658") + geom_smooth(method="lm", colour="#CC6677") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("Methylation residuals") + ylab("Expression residuals") +
  annotate("text", y=max(scatter_map4k1$Expression), x=min(scatter_map4k1$Methylation), 
           label=paste0("r=", correlation, ", p-value=", p_value), hjust=0.03, vjust=-0.1, size=2.95)

pdf(file = "figures/figure_4/correlation_AHRR.pdf", w = 2.5, h = 2.2)
g
dev.off()

#Last plot:
to_plot <- matrix(c(89.81481, 10.18519, "Anticorrelated", "Correlated"), nrow=2)
to_plot <- as.data.frame(to_plot)
colnames(to_plot) <- c("value", "variable")
to_plot$value <- as.numeric(to_plot$value)
to_plot$variable <- factor(to_plot$variable, levels=c("Correlated", "Anticorrelated"))
to_plot$dummy <-"dummy"
size <- 16
# g <- ggplot(to_plot) + geom_col(aes(dummy, value, fill=variable), width = 0.9) + xlab("") + 
#   ylab("Percentage of\ncorrelations") +
#   scale_fill_manual(values=c("#88CCEE", "#CC6677")) + theme_classic() +
#   theme(legend.title = element_blank(),
#         axis.title.y = element_text(margin = margin(r = 2), size = size),
#         axis.text.y = element_text(size = size-1, colour = "black"), 
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         legend.text = element_text(size=size),
#         legend.position = "bottom")
g <- ggplot(to_plot) + geom_col(aes(dummy, value, fill=variable), width = 0.9) + xlab("") + 
  ylab("Percentage of correlations") +
  scale_fill_manual(values=c("#88CCEE", "#CC6677")) + theme_classic() +
  coord_flip() +
  theme(legend.title = element_blank(),
        axis.title.x = element_text(margin = margin(r = 2), size = size),
        axis.text.x = element_text(size = size-1, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.text = element_text(size=size),
        legend.position = "bottom")
g
# pdf(file = "figures/figure_5/hypomethylated_enhancers.pdf", w = 1.8, h = 3)
pdf(file = "figures/figure_5/hypomethylated_enhancers.pdf", w = 3.2, h = 1.8)
g
dev.off()
