#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Generate publication Figure S8 (Reversibility on Splicing)


#Set path 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Loading libraries
suppressPackageStartupMessages(library(ComplexHeatmap))
library(RColorBrewer)
suppressPackageStartupMessages(library(circlize))
library(ggplot2)

source("safe_colourblind_pallete.R")

#Load data
figure_data <- readRDS(file = "data/figure_S8.rds") #list of expression_barplot, expression_pirate, splicing_barplot and splicing_pirate
for(i in 1:length(figure_data)) assign(names(figure_data)[i], figure_data[[i]])

# Figure S8A ----
#Heatmap
pushViewport(viewport(gp = gpar(fontfamily = "Helvetica")))

pdf("figures/figure_S8/Heatmaps.pdf",  width = 6.8 , height = 3)
# rownames(data) <- gsub("\\(.*","",rownames(data))
Heatmap(t(as.matrix(data)),
        # col = c("white", brewer.pal(9,"Greys")[4:6]),
        # col = colorRamp2(c(0, 1, 62), brewer.pal(9,"BuPu")[c(1, 2:7)]),
        col = brewer.pal(9,"BuPu")[c(1:7)],
        # na_col = "white",
        cluster_rows = F,
        cluster_columns = F,
        name = "Number of DSEs",
        row_names_side = "left",
        column_names_side = "bottom",
        row_names_gp = gpar(fontsize = 12),
        column_names_gp = gpar(fontsize = 12),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(t(data)[i, j], x, y, gp = gpar(fontsize = 11))},
        heatmap_legend_param = list(direction = "vertical"))
dev.off()


# Figure S8B ----
figure_data <- readRDS(file = "data/figure_S8_logFC_reversible.rds") #list of expression_barplot, expression_pirate, splicing_barplot and splicing_pirate
for(i in 1:length(figure_data)) assign(names(figure_data)[i], figure_data[[i]])


#Boxplot showing that ex-smokers are closer to never smokers:
get_box_stats <- function(y, upper_limit) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "N =", length(y), "\n"
    )
  ))
}

p_val <- test$p.value
p_val <- formatC(p_val, format = "e", digits = 2)

to_plot <- as.data.frame(rbind(cbind(logFC[,1], "Never smokers\nvs\nex-smokers"), cbind(logFC[,2], "Ex-smokers\nvs\nsmokers")))
to_plot$V1 <- as.numeric(to_plot$V1)
to_plot$V2 <- factor(to_plot$V2, levels=c("Never smokers\nvs\nex-smokers", "Ex-smokers\nvs\nsmokers"))
g <- ggplot(data = to_plot, aes(V2, V1)) + 
  geom_violin(aes(fill = V2),
              col = "black") +
  geom_boxplot(col = "black",
               outlier.shape = NA,
               notch = T,
               width = 0.25) +
  theme_bw() + ggtitle(paste0("p-value = ", p_val)) +
  xlab("") +
  ylab("logFC") +
  scale_fill_manual(values = c("#88CCEE",  "#CC6677")) +
  stat_summary(fun.data = get_box_stats, fun.args = list(upper_limit = max(to_plot$V1) * 1.15), 
               geom = "text", hjust = 0.5, vjust = 0.9, size = 4) +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13, colour = "black"), 
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "grey"),
        legend.position = "none", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 13))

pdf(file = "figures/figure_s8/logFC_methylation.pdf", w = 3.5, h = 3.5)
g
dev.off()



# Figure S8C-d
figure_data <- readRDS(file = "data/figure_S8_logFC_subset.rds") #list of expression_barplot, expression_pirate, splicing_barplot and splicing_pirate
for(i in 1:length(figure_data)) assign(names(figure_data)[i], figure_data[[i]])


p_val <- test_subset$p.value
p_val <- formatC(p_val, format = "e", digits = 2)

boxplot_subset$Variable[boxplot_subset$Variable=="Never vs Ex"] <- "Never smokers\nvs\nex-smokers"
boxplot_subset$Variable[boxplot_subset$Variable=="Ex vs Smokers"] <- "Ex-smokers\nvs\nsmokers"

boxplot_subset$Variable <- factor(boxplot_subset$Variable, levels=c("Never smokers\nvs\nex-smokers", "Ex-smokers\nvs\nsmokers"))

g <- ggplot(data = boxplot_subset, aes(Variable, logFC)) + 
  geom_violin(aes(fill = Variable),
              col = "black") +
  geom_boxplot(col = "black",
               outlier.shape = NA,
               notch = T,
               width = 0.25) +
  theme_bw() + ggtitle(paste0("p-value = ", p_val)) +
  xlab("") +
  ylab("logFC") +
  scale_fill_manual(values = c("#88CCEE",  "#CC6677")) +
  stat_summary(fun.data = get_box_stats, fun.args = list(upper_limit = max(boxplot_subset$logFC) * 1.15), 
               geom = "text", hjust = 0.5, vjust = 0.9, size = 4) +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13, colour = "black"), 
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "grey"),
        legend.position = "none", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 13))

pdf(file = "figures/figure_s8/logFC_methylation_subset.pdf", w = 3.5, h = 3.5)
g
dev.off()

p_val <- test_expression$p.value
p_val <- formatC(p_val, format = "e", digits = 2)

boxplot_expression$Variable <- factor(boxplot_expression$Variable, levels=c("Never smokers\nvs\nex-smokers", "Ex-smokers\nvs\nsmokers"))

g <- ggplot(data = boxplot_expression, aes(Variable, logFC)) + 
  geom_violin(aes(fill = Variable),
              col = "black") +
  geom_boxplot(col = "black",
               outlier.shape = NA,
               notch = T,
               width = 0.25) +
  theme_bw() + ggtitle(paste0("p-value = ", p_val)) +
  xlab("") +
  ylab("logFC") +
  scale_fill_manual(values = c("#88CCEE",  "#CC6677")) +
  stat_summary(fun.data = get_box_stats, fun.args = list(upper_limit = max(boxplot_expression$logFC) * 1.15), 
               geom = "text", hjust = 0.5, vjust = 0.9, size = 4) +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13, colour = "black"), 
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "grey"),
        legend.position = "none", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 13))

pdf(file = "figures/figure_s8/logFC_expression_subset.pdf", w = 3.5, h = 3.5)
g
dev.off()


# Figure S8E

#Methylation figure
methylation_data <- readRDS(file = "data/methylation_ml.rds")

pdf("figures/figure_S8/feature_importance.pdf", w = 8, h = 8)

ggplot(methylation_data$feature_importance, aes(y = Gain, x = reorder(Gene_name, Gain))) +
  geom_bar(stat="identity", fill = safe_colorblind_palette[6]) +
  coord_flip() +
  scale_y_continuous(breaks = seq(0, 0.6, by = 0.1)) +
  theme_classic() +
  xlab("") +
  theme(axis.title.x = element_text(size = 17, colour = "black", margin = margin(t = 20)),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.text.x = element_text(size = 16, colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 20))

dev.off()
