#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to plot expression downsampling
# @software version: R=4.2.2


#Set path 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")


load("Downsampling/downsampling_expression.Rdata")
medians <- round(medians)

#
#
# #plot
# library(ComplexHeatmap)
# library(RColorBrewer)
# Heatmap(apply(to_plot, 2,function(x) x/max(x,na.rm=T)), col = brewer.pal(9,"BuPu")[1:7],
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(prettyNum(to_plot[i, j], big.mark = ","), x, y, gp = gpar(fontsize = 10))},
#         na_col = "white",
#         cluster_rows = F,
#         cluster_columns = F,
#         name = "DE signal",
#         row_names_side = "left",)

tissue_info <- read.csv("data/public/tissue_abreviation.txt")
rownames(medians) <- sapply(rownames(medians), function(name) tissue_info$Name[tissue_info$X==name])

library(ComplexHeatmap)
library(RColorBrewer)
medians <- medians[rev(order(medians[,3], medians[,4], medians[,5], na.last = F)),]
no_nas <- medians
no_nas <- apply(apply(apply(no_nas, 2, as.numeric), 2, round), 2, prettyNum, big.mark=",")
no_nas[no_nas=="NA"] <- ""

colnames(medians) <- as.numeric(colnames(medians))*2
pdf("../figures/figures/figure_s1/downsampling.pdf", height = 6.5, width = 4.8)
Heatmap(apply(medians, 2,function(x) x/max(x,na.rm=T)), col = brewer.pal(9,"BuPu")[c(1,4:7)],
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(no_nas[i, j], x, y, gp = gpar(fontsize = 9))},
        na_col = "white",
        cluster_rows = F,
        column_names_side = "bottom",
        cluster_columns = F,
        name = "Nº DEGs\n(Scaled)",
        column_names_rot = 0,
        column_title_side = "bottom",
        column_title = "Number of samples",
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        row_names_side = "left",)
dev.off()
#Sort tissues


#To main:
library(reshape2)
library(ggplot2)

df <- melt(medians)
df$Var2 <- as.factor(df$Var2)
cols <- tissue_info$color
names(cols) <- tissue_info$Name

g <- ggplot(df, aes(Var2, value, col=Var1)) + geom_point() +
  scale_color_manual(values=cols) +
  geom_line(aes(group = Var1)) + theme_bw() +
  theme(legend.title = element_blank(),
        axis.text = element_text(colour="black", size=13),
        legend.text = element_text(colour="black", size=13),
        axis.title = element_text(size=16),
        # legend.spacing.y = unit(-0.05, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth=1)) +
  ylab("Number of DEGs") + xlab("Number of samples")

pdf("../figures/figures/figure_s1/downsampling_main.pdf", height = 2.5, width = 10.2)
g
dev.off()
