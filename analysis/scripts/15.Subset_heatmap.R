#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Methylation results analysis on reversibility
# @software version: R=4.2.2

#Set path 
# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")

tissues <- c("BreastMammaryTissue", "ColonTransverse", "Lung", "MuscleSkeletal", "Ovary", "Prostate", "Testis", "WholeBlood")
metadata <- list()
table <- data.frame("tissue"=1, "DEGs"=1, "DMPs"=1)

for(tissue in tissues){
  print(tissue)
  metadata_s <- readRDS(paste0("tissues/", tissue, "/metadata_subset.rds"))
  metadata[[tissue]] <- metadata_s
  methylation <- readRDS(paste0("tissues/", tissue, "/DML_results_subset.rds"))
  expression <- readRDS(paste0("tissues/", tissue, "/voom_limma_results_subset.rds"))
  table <- rbind(table, c(tissue, sum(expression$Smoking2$adj.P.Val<0.05), sum(results$Smoking2$adj.P.Val<0.05)))
}
table <- table[-1,]
table <- as.data.frame(table)


#Plotting
library(RColorBrewer) 
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))

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
                   row_names_gp = gpar(fontsize = 9),
                   column_names_gp = gpar(fontsize = 10),
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(prettyNum(to_plot[,2:3][i, j], big.mark = ","), x, y, gp = gpar(fontsize = 10))},
                   heatmap_legend_param = list(direction = "horizontal")
)

pdf(file = "figures/figure_6/heatmap.pdf", w = 2.8, h = 4)
draw(ht_DSGs, heatmap_legend_side = "bottom")
dev.off()