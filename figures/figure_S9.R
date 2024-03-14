# @Author: Rogerio Ribeiro and Jose Miguel Ramirez
# @E-mail: rogerio.e.ramos.ribeiro@gmail.com and jose.ramirez1@bsc.es
# @Date:   2022-03-28
# @Description: Generate publication Figure S9 (Machine learning in gexp)


library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)

#Set working dir 
#setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Colours
source("safe_colourblind_pallete.R")

#Tissue names
tissue_names <- read.csv(file = "tissue_abreviation.txt")

# Load data
figure_data <- readRDS(file = "data/machine_learning_gene_expression.rds")
for(i in 1:length(figure_data)) assign(names(figure_data)[i], figure_data[[i]])

## Figure S7B - Plot Heatmap with AUC across tissues

model_classification <- merge(model_classification, tissue_names, by= "tissue") %>% 
  arrange(-value)

ht_ml <- Heatmap(t(model_classification$value),
                 col = brewer.pal(9,"BuPu")[1:7],
                 na_col = "white",
                 cluster_rows = F,
                 cluster_columns = F,
                 name = "AUC",
                 column_names_side = "bottom",
                 column_labels = model_classification$TISSUENAMEABREV,
                 row_names_side = "right",
                 row_names_gp = gpar(fontsize = 13),
                 column_names_gp = gpar(fontsize = 15),
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(round(t(model_classification[3])[i, j], 2), x, y, gp = gpar(fontsize = 10))},
                 heatmap_legend_param = list(direction = "vertical")
)

pdf(file = "figures/figure_S10/ml_classification.pdf", w = 15, h = 3.5)
draw(ht_ml,
     heatmap_legend_side = "right")
dev.off()
