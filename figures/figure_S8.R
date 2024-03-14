# @Author: Rogerio Ribeiro and Jose Miguel Ramirez
# @E-mail: rogerio.e.ramos.ribeiro@gmail.com and jose.ramirez1@bsc.es
# @Date:   2022-03-28
# @Description: Generate publication Figure s7

#Set path 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Load libs
library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggrepel)


# Load data 
source("safe_colourblind_pallete.R")
tissue.data <- read.csv(file = "tissue_abreviation.txt")


if (!dir.exists("figures/figure_s6/")){
  dir.create("figures/figure_s6/")
}

#Load figure data 
degs <- readRDS(file = "data/degs_summary.rds")
reverse <- readRDS(file = "data/reversible_genes.rds")
ml_data <- readRDS(file = "data/methylation_data.rds")
splicing_data <- readRDS(file = "data/splicing_data2.rds")

#Figure S6A ----

degs_ns_ex <- degs$DEGS.per.tissue_ex_ns
degs_s_ex <- degs$DEGS.per.tissue_ex_s

#Total number of DEGS in each comparison
degs_ex <- cbind(degs_ns_ex[,3], degs_s_ex[,3])
degs_ex <- merge(degs_ex, tissue.data, by.x = 0, by.y = "X") %>%
  arrange(-V2)

ht_DEGs <- Heatmap(degs_ex[,c(2,3)],
                   col = brewer.pal(9,"BuPu")[1:7],
                   na_col = "white",
                   cluster_rows = F,
                   cluster_columns = F,
                   name = "Nº DEGs\n(scaled)",
                   row_names_side = "left",
                   row_labels = degs_ex$Name,
                   column_names_side = "bottom",
                   column_labels=  c("Never vs Ex", "Ex vs Smoker"), 
                   row_names_gp = gpar(fontsize = 10),
                   column_names_gp = gpar(fontsize = 10),
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(prettyNum(degs_ex[,c(2,3)][i, j], big.mark = ","), x, y, gp = gpar(fontsize = 10))},
                   heatmap_legend_param = list(direction = "vertical")
)

pdf(file = "figures/figure_s6/ex_smoker_DEGS.pdf", w = 3.5, h = 8)
draw(ht_DEGs,
     heatmap_legend_side = "right")
dev.off()

# FIgure S6B ----
splicing_data <- 

#Figure S6B ----
n_degs_ns_vs_ex <- lapply(ml_data, function(x) sum(x[["Smoking1"]]$adj.P.Val < 0.05))
n_degs_ex_vs_s <- lapply(ml_data, function(x) sum(x[["Smoking1-Smoking2"]]$adj.P.Val < 0.05))

n_dml <- cbind(n_degs_ns_vs_ex, n_degs_ex_vs_s)
n_dml <- data.frame(n_dml)
n_dml[,1] <- as.numeric(n_dml[,1])
n_dml[,2] <- as.numeric(n_dml[,2])


ht_DMP <- Heatmap(data.frame(n_dml),
                   col = brewer.pal(9,"BuPu")[1:7],
                   na_col = "white",
                   cluster_rows = F,
                   cluster_columns = F,
                   name = "Nº Positions",
                   row_names_side = "left",
                   row_labels = row.names(n_dml),
                   column_names_side = "bottom",
                   column_labels=  c("Never vs Ex", "Ex vs Smoker"), 
                   row_names_gp = gpar(fontsize = 10),
                   column_names_gp = gpar(fontsize = 10),
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(prettyNum(n_dml[i, j], big.mark = ","), x, y, gp = gpar(fontsize = 10))},
                   heatmap_legend_param = list(direction = "vertical")
)

pdf(file = "figures/figure_s6/ex_smoker_DMP.pdf", w = 3.5, h = 2.6)
draw(ht_DMP,
     heatmap_legend_side = "right")
dev.off()

## Figure S6C ----

figure_data <- readRDS(file = "data/figure_s9.rds") #list of expression_barplot, expression_pirate, splicing_barplot and splicing_pirate
for(i in 1:length(figure_data)) assign(names(figure_data)[i], figure_data[[i]])

pushViewport(viewport(gp = gpar(fontfamily = "Helvetica")))

pdf("figures/figure_S6/Heatmaps.pdf",  height = 5.7 , w = 5)
# rownames(data) <- gsub("\\(.*","",rownames(data))
Heatmap(as.matrix(data),
        col = brewer.pal(9,"BuPu")[c(1:7)],
        cluster_rows = F,
        cluster_columns = F,
        name = "Number of DSEs",
        row_names_side = "left",
        column_names_side = "bottom",
        row_names_gp = gpar(fontsize = 12),
        column_names_gp = gpar(fontsize = 12),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(data[i, j], x, y, gp = gpar(fontsize = 11))},
        heatmap_legend_param = list(direction = "vertical"))
dev.off()



# Figure S8C-d

#Boxplot showing that ex-smokers are closer to never smokers:
get_box_stats <- function(y, upper_limit) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "N =", length(y), "\n"
    )
  ))
}


figure_data <- readRDS(file = "data/figure_s9_logFC_subset.rds") #list of expression_barplot, expression_pirate, splicing_barplot and splicing_pirate
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

pdf(file = "figures/figure_s6//logFC_methylation_subset.pdf", w = 4.5, h = 3.5)
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

pdf(file = "figures/figure_s6/logFC_expression_subset.pdf", w = 4.5, h = 3.5)
g
dev.off()
