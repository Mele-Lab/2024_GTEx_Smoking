# @Author: Rogerio Ribeiro and Jose Miguel Ramirez
# @E-mail: rogerio.e.ramos.ribeiro@gmail.com
# @Date:   2023-08-14
# @Description: Generate publication Figure S1

#Set path 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


#Load libs
library(tidyverse)
library(factoextra)
library(UpSetR)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)

# Load data 
source("safe_colourblind_pallete.R")
tissue.data <- read.csv(file = "tissue_abreviation.txt")

# Load figure data
data_prev <- readRDS(file = "data/reproduction_analysis.rds")
data_degs <- readRDS(file = "data/degs_summary.rds")
#metadata <- readRDS(file = "data/metadata.rds")

for(i in 1:length(data_prev)) assign(names(data_prev)[i], data_prev[[i]])
for(i in 1:length(data_degs)) assign(names(data_degs)[i], data_degs[[i]])

# Since this figure is simple I can use the ggpubr to generate the final pannel

#Figure S1A ---- 
## Heatmap of the up + down genes 

DEGS.per.tissue <- DEGS.per.tissue %>% 
  as.data.frame() %>% 
  merge(tissue.data, by.x = 0, by.y = "X") %>%
  arrange(-Total)

genes <- c(DEGS.per.tissue$Up, DEGS.per.tissue$Down)
# genes_Scale <- genes / max(genes)
# 
# genes_scaled <- matrix(genes_Scale, byrow = F, ncol = 2)

# ht_DEGs <- Heatmap(genes_scaled,
ht_DEGs <- Heatmap(matrix(genes, byrow = F, ncol = 2),
                   col = brewer.pal(9,"BuPu")[1:7],
                   na_col = "white",
                   cluster_rows = F,
                   cluster_columns = F,
                   name = "Number of DEGs",
                   row_names_side = "left",
                   row_labels = DEGS.per.tissue$Name,
                   column_names_side = "bottom",
                   column_names_rot = 0, #I will rotate in inkscape
                   column_names_centered = T,
                   column_labels = c("Up", "Down"),
                   row_names_gp = gpar(fontsize = 10),
                   column_names_gp = gpar(fontsize = 10),
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(prettyNum(DEGS.per.tissue[2:3][i, j], big.mark = ","), x, y, gp = gpar(fontsize = 10))},
                   heatmap_legend_param = list(direction = "vertical", title_position = "topcenter")
)

pdf(file = "figures/figure_s1/figure_S1A.pdf", height = 6.5, width = 3.5)
draw(ht_DEGs,
     heatmap_legend_side = "right")
dev.off()


#Figure S1B----
load("data/downsampling_expression.Rdata")
medians <- round(medians)
rownames(medians) <- sapply(rownames(medians), function(name) tissue.data$Name[tissue.data$X==name])

medians <- medians[rev(order(medians[,3], medians[,4], medians[,5], na.last = F)),]
no_nas <- medians
no_nas <- apply(apply(apply(no_nas, 2, as.numeric), 2, round), 2, prettyNum, big.mark=",")
no_nas[no_nas=="NA"] <- ""

colnames(medians) <- as.numeric(colnames(medians))*2
pdf("figures/figure_s1/figure_S1B.pdf", height = 6.5, width = 4.8)
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


# Figure S1C ----

degs_summary_per_logFC <- readRDS("data/degs_summary_per_logFC.rds")

DEGS.per.tissue  <- dcast(degs_summary_per_logFC, tissue ~ logFC, value.var = "n_degs", fill = 0)
DEGS.per.tissue <- merge(DEGS.per.tissue, tissue.data,by.x = "tissue", by.y = "SMTSDNOSEP")
DEGS.per.tissue <- merge(DEGS.per.tissue, n_samples_per_tissue,by = "tissue")


DEGS.per.tissue <- DEGS.per.tissue  %>% 
  arrange(-`0`)

tissues <- DEGS.per.tissue$tissue

names(DEGS.per.tissue)[2:6] <- c("logFC 0", "logFC 0.5", "logFC 0.7", "logFC 1", "logFC 1.5")

DEGS.per.tissue_scaled <- DEGS.per.tissue

DEGS.per.tissue_scaled[2:6] <- apply(DEGS.per.tissue[2:6], 2, scale)
ht_DEGs <- Heatmap(as.matrix(DEGS.per.tissue_scaled[2:6]),
                   col = brewer.pal(9,"BuPu")[1:7],
                   na_col = "white",
                   show_heatmap_legend = F,
                   cluster_rows = F,
                   cluster_columns = F,
                   name = "Number of DEGs",
                   row_names_side = "left",
                   row_labels = DEGS.per.tissue$SMTSCUSTOM,
                   column_names_side = "bottom",
                   row_names_gp = gpar(fontsize = 9),
                   column_names_gp = gpar(fontsize = 10),
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(prettyNum(DEGS.per.tissue[2:6][i, j], big.mark = ",",), x, y, gp = gpar(fontsize = 10))},
                   heatmap_legend_param = list(direction = "vertical", title_position = "topcenter")
)


pdf("figures/figure_s1/ht_degs_logFC.pdf", w = 6, h = 10)
plot(ht_DEGs)
dev.off()

# Figure S1d-----

colnames(reproduced_data) <- c("Replicated", "Total", "study", "tissue", "OR", "ConfIntL", "ConfIntUpper", "p_vals")


reproduced_data[6,4] <- c("Adipose Subcutaneous")

reproduced_data<- reproduced_data %>%
  as_tibble() %>%
  merge( tissue.data, by.x = "tissue", by.y = "TISSUENAMEABREV") %>%
  mutate_at(.vars = 5, as.numeric) %>%
  mutate(study_name = paste0(tissue, "; ", study)) %>%
  dplyr::select(study_name, OR, color, ConfIntL, ConfIntUpper, Total, p_vals)

reproduced_data <- reproduced_data %>% 
  mutate_at(
    .vars = c(2,4:6), 
    as.numeric
  ) %>% 
  mutate(log_or = log(OR)) %>% 
  mutate(log_ConfIntL = log(ConfIntL)) %>% 
  mutate(log_ConfIntUpper = log(ConfIntUpper))

reproduced_data$study_name <- c("Adipose subcutaneous: Tsai et al. 2018"," Lung: Bossé et al. 2012",
                                "Lung: Landi et al. 2008", "Lung: Pintarelli et al. 2019",
                                "Whole blood: Huan et al. 2016", "Whole blood: Vink et al. 2017")

reproduced_data <- reproduced_data[c(1,5,6,2,3,4),]
reproduced_data$study_name <- factor(reproduced_data$study_name, levels=c(reproduced_data$study_name))
reproduced_data$p_vals <- as.numeric(reproduced_data$p_vals)
a <- ggplot(reproduced_data, aes(x = study_name, y = log_or, colour = color)) +
  geom_errorbar(aes(ymin=log_ConfIntL, ymax=log_ConfIntUpper), width=.1) +
  geom_point(aes(size = -log10(p_vals))) + ylim(0,6.5) +
  scale_color_manual(values = c("#CC6677", "#88CCEE", "#DDCC77"), guide = "none") +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 9, vjust = 0.5, colour="black"),
        axis.text.y = element_text(size = 8, colour="black"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10, colour="black")) +
  ylab("Log(Odds ratio)") +
  coord_flip() +
  geom_hline(yintercept=0, linetype="dashed") +
  guides(size=guide_legend(title="-log10(p-value)"))



b <- ggplot(reproduced_data, aes(y = Total, x = study_name, fill = color)) +
  geom_col() +
  scale_fill_manual(values = c("#CC6677", "#88CCEE", "#DDCC77")) +
  coord_flip() + theme_bw() +
  scale_y_continuous(breaks=c(0, 1000, 2000), limits = c(0, 2300)) +
  geom_text(aes(label=Total), size=3, hjust = -0.1) +
  ylab("Number of DEGs") + xlab("") + 
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        # axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 10, colour="black"),
        axis.text.x = element_text(size = 9, colour="black"), 
        legend.position = "none")

pdf(file =  "figures/figure_s1/figure_S1C.pdf", w = 6, h = 2)
ggpubr::ggarrange(a, b, widths = c(0.65,0.35), common.legend = T, legend = "right")
dev.off()


#Figure S1e ----
#NOTE: I am grouping by gene name because there are no duplicates names in the subset of genes I want to plot
DEGS.recurrent <- DEGS %>%
  group_by(gene_name) %>% 
  summarise(n_tissues = n()) %>% 
  filter(n_tissues >= 9)

# SLITRK2 is up in 9 tissues and down in 1.

DEGS_recurrent_across_tissues <- DEGS %>% 
  filter(gene_name %in% DEGS.recurrent$gene_name) %>% 
  merge(tissue.data, by = "tissue") %>% 
  arrange(TISSUENAMEABREV) %>%
  # mutate(gene_name = factor(gene_name, levels = unique(gene_name)[length(gene_name):1])) %>%
  mutate(gene_name = factor(gene_name, levels = DEGS.recurrent$gene_name[order(DEGS.recurrent$n_tissues)])) %>%
  mutate(TISSUENAMEABREV = factor(TISSUENAMEABREV, levels = unique(TISSUENAMEABREV)))

max <- max(abs(DEGS_recurrent_across_tissues$logFC))

#Only first word of the tissue is uppercase 

figureS1D <- ggplot(DEGS_recurrent_across_tissues, aes(y = Name, x = gene_name)) + 
  geom_tile(aes(fill = logFC), colour = "black") + 
  scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red", space = "Lab",  #Using these coulours instead of the colorblind pallete as we can distinguish better between up and down
                       limits = c(-min(max), max((max)))) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 20, angle = 90, colour = "black", hjust = 1, vjust = 0.3), 
        # axis.text.y =  element_blank(), 
        axis.text.y =  element_text(size = 18, colour = "black"), 
        legend.position = "right", 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.5, "cm"),
        axis.ticks = element_blank(), 
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "grey")) + 
  labs(x="", y = "") + 
  coord_fixed(ratio = 0.6) + 
  coord_flip()

figureS1D

ggsave(figureS1D, file = "figures/figure_s1/figure_S1D.pdf", w = 12, h = 5.8, device="pdf")
