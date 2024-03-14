# @Author: Rogerio Ribeiro and Jose Miguel Ramirez
# @E-mail: rogerio.e.ramos.ribeiro@gmail.com
# @Date:   2023-08-14
# @Description: Generate publication Figure 2 (Study overview and RNA-seq results)

#Set path 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Load libs
library(reshape2)
library(tidyverse)
library(DOSE)
library(factoextra)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggrepel)

#Load colours 
source("safe_colourblind_pallete.R")
tissue.data <- read.csv(file = "tissue_abreviation.txt")


# Load data
figure_data <- readRDS(file = "data/degs_summary.rds")
enrichement_analysis_data <- readRDS(file = "data/enrichement_results.summary.rds")

splicing_data <- readRDS(file = "data/splicing_data.rds")
splicing_data_2 <- readRDS(file = "data/splicing_data2.rds")


for(i in 1:length(figure_data)) assign(names(figure_data)[i], figure_data[[i]])
for(i in 1:length(enrichement_analysis_data)) assign(names(enrichement_analysis_data)[i], enrichement_analysis_data[[i]])
for(i in 1:length(splicing_data)) assign(names(splicing_data)[i], splicing_data[[i]])
for(i in 1:length(splicing_data_2)) assign(names(splicing_data_2)[i], splicing_data_2[[i]])


# Figure 2A 
DEGS.per.tissue <- DEGS.per.tissue %>% 
  as.data.frame() %>% 
  merge(tissue.data, by.x = 0, by.y = "X")


metadata.filtered <- metadata[DEGS.per.tissue$tissue]


if (all(DEGS.per.tissue$tissue == names(metadata.filtered))) names(metadata.filtered) = DEGS.per.tissue$TISSUENAMEABREV

#Define order by the number of DEGs 
DEGS.per.tissue <- DEGS.per.tissue  %>% 
  arrange(-Total)

tissues <- DEGS.per.tissue$TISSUENAMEABREV

row_ha_left <- HeatmapAnnotation("Number of\nsamples" = anno_barplot(t(sapply(tissues, function(tissue) table(metadata.filtered[[tissue]]$Smoking))),
                                                          gp = gpar(fill = safe_colorblind_palette[c(1,3,2)]),
                                                          border=F),
                                                         gap = unit(0.25,"cm"),
                                                         show_legend = T, 
                                                         show_annotation_name = T,
                                                         annotation_name_rot = 90,
                                                         annotation_name_gp = gpar(fontsize = 10),
                                                         which = "row")


# ht_DEGs <- Heatmap(apply(DEGS.per.tissue[4], 2,function(x) x/max(x,na.rm=T)),
names(DEGS.per.tissue)[4] <- "Smoking-\nDEGs"
ht_DEGs <- Heatmap(as.matrix(DEGS.per.tissue[4]),
                   col = brewer.pal(9,"BuPu")[1:7],
                   na_col = "white",
                   cluster_rows = F,
                   cluster_columns = F,
                   name = "Number of DEGs",
                   row_names_side = "left",
                   row_labels = DEGS.per.tissue$Name,
                   column_names_side = "bottom",
                   left_annotation = row_ha_left,
                   row_names_gp = gpar(fontsize = 9),
                   column_names_gp = gpar(fontsize = 10),
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(prettyNum(DEGS.per.tissue[4][i, j], big.mark = ","), x, y, gp = gpar(fontsize = 10))},
                   heatmap_legend_param = list(direction = "vertical", title_position = "topcenter")
)




pdf(file = "figures/figure_2/figure_2a", w = 3.9, h = 9)
draw(ht_DEGs,
     heatmap_legend_side = "right")
dev.off()


# Figure 2B

load("data/downsampling_expression.Rdata")
medians <- round(medians)

df <- melt(medians)
df$Var2 <- as.factor(df$Var2)
cols <- tissue.data$color
names(cols) <- tissue.data$Name

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

pdf("figures/figure_2/figure_2b.pdf", height = 2.5, width = 9.5)
g
dev.off()


# Figure 2C
tissue.per.gene$direction <- factor(tissue.per.gene$direction, levels = c("Upregulated", "Downregulated"))

#Adding labels
top <- tissue.per.gene[1:8,]
top$label <- sapply(top$ensembl, function(gene) unique(DEGS$gene_name[gsub("\\..*","", rownames(DEGS))==gene]))

pdf("figures/figure_2/figure_2c.pdf", w = 3.5, h = 4)
ggplot(tissue.per.gene, aes(y = n_tissue, x = direction, color = direction)) + #aes(label = label) 
  geom_jitter(width = .35, height = 0, size = 3, alpha=0.8) +
  xlab("") + ylab("Number of Tissues") +
  scale_color_manual(name = "",values=safe_colorblind_palette[c(2,1)]) +
  geom_text_repel(data=top, aes(y = n_tissue, x = direction, label=label), col="black", size=4.5) + #The names will be correctly adjusted in inkscape manually
  #geom_text_repel(size = 5) + 
  theme_bw() + 
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        axis.text.x= element_text(size = 13, colour = "black"),
        axis.ticks.x=element_blank(), 
        axis.text.y =  element_text(size = 18, colour = "black"), 
        axis.title.x =  element_blank(),
        axis.title.y =  element_text(size = 20, colour = "black"), 
        legend.position = "none")
dev.off()



#Figure 2D----
terms_in_n_tissues.up <- tissue_per_terms.up %>% 
  group_by(n_tissues) %>% 
  summarise(n_terms_per_tissue = n()) %>% 
  mutate(group = factor(cut(n_tissues, c(0,1,2,5), labels = c("1", "2", "3-4")))) %>%
  group_by(group) %>% 
  summarise( n_terms_per_tissue = sum( n_terms_per_tissue)) %>%
  mutate(x = "Up\nTerms")


terms_in_n_tissues.down <- tissue_per_terms.down %>% 
  group_by(n_tissues) %>% 
  summarise(n_terms_per_tissue = n()) %>% 
  mutate(group = factor(cut(n_tissues, c(0,1,2,5), labels = c("1", "2", "3-4")))) %>%
  group_by(group) %>% 
  summarise( n_terms_per_tissue = sum( n_terms_per_tissue)) %>%
  mutate(x = "Down\nTerms")

terms_in_n_tissues <- rbind(terms_in_n_tissues.up, terms_in_n_tissues.down) %>% 
  mutate(x = factor(x, levels = c("Up\nTerms", "Down\nTerms")))

#Based on the colours from safe_color_blind_pallete
colour_pallete <- c( "#DCF0FA", "#F0D1D6", "#C9E7F7", "#DB94A0", "#88CCEE", "#CC6677") #Order is important

terms_in_n_tissues$group.2 <- factor(paste0(terms_in_n_tissues$group, terms_in_n_tissues$x))

pdf("figures/figure_2/figure_2D.pdf", w = 3.4, h = 4)
ggplot(terms_in_n_tissues, aes(fill = group.2, y=n_terms_per_tissue, x=x, label = n_terms_per_tissue)) +
  geom_col(position="stack") + ylab("Number of Terms") + 
  xlab("") +
  scale_fill_manual(name = "Nº Tissues", values = colour_pallete) +
  guides(fill = guide_legend()) + 
  geom_text(position=position_stack(vjust = .5), size = 5, colour = c("black", "black", "black", "black", "black", "black")) + 
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(size = 14, colour = "black"),
        axis.ticks.x=element_blank(), 
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"), 
        legend.text = element_text(size = 13, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.7)) +
  labs(fill='Number of tissues')
dev.off()

# Figure 2E ----

enrichement_analysis_orsum <- enrichement_analysis_data$orsum_top_terms %>% 
  filter(`Representing term rank` <=5)

#Modify tissue and term order to highlight the shared terms. Remove similar terms in order to simplify the image. Also remove some terms that do not make sense
#i.e brain related stuff in lung

tissue_order <- c("Lung", "Thyroid", "Pancreas", "Skin Lower Leg", 
                  "Artery Aorta", "Artery Tibial",
                  "Adipose Subcutaneous",  "Stomach",
                  "Esophagus Mucosa")

# Only keep a subset of terms (common enriched terms)
terms2keep <- enrichement_analysis_orsum %>% 
  group_by(`Representing term name`) %>% 
  summarise(n_tissue = n()) %>% 
  filter(n_tissue >= 2) %>% 
  pull(`Representing term name`)


enrichement_analysis_orsum <- enrichement_analysis_orsum %>% 
  filter(`Representing term name` %in% terms2keep) %>%
  mutate(`Representing term name` = as.character(`Representing term name`)) %>%
  mutate(TISSUENAMEABREV =  factor(TISSUENAMEABREV, levels = tissue_order[length(tissue_order):1])) %>%
  arrange(direction, TISSUENAMEABREV) %>% 
  mutate(`Representing term name` = factor(`Representing term name`, levels = unique(`Representing term name`))) 


pdf("figures/figure_2/figure_2E.pdf", w = 11, h = 8.5)
ggplot(enrichement_analysis_orsum, aes(y = TISSUENAMEABREV, x = `Representing term name`, color = direction)) +
  geom_point(size = 4) +
  scale_color_manual(name = "Gene Set", values = safe_colorblind_palette[c(2,1)]) +
  theme_dose(10) +
  xlab("") + ylab("") + 
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, angle = 80, hjust = 1, vjust = 1), 
        axis.text.y = element_text(size = 16), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 22), 
        legend.position = "top", 
        legend.direction="horizontal")
dev.off()


# Figure 2F----

to_plot <- splicing_table[splicing_table$DSE>0,]
to_plot <- to_plot[rev(order(to_plot$DSE)),]

to_plot$tissue <- sapply(to_plot$tissue, function(name) tissue_info$Name[tissue_info$TISSUENAMEABREV==name])
names(metadata) <- sapply(names(metadata), function(name) tissue_info$Name[tissue_info$tissue==name])
tissues <- to_plot$tissue


ht_DSGs <- Heatmap(apply(to_plot[2], 2,function(x) x/max(x,na.rm=T)),
                   col = brewer.pal(9,"BuPu")[1:7],
                   na_col = "white",
                   cluster_rows = F,
                   cluster_columns = F,
                   name = "DS signal",
                   row_names_side = "left",
                   row_labels = to_plot$tissue,
                   column_names_side = "bottom",
                   left_annotation = row_ha_left,
                   row_names_gp = gpar(fontsize = 9),
                   column_names_gp = gpar(fontsize = 10),
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(prettyNum(to_plot[2][i, j], big.mark = ","), x, y, gp = gpar(fontsize = 10))},
                   heatmap_legend_param = list(direction = "horizontal")
)

pdf(file = "figures/figure_2/figure_2F.pdf", w = 2.6, h = 3.8)
draw(ht_DSGs, heatmap_legend_side = "bottom")
dev.off()

### Figure 2F
get_box_stats <- function(y, upper_limit ) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "N =", length(y), "\n"
    )
  ))
}

data_cldn7 <- data_cldn7[data_cldn7$Smoking!=1,]
data_cldn7$Smoking <- droplevels(data_cldn7$Smoking)
# data_cldn7$title <- "CLDN7 (Lung)\nSkipping exon:chr17:7260726-7260642"
data_cldn7$title <- "Skipping exon\nchr17:7260726-7260642"

g4 <- ggplot(data_cldn7, aes(Smoking, V1)) +
  geom_violin(aes(fill = Smoking),
              col = "black") +
  geom_boxplot(col = "black",
               outlier.shape = NA,
               notch = T,
               width = 0.25) +
  geom_jitter(col = "black",
              alpha = 0.1,
              size = 0.8) +
  theme_minimal() +
  scale_fill_manual(values=c("#88CCEE", "#CC6677")) +
  ylab("Percentage Spliced In") + ggtitle("CLDN7 (Lung)") +
  facet_grid(. ~ title) + # Add the title strip
  theme(axis.title.y = element_text(margin = margin(r = 20), size = 18),
        axis.text.y = element_text(size = 15, colour="black"),
        axis.text.x = element_text(size = 15, colour="black"),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "grey"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 20),
        strip.background = element_rect(fill="#88CCEE"),
        strip.text = element_text(size = 16)) +
  stat_summary(aes(x = Smoking, y = V1),
               fun.data = get_box_stats, fun.args = list(upper_limit = max(data_cldn7$V1) * 1.08),
               geom = "text",
               hjust = 0.5, vjust = 0.9, size = 5) +
  scale_x_discrete(labels=c("Never smokers", "Smokers")) + xlab("")
g4

ggsave("figures/figure_2/figure_2g_cldn7_boxplot.pdf", g4, device="pdf", width = 4, height = 3.5)

