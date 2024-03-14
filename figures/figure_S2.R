# @Author: Rogerio Ribeiro and Jose Miguel Ramirez
# @E-mail: rogerio.e.ramos.ribeiro@gmail.com
# @Date:   2022-05-09
# @Description: Generate publication Figure S2

#Set path 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


#Load libs
library(tidyverse)
library(UpSetR)
library(DOSE)
library(igraph)

# Load data 
source("safe_colourblind_pallete.R")
tissue.data <- read.csv(file = "tissue_abreviation.txt")

#Load Enrichment analysis data 
data <- readRDS(file = "data/enrichement_results.summary.rds")

for(i in 1:length(data)) assign(names(data)[i], data[[i]])

if (!dir.exists("figures/figure_s2/")){
  dir.create("figures/figure_s2/")
}


tissue.data <- read.csv(file = "tissue_abreviation.txt")

# Figure S2A ----

lung <- orsum_top_terms  %>%
  filter(`Representing term rank` <= 5) %>%
  filter(tissue  == "Lung")  %>%
  mutate(direction = factor(direction, levels = c("Upregulated", "Downregulated"))) %>% 
  arrange(direction) %>% 
  mutate(term = factor(`Representing term name`)) %>% 
  mutate(term = factor(term, levels = term[length(term):1]))

pdf(file = "figures/figure_s2/lung_GO_BP.pdf", width = 9, height = 8)
ggplot(lung, aes(y = term, x = `Representing term size`, fill = direction)) + 
  geom_bar(stat="identity", width = 0.7) +
  scale_fill_manual(name = "Gene Set (GO:BP)", values = safe_colorblind_palette[c(2,1)]) +
  theme_bw() + ggtitle("Gene set (GO:BP)") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=22, colour = "black"),
        axis.text.x = element_text(size = 14, colour = "black"), 
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.title = element_text(size = 18, colour = "black"),
        legend.position = "none") + 
  xlab("Number of annotated Genes\n(Lung)") + ylab("") + 
  scale_x_continuous(expand = c(0,0))
dev.off()


# Figure S2B ----

KEGG_top_terms <- KEGG_top_terms %>% merge(tissue.data, by.x = "sample", by.y = "tissue")

pdf(file = "figures/figure_s2/KEGG_enrichement.pdf",width = 9, height = 8.5)
KEGG_top_terms  %>% 
  mutate(Direction = factor(Direction, levels = c("Upregulated", "Downregulated"))) %>%
  arrange(Direction, sample) %>%
  mutate(Description = factor(Description, levels = unique(Description)[length(Description):1]))%>%
  ggplot(aes(x=TISSUENAMEABREV, y=Description, color = Direction)) + 
  geom_point(size = 4) +
  scale_color_manual(name = "Gene set (KEGG)", values = safe_colorblind_palette[c(2,1)]) +
  theme_bw() + ggtitle("Gene set (KEGG)") +
  xlab("") + ylab("") + 
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=22, colour = "black"),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, colour = "black"), 
        axis.text.y = element_text(size = 16, colour = "black"), 
        legend.text = element_text(size = 22, colour = "black"), 
        legend.position = "none")
dev.off()



# Plot figure S2C----
thyroid_do_down$tissue <- "Thyroid"
lung_do_up$tissue <- "Lung"

thyroid_do_down$direction <- "Downregulated"
lung_do_up$direction <- "Upregulated"

lung_do_up <- lung_do_up %>% 
  slice_head(n = 15)


do_lung_thyroid <- rbind(lung_do_up, thyroid_do_down) %>% 
  mutate(Description = factor(Description, levels = unique(Description)[length(unique(Description)):1])) %>% 
  mutate(direction = factor(direction, levels = c("Upregulated", "Downregulated")))


pdf(file = "figures/figure_s2/DO_enrichement_lung_thyroid.pdf", width = 8, height = 8.5)
ggplot(do_lung_thyroid, aes(y = Description, x = tissue, colour = direction)) + 
  geom_point(size = 4) +
  scale_color_manual(values = safe_colorblind_palette[c(2,1)]) +
  theme_bw() + ggtitle("Gene set (Disease Ontology)") +
  xlab("") + ylab("") + 
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=22, colour = "black"),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, colour = "black"), 
        axis.text.y = element_text(size = 16, colour = "black"), 
        legend.text = element_text(size = 22, colour = "black"),
        legend.title = element_blank())
dev.off()
