#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Generate publication Figure S7


#Set path 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Loading libraries
library(ggplot2)


#LogFC:
figure_data <- readRDS(file = "data/logFC_pos_neg.rds")
for(i in 1:length(figure_data)) assign(names(figure_data)[i], figure_data[[i]])

get_box_stats <- function(y, upper_limit = max(y) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "n =", length(y)
    )
  ))
}

g <-  ggplot(smoking, aes(dummy, abs(logFC))) + geom_violin(aes(fill=ecpg)) + 
  geom_boxplot(outlier.shape = NA, notch = T,
               width = 0.25) + theme_bw() + 
  stat_summary(fun.data = get_box_stats, fun.args = list(upper_limit = max(abs(smoking$logFC)) * 1.5), 
               geom = "text", hjust = 0.5, vjust = 0.9, size = 3.8) +
  scale_fill_manual(values=c("#88CCEE", "#CC6677"))  +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(colour="black", size=11),
        axis.title = element_text(size=13)) + xlab("")
g

pdf(file = "figures/figure_s5/logFC_pos_neg.pdf", w = 4.5, h = 3)
g
dev.off()


#AHRR correlations:
to_plot_correlations <- readRDS(file = "data/ahrr_correlations.rds")

#Improve annotation based on ChromHMM
table(to_plot_correlations$category) #Only 4 annotated as enhancer and 1 as promoter 
to_plot_correlations$category[to_plot_correlations$region_chromhmm %in% c("Weak transcription", "Strong transcription")] <- "gene_body"
to_plot_correlations$category[to_plot_correlations$region_chromhmm %in% c("Active enhancer", "Weak enhancer", "Genic enhancer")] <- "enhancer"
to_plot_correlations$category[to_plot_correlations$region_chromhmm %in% c("Active TSS", "Flanking TSS")] <- "promoter"
table(to_plot_correlations$category) #18 annotated as enhancer and 4 as promoter 
to_test <- to_plot_correlations[to_plot_correlations$p.adj<0.05,]
table(to_test$category, to_test$cor>0)

#39/49
#7/10

# to_plot_correlations$category[to_plot_correlations$region_chromhmm %in% c("Genic enhancer")] <- "genic enhancer" #3 of the positively correlated with enhancers are genic enhancers
to_plot_correlations$New[to_plot_correlations$New=="TRUE"] <- "Known"
to_plot_correlations$New[to_plot_correlations$New=="FALSE"] <- "New"
table(to_plot_correlations$New[to_plot_correlations$p.adj<0.05 & to_plot_correlations$cor>0])
g <- ggplot(to_plot_correlations) + geom_point(aes(cor, -log10(p.adj), color=category, shape=New), size=3) + theme_bw() + #p value of correlation
  ylab("-log10(p-value)") + geom_hline(yintercept=-log10(0.05)) +
  scale_color_manual(values=c("#CC6677", "#88CCEE", "#DDCC77"), labels=c("Enhancer", "Gene body", "Promoter")) +
  theme(legend.title = element_blank(),
        legend.position="bottom",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.box="vertical") + xlab("Correlation (rho)") +
  scale_shape_manual(values=c(1, 2))
  # scale_shape_manual(values=c(16, 17))
g
ggsave("figures/figure_s7/ahrr_correlations.pdf", g, width = 3.3, height = 4)


MOFAobject.trained <- load("..") # Load the trained MOFA object
var_explained <- get_variance_explained(MOFAobject.trained)
data_var_explained <- var_explained$r2_per_factor$group1

min_r2 <- min(data_var_explained)
max_r2 <- max(data_var_explained)

library(reshape2)
data_var_explained <- melt(data_var_explained, 
                           varnames = c("Factor", "Feature"), value.name = "Value"
                           )
 
data_var_explained$Factor <- factor(data_var_explained$Factor, levels = rev(unique(data_var_explained$Factor)))


pdf("figures/figure_s7/lung_variance_per_modality_per_view.v2.pdf", h = 6, w = 5)
ggplot(data_var_explained, aes(y = Factor, x = Feature)) + 
  geom_tile(aes(fill = Value), color = "black") +
  geom_text(aes(label = round(Value, 3)), size = 4, color = "black") + 
  labs(x = "", y = "", title = "") + 
  scale_fill_gradientn(colors = c("gray97", "royalblue"), guide = "colorbar", limits = c(min_r2, max_r2)) + guides(fill = guide_colorbar("Var. (%)")) + 
  theme(axis.text.x = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 14, color = "black"), 
        axis.line = element_blank(), axis.ticks = element_blank(), 
        panel.background = element_blank(), strip.background = element_blank(), 
        strip.text = element_text(size = rel(1)))
dev.off()



# MOFA  heamap 

corr_cor <- readRDS("data/corr_mofa_factors.rds")[[["Lung"]]

heatmap_object <- Heatmap(as.matrix(corr_cor),
                          name = "Correlation",   # Name of the heatmap
                          col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),  # Color scheme
                          show_column_names = TRUE,  # Show column names
                          show_row_names = TRUE,     # Show row names
                          column_title = str_to_title(tissue),  # Column title
                          row_title = "",     # Row title
                          cluster_rows = FALSE,        # Do not cluster rows
                          cluster_columns = FALSE,     # Do not cluster columns
                          heatmap_legend_param = list(title = "Correlation", at = seq(-1, 1, 0.2)),
                          cell_fun = function(j, i, x, y, width, height, fill) {
                            grid.text(ifelse(corr_fdr[i,j] < 0.05,
                                             prettyNum(corr_cor[i, j], digits = 2),
                                             ""), x, y, gp = gpar(fontsize = 10))
                          }
)

pdf(paste0(base_image_folder, "/", tissue, "_heatmap_correlation_with_covariates.pdf"), h  = 4, w = 6)
draw(heatmap_object, heatmap_legend_side = "right")
dev.off()




# GO results

# Load ORSUM results
orsum <- read.csv("MOFA/ORSUM/all/filteredResult-Summary.tsv", sep = "\t")

tissue_order <- c("Gene Expression ", "Promoters ","Enhancers ", "Gene Body ")

colnames(orsum)[6:13] <- c("Enhancers (positive)", "Gene Body (positive)", "Gene Expression (positive)", "Promoters (positive)",  
                           "Enhancers (negative)", "Gene Body (negative)", "Gene Expression (negative)", "Promoters (negative)")

analysis <- colnames(orsum)[6:13]

enrichement_analysis_orsum <- orsum %>% 
  mutate(Representing.term.rank = as.numeric(Representing.term.rank)) %>%
  filter(Representing.term.rank <= 5) %>%
  pivot_longer(cols = all_of(analysis)) %>% 
  filter(value != "None") %>% 
  mutate(direction = gsub(".*\\((positive|negative)\\).*", "\\1", name)) %>%
  mutate(data_layer = gsub("\\((positive|negative)\\)", " ", name)) %>% 
  mutate(data_layer = gsub(" $", "", data_layer)) %>%
  mutate(`Representing.term.name` = factor(`Representing.term.name`, levels = unique(`Representing.term.name`[nrow(.):1])))


duplicated_table <- enrichement_analysis_orsum %>% 
  group_by(Representing.term.name, data_layer) %>%
  filter(n() > 1) %>%
  mutate(direction = "both") %>%
  distinct(Representing.term.name, data_layer,  .keep_all = TRUE)

enrichement_analysis_orsum <- enrichement_analysis_orsum %>% 
  group_by(Representing.term.name, data_layer) %>%
  filter(n() == 1) %>%
  bind_rows(duplicated_table)

# Organize the data
enrichement_analysis_orsum <- enrichement_analysis_orsum %>%  
  mutate(Representing.term.name = as.character(Representing.term.name)) %>%
  mutate(direction = factor(direction, levels = c("positive", "negative", "both"))) %>%
  mutate(data_layer = factor(data_layer, levels = tissue_order)) %>% 
  arrange(data_layer, direction)

enrichement_analysis_orsum$Representing.term.name <- factor(enrichement_analysis_orsum$Representing.term.name, levels = rev(unique(enrichement_analysis_orsum$Representing.term.name)))  

pdf("MOFA/enrichement_plot_permutations.pdf", w = 10, h = 10)
ggplot(enrichement_analysis_orsum, aes(x = data_layer, y = Representing.term.name, colour = direction)) +
  geom_point(size = 4) +
  scale_color_manual(name = "Gene Set", values = safe_colorblind_palette[c(2,1,13)]) +
  theme_dose(10) +
  xlab("") + ylab("") + 
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, angle = 80, hjust = 1, vjust = 1), 
        axis.text.y = element_text(size = 16), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 16), 
        legend.position = "right", 
        legend.direction="vertical")
dev.off()


#TFBS analysis is in figure_S6.R
