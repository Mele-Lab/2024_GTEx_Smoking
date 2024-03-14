# @Author: Rogerio Ribeiro and Jose Miguel Ramirez
# @E-mail: rogerio.e.ramos.ribeiro@gmail.com and jose.ramirez1@bsc.es
# @Date:   2022-03-24
# @Description: Generate publication Figure 6 (Gene expression reversability and Machine learning results)


#Set path 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Load libraries for ploting 
library(tidyverse)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(dplyr)

#Load data 
degs <- readRDS(file = "data/degs_summary.rds")

ex_smoker_data <- readRDS(file = "data/reversible_genes.rds")
for(i in 1:length(ex_smoker_data)) assign(names(ex_smoker_data)[i], ex_smoker_data[[i]])

ml_data <- readRDS(file = "data/machine_learning_gene_expression.rds")
for(i in 1:length(ml_data)) assign(names(ml_data)[i], ml_data[[i]])

ml_data_meth <- readRDS(file = "data/methylation_ml.rds")
for(i in 1:length(ml_data_meth)) assign(names(ml_data)[i], ml_data[[i]])

reverse <- readRDS(file = "data/reversible_genes.rds")

# Colours pallete
source("safe_colourblind_pallete.R")

#Load tissues data
tissue.data <- read.csv(file = "tissue_abreviation.txt")


if(!dir.exists("figures/figure_7")){
  dir.create("figures/figure_7")
}


#Figure 6A heatmap ----
reversability_per_tissue <- reversability %>% 
  group_by(tissue, reversability) %>% 
  summarise(n_genes = n()) %>% 
  pivot_wider(names_from = "reversability", values_from = "n_genes") %>% 
  mutate(total = sum(`partially reversible`, reversible, `non-reversible`, na.rm = T)) %>% 
  merge(tissue.data, by= "tissue") %>% 
  arrange(-total)

reversability_per_tissue[is.na(reversability_per_tissue)] <- 0

colnames(reversability_per_tissue)[c(2,3,4)] <- c("Partially Reversible", "Reversible", "Non Reversible")

colours <- gray.colors(200, start = 0, end = 1)[c(200:70)]


pdf("figures/figure_7/heatmap_reverse_genes.pdf",  w = 7 , h = 14)
Heatmap(as.matrix(reversability_per_tissue[,c(3,2,4)] + 1),
        col = brewer.pal(9,"BuPu")[1:7],
        width = 3*unit(15, "mm"), 
        height = ncol(reversability_per_tissue)*unit(15, "mm"), 
        na_col = "white",
        cluster_rows = F,
        cluster_columns = F,
        name = "Nº genes)",
        row_labels  = reversability_per_tissue$TISSUENAMEABREV,
        row_names_side = "left",
        column_names_side = "bottom",
        row_names_gp = gpar(fontsize = 12),
        column_names_gp = gpar(fontsize = 13),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(reversability_per_tissue[,c(3,2,4)][i, j], x, y, gp = gpar(fontsize = 13))})
dev.off()

# Figure 6B----
meth_reverse_data <- readRDS(file = "data/figure_S8_logFC_reversible.rds")$data
colnames(meth_reverse_data) <- c("Reversible", "Partially reversible", "Non reversible")

row.names(meth_reverse_data)[2] <-  "Colon transverse"


pdf("figures/figure_7/heatmap_reverse_loci.pdf",  w = 5, h = 3)
Heatmap(as.matrix(meth_reverse_data),
        col = brewer.pal(9,"BuPu")[1:7],
        width = 3*unit(15, "mm"), 
        height =  unit(12, "mm"),  
        na_col = "white",
        cluster_rows = F,
        cluster_columns = F,
        name = "Nº genes",
        row_labels  = row.names(meth_reverse_data),
        row_names_side = "left",
        column_names_side = "bottom",
        row_names_gp = gpar(fontsize = 12),
        column_names_gp = gpar(fontsize = 13),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(meth_reverse_data[i, j], x, y, gp = gpar(fontsize = 13))})
dev.off()



#Figure 6 C-E ----
#Takes as input the exp of Raw TPM's in log 
plot_genes_TPMs <- function(expression_data, title){
  plot <- expression_data %>% 
    mutate(label = factor(label, levels = c("Never Smoker", "Ex Smoker", "Smoker"))) %>%
    mutate(title = title) %>%
    ggplot( aes(x=label, y=exp)) +
    geom_violin(aes(fill = label), col = "black") +
    geom_boxplot(col = "black", outlier.shape = NA, notch = T, width = 0.25) +
    geom_jitter(col = "black", alpha = 0.1, size = 1.2) +
    theme_minimal() +
    scale_fill_manual(values=c(safe_colorblind_palette[c(1,3,2)])) +
    scale_x_discrete(labels=c("Never\nsmokers", "Ex\nsmokers", "smokers")) + ylab("Expression") + xlab("") + 
    facet_grid(. ~ title) + # Add the title strip
    stat_summary(aes(x = label, y = exp),
                 fun.data = get_box_stats, geom = "text",
                 hjust = 0.5, vjust = 0.9, size = 6) + 
    theme(axis.title.y = element_text(margin = margin(r = 20), size = 18),
          axis.text.y = element_text(size = 16, colour = "black"), 
          axis.text.x = element_text(size = 16, colour = "black"),
          panel.background = element_rect(fill = "white"),
          axis.line = element_line(colour = "black"),
          legend.position = "none", 
          plot.title = element_text(hjust = 0.5, size = 24),
          strip.background = element_rect(fill="#88CCEE"), 
          strip.text = element_text(size = 20)) + 
    xlab("") + 
    ylab("Log(TPM)")
  
  print(plot)
}

#Set function to help ploting the boxplots
get_box_stats <- function(y, upper_limit = max(y) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "N =", length(y), "\n"
    )
  ))
}



pdf("figures/figure_7/example_reversible_genes.pdf", w = 6, h = 6)
plot_genes_TPMs(reversible_example, "CYP1A1 (Lung)")
dev.off()

pdf("figures/figure_7/example_pr_genes.pdf", w = 6, h = 6)
plot_genes_TPMs(pr_gene_example, "GPR15 (Lung)")
dev.off()

pdf("figures/figure_7/example_nr_genes.pdf", w = 6, h = 6)
plot_genes_TPMs(non_reversible_genes_example, "TMEM59L (Thyroid)")
dev.off()


#Figure 6F ----
logFC <- reverse$logFC
sig <-  wilcox.test(abs(logFC$logFC_NS_EX), abs(logFC$logFC_EX_VS_S), paired = T)$p.value


# logFC %>%
#   dplyr::select(logFC_NS_EX,logFC_EX_VS_S, tissue) %>%
#   pivot_longer(!tissue, names_to = "comparison", values_to = "logFC") %>%
#   ggplot(aes(x = abs(logFC), color  = comparison)) +
#   stat_ecdf(linewidth = 1.5) +
#   #geom_line() +
#   scale_color_manual(labels = c("Ex-Smokers vs Smokers", "Never Smokers vs Ex-Smokers"), values = safe_colorblind_palette[c(2,1)], name = "Smoking") +
#   theme_classic() +
#   theme(
#     legend.position="top",
#     legend.title = element_blank(),
#     legend.text = element_text(size = 14),
#     axis.text=element_text(size=15),
#     axis.title.y = element_text(margin = margin(r = 20), size = 18),
#     axis.text.x = element_text(size = 16, face = "bold"),
#     axis.text.y = element_text(size = 18, face = "bold"),
#     axis.title = element_text(size = 18),
#     panel.background = element_rect(fill = "white"),
#     axis.line = element_line(colour = "grey"),
#   )


logFC_ns_ex <- logFC %>% select(logFC_NS_EX) %>% arrange(abs(logFC_NS_EX)) %>% mutate(rank = 1:nrow(.))
logFC_es_s <- logFC %>% select(logFC_EX_VS_S) %>% arrange(abs(logFC_EX_VS_S)) %>% mutate(rank = 1:nrow(.))

logFC_ranked <- merge(logFC_ns_ex, logFC_es_s, by = "rank")

pdf("figures/figure_7//logFC_cumulative_plot.pdf", w = 10, h = 10)
logFC_ranked %>% 
  dplyr::select(logFC_NS_EX,logFC_EX_VS_S, rank) %>% 
  ggplot() + 
  geom_line(aes(y = cumsum(abs(logFC_NS_EX)), x = rank, color = 'Never Smokers\nvs Ex Smokers'), size = 1.5) + 
  geom_line(aes(y = cumsum(abs(logFC_EX_VS_S)), x = rank, color = 'Ex Smokers\nvs Smokers'), size = 1.5) + 
  theme_classic() +
  scale_color_manual(name='Comparison',
                     breaks=c('Never Smokers\nvs Ex Smokers', 'Ex Smokers\nvs Smokers'),
                     values=c('Never Smokers\nvs Ex Smokers'= safe_colorblind_palette[1],  'Ex Smokers\nvs Smokers'=  safe_colorblind_palette[2])) + 
  theme(
    legend.position="none",
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    axis.text=element_text(size=15),
    axis.title.y = element_text(margin = margin(r = 20), size = 30),
    axis.text.x = element_text(size = 26, colour = "black"),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title = element_text(size = 18), 
    panel.background = element_rect(fill = "white"),
    axis.line = element_line(colour = "black"),
    axis.ticks.x.top = element_blank(), 
    axis.text.x.top = element_blank(),
    axis.ticks.y.right = element_blank(), 
    axis.text.y.right = element_blank(), 
  )+ 
  scale_x_continuous(sec.axis = sec_axis( ~ . * 1, name = " ")) +  #Used to make the box
  scale_y_continuous(sec.axis = sec_axis( ~ . * 1, name = " ")) + #Used to make the box
  ylab("Cumulative LogFC") 
dev.off()



# Figure 6g ----

figure_data_ml <- readRDS(file = "data/figure_S8_logFC_reversible.rds") #list of expression_barplot, expression_pirate, splicing_barplot and splicing_pirate

logFC_methyl <- figure_data_ml$logFC %>% as.data.frame()
p_val <-  wilcox.test(abs(logFC_methyl$logFC_never_ex), abs(logFC_methyl$logFC_ex_smokers), paired = T)$p.value

#p_val <- test$p.value
#p_val <- formatC(p_val, format = "e", digits = 2)

never_ex <- data.frame(logFC_methyl[,1])
colnames(never_ex) <- "logFC"
never_ex <- never_ex %>% 
  arrange(logFC )


ex_smk <- data.frame(logFC_methyl[,2])
colnames(ex_smk) <- "logFC"
ex_smk <- ex_smk %>% 
  arrange(logFC )

logFC_methyl_2 <- cbind(never_ex, ex_smk)
logFC_methyl_2$rank <- 1:nrow(logFC_methyl_2)
colnames(logFC_methyl_2) <- c("logFC_never_vs_ex", "logFC_ex_vs_smk", "rank")


pdf("figures/figure_7//logFC_cumulative_plot_ml.pdf", w = 10, h = 10)
logFC_methyl_2 %>%
  ggplot() + 
  geom_line(aes(y = cumsum(abs(logFC_never_vs_ex)), x = rank, color = 'Never Smokers\nvs Ex Smokers'), size = 1.5) + 
  geom_line(aes(y = cumsum(abs(logFC_ex_vs_smk)), x = rank, color = 'Ex Smokers\nvs Smokers'), size = 1.5) + 
  theme_classic() +
  scale_color_manual(name='Comparison',
                     breaks=c('Never Smokers\nvs Ex Smokers', 'Ex Smokers\nvs Smokers'),
                     values=c('Never Smokers\nvs Ex Smokers'= safe_colorblind_palette[1],  'Ex Smokers\nvs Smokers'=  safe_colorblind_palette[2])) + 
  theme(
    legend.position="none",
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    axis.text=element_text(size=15),
    axis.title.y = element_text(margin = margin(r = 20), size = 30),
    axis.text.x = element_text(size = 26, colour = "black"),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title = element_text(size = 18), 
    panel.background = element_rect(fill = "white"),
    axis.line = element_line(colour = "black"),
    axis.ticks.x.top = element_blank(), 
    axis.text.x.top = element_blank(),
    axis.ticks.y.right = element_blank(), 
    axis.text.y.right = element_blank(), 
  )+ 
  scale_x_continuous(sec.axis = sec_axis( ~ . * 1, name = " ")) +  #Used to make the box
  scale_y_continuous(sec.axis = sec_axis( ~ . * 1, name = " ")) + #Used to make the box
  ylab("Cumulative LogFC") + 
  xlab("Number of positions")
dev.off()

# ChrommHMM figure 


results <- readRDS("data/reversibile_chromHMM.rds")


g1 <- ggplot(results, aes(x=log(OR), y=region, colour=subset, alpha=sig)) + 
  geom_errorbar(aes(xmin=log(CI_lower), xmax=log(CI_upper)), width=.4) + 
  geom_vline(xintercept = 0) + 
  geom_point(size=2.5) + ylab('') + theme_bw() +
  scale_colour_manual(name="CpGs closer to", 
                      values=c("#88CCEE", "#CC6677"), 
                      labels = c("Never smokers", "Smokers")) +
  xlab("Log(Odds ratio)") +
  scale_alpha_discrete(range = c(0.4, 1), drop = FALSE) +
  theme(axis.text.x = element_text(colour="black", size=11),
        axis.text.y = element_text(colour="black", size=13),
        legend.text = element_text(colour="black", size=12),
        axis.title.x = element_text(size=13),
        legend.spacing.y = unit(-0.05, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth=1))

g2 <- ggplot(results) + geom_col(aes(sample_size, region, fill=subset), width = 0.6) +
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
  scale_x_continuous(breaks=c(0, 10000, 20000)) #Only for lung

pdf("figures/figure_7/chromHMM.pdf", w = 7.5, h = 6)
ggarrange(g1, g2, ncol=2, common.legend = TRUE, legend="right", widths = c(1, 0.3))
dev.off()


# Figure 6I----
ml_prediction_meth <- ml_data_meth$ex_smoker_predictions %>% 
  group_by(ex_smoker_pred) %>% 
  summarise(n_ex = n())

ml_prediction_meth$percentage <- ml_prediction_meth$n_ex / sum(ml_prediction_meth$n_ex)
ml_prediction_meth$tissue <- c("Lung (40)")
ml_prediction_meth$total_n_ex <- sum(ml_prediction_meth$n_ex)
ml_prediction_meth$pred <- ml_prediction_meth$ex_smoker_pred
ml_prediction_meth <- ml_prediction_meth %>%  select("tissue", "n_ex", "pred")

ex_smoker_pred <- ex_smoker_pred %>% select("tissue", "n_ex", "pred")
ex_smoker_pred$pred <- ifelse(ex_smoker_pred$pred == "Never Smoker", "Never smoker", ml_prediction_meth$pred)

ex_smoker_pred$type <- "GeneExpression"
ml_prediction_meth$type <- "Methylation"

ex_smoker_pred <- rbind(ex_smoker_pred, ml_prediction_meth) 

#Histology figure ex-smokers
hist_ex_smokers_summarised <- readRDS(file = "data/ex_smoker_hist.rds")

total_number_of_ex_smokers <- hist_ex_smokers_summarised %>% 
  group_by(tissue) %>%
  summarise(n_total = sum(n))


total_number_of_ex_smokers[1,1] <- "Esophagus\nmucosa"


hist_ex_smokers_summarised$tissue <- ifelse(hist_ex_smokers_summarised$tissue == "Esophagus...Mucosa",  "Esophagus\nmucosa", hist_ex_smokers_summarised$tissue )
hist_ex_smokers_summarised$prediction <- ifelse(hist_ex_smokers_summarised$prediction == "Never-Smoker",  "Never Smoker", "Smoker")

hist_ex_smokers_summarised <- merge(hist_ex_smokers_summarised, total_number_of_ex_smokers, by = "tissue")

hist_ex_smokers_summarised$tissue <- paste0(hist_ex_smokers_summarised$tissue, "(", hist_ex_smokers_summarised$n_total, ")")

hist_ex_smokers_summarised$n_ex <- hist_ex_smokers_summarised$n
hist_ex_smokers_summarised$pred <- hist_ex_smokers_summarised$prediction

hist_ex_smokers_summarised <- hist_ex_smokers_summarised %>% select(tissue, n_ex, pred)
ex_smoker_pred <- ex_smoker_pred %>% select(tissue, n_ex, pred, type)

hist_ex_smokers_summarised$type <- "Hist"

#Reoder by smoker percent
hist_ex_smokers_summarised <- hist_ex_smokers_summarised[c(3:4, 7:8, 1:2, 5:6), ]

# Merge all the datasets
ex_smoker_pred <- rbind(ex_smoker_pred,
                        hist_ex_smokers_summarised)

# Put methylation at the end
ex_smoker_pred <- ex_smoker_pred[c(1:22, 25:32, 23,24), ]

ex_smoker_pred$tissue2 <- paste0(ex_smoker_pred$tissue, ex_smoker_pred$type)
ex_smoker_pred$tissue2 <- factor(ex_smoker_pred$tissue2, levels = unique(ex_smoker_pred$tissue2)[length(unique(ex_smoker_pred$tissue2)):1])

number_of_ex_smokers <- gsub("\\(|\\)", "", str_extract(ex_smoker_pred$tissue, "\\((\\d+)\\)"))
number_of_ex_smokers <- data.frame(n = as.numeric(number_of_ex_smokers), tissue = ex_smoker_pred$tissue2) %>% distinct(tissue, .keep_all = TRUE)

pdf("figures/figure_7/ex_smoker_classification_percentage.pdf", w = 9, h =12)
ex_smoker_pred %>% 
  ggplot(aes(y = tissue2, x  = n_ex, fill = pred)) + 
  geom_bar(position = "fill", stat = "identity", width = 0.70) + 
  theme_classic() + 
  ylab("") + 
  xlab("Percentage ex-smokers (%)") +
  scale_fill_manual(values = safe_colorblind_palette[c(1,2)]) +
  #scale_y_continuous(expand = c(0,0)) + 
  geom_vline(xintercept=0.5) + 
  theme(legend.position =  "right", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 16, colour = "black"),
        axis.title.y = element_text(size = 18, colour = "black", margin = margin(t = 20)), 
        axis.text.x = element_text(size = 14, colour = "black"), 
        axis.text.y = element_text(size = 16, colour = "black", hjust = 1, vjust = 0.3), 
        axis.title.x = element_text(size = 17, colour = "black", margin = margin(t = 0, r = 20, b = 0, l = 0)))
dev.off()  


pdf("figures/figure_7/ex_smoker_number.pdf", w = 6.5, h =12)
ggplot(number_of_ex_smokers, aes(y = tissue, x  = n)) + 
  geom_bar(stat = "identity", width = 0.70, fill = safe_colorblind_palette[3]) + 
  theme_classic() + 
  ylab("") + 
  xlab("Number of ex-smokers") +
  #scale_y_continuous(expand = c(0,0)) + 
  geom_vline(xintercept=0.5) + 
  theme(legend.position =  "right", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 16, colour = "black"),
        axis.title.y = element_text(size = 18, colour = "black", margin = margin(t = 20)), 
        axis.text.x = element_text(size = 14, colour = "black"), 
        axis.text.y = element_text(size = 16, colour = "black", hjust = 1, vjust = 0.3), 
        axis.title.x = element_text(size = 17, colour = "black", margin = margin(t = 0, r = 20, b = 0, l = 0)))
dev.off()
