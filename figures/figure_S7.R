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

#TFBS analysis is in figure_S6.R