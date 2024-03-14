# @Author: Rogerio Ribeiro
# @E-mail: rogerio.e.ramos.ribeiro@gmail.com
# @Date:   2022-03-24
# @Description: Generate publication Figure 7 (Histology analysis) 



#Set path 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Load libraries for ploting 
library(tidyverse)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(pROC)
library(plotROC)

tissue.data <- read.csv(file = "tissue_abreviation.txt")


#Load colours 
source("safe_colourblind_pallete.R")

## Figure 7A
figure_7a <- read.csv(file = "data/figure_7_a_data.csv", sep = "\t")

figure_7a$smoker_status <- factor(figure_7a$smoker_status, levels = c("Smoker", "Ex-smoker", "Never Smoker"))
figure_7a$tissue <- factor(figure_7a$tissue, levels = c("Pancreas", "Esophagus mucosa", "Thyroid", "Lung")[4:1])

figure_7a <- figure_7a %>% 
  group_by(tissue, smoker_status) %>% 
  summarise(n_subjects = sum(n_subjects), n_images = sum(n_images))

pdf("figures//figure_7/figure_7A.pdf", w = 7, h = 6)
ggplot(figure_7a, aes(x = n_subjects, y = tissue, fill = smoker_status)) + 
  geom_bar(position="stack", stat="identity") + 
  coord_flip() + 
  theme_bw() + 
  scale_fill_manual(name = "",values=safe_colorblind_palette[c(2,3,1)]) +
  ylab("") + xlab("Number of Images") + 
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        axis.text.y= element_text(size = 15, colour = "black"),
        axis.ticks.y=element_blank(), 
        axis.text.x =  element_text(size = 18, colour = "black", angle = 90, vjust = 0.3, hjust = 1), 
        axis.title.y =  element_text(size = 18, colour = "black"), 
        legend.text = element_text(size = 16))

dev.off()

## Figure 7B

eval_testSet_allTissues2 <- read_csv("D:/PhD/smoking_analysis.v4/histology/eval_testSet_allTissues2.csv")

eval_all_tissues <- rbind(eval_testSet_allTissues2 %>% select( true_smoking, pred_Lung) %>% rename(pred = "pred_Lung") %>% mutate(tissue = "Lung"), 
                          rbind(eval_testSet_allTissues2 %>% select( true_smoking, pred_Thyroid)  %>% rename(pred = "pred_Thyroid") %>% mutate(tissue = "Thyroid"),
                                rbind(eval_testSet_allTissues2 %>% select( true_smoking, pred_Pancreas)  %>% rename(pred = "pred_Pancreas") %>% mutate(tissue = "Pancreas"),
                                      eval_testSet_allTissues2 %>% select( true_smoking, pred_Esophagus_Mucosa )  %>% rename(pred = "pred_Esophagus_Mucosa") %>% mutate(tissue = "Esophagus mucosa")
                                )
                          )
)


eval_all_tissues$tissue <- factor(eval_all_tissues$tissue, levels = c("Pancreas", "Esophagus mucosa", "Thyroid", "Lung")[4:1])


pdf("figures/figure_7/figure_7B.pdf", w = 9, h = 8)
ggplot(eval_all_tissues, aes(d = true_smoking, m = pred, colour = tissue)) +
  geom_roc(n.cuts = 0, size = 2) +
  labs(x = "False Positive Rate",
       y = "True Positive Rate") +
  geom_abline(slope=1, intercept=  0, size = 1, linetype  =3) + 
  scale_color_manual(values = c("Lung" = safe_colorblind_palette[5],
                                "Thyroid" = safe_colorblind_palette[6],
                                "Pancreas" = safe_colorblind_palette[8], 
                                "Esophagus mucosa" = safe_colorblind_palette[4])) +
  theme_classic() + 
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    #panel.border = element_rect(colour = "black"),
    axis.ticks.x=element_blank(), 
    axis.text =  element_text(size = 24, colour = "black"), 
    axis.title =  element_text(size = 20, colour = "black"),
    legend.position = "none"
  )
dev.off()

pdf("figures//figure_7/figure_7BLegend.pdf", w = 8, h = 8)
ggplot(eval_all_tissues, aes(d = true_smoking, m = pred, colour = tissue)) +
  geom_roc(n.cuts = 0) +
  labs(x = "False Positive Rate",
       y = "True Positive Rate") +
  scale_color_manual(values = c("Lung" = safe_colorblind_palette[5],
                                "Thyroid" = safe_colorblind_palette[6],
                                "Pancreas" = safe_colorblind_palette[8], 
                                "Esophagus mucosa" = safe_colorblind_palette[4])) +
  theme_classic() + 
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    #panel.border = element_rect(colour = "black"),
    axis.ticks.x=element_blank(), 
    axis.text =  element_text(size = 18, colour = "black"), 
    axis.title =  element_blank(),
    legend.text = element_text(size = 18), 
    legend.title = element_blank()
  )
dev.off()




## Figure 7D
bayes <- readRDS(file = "data/bayesPrism.rds")

data_macro_smk_never_smokers <- bayes$dataSig %>% 
  mutate(Smoking = ifelse(Smoking == 0, "Never smoker", ifelse(Smoking == 1, "Ex smoker", "Smokers")))

data_macro_smk_never_smokers$Smoking <- factor(data_macro_smk_never_smokers$Smoking, levels = c("Never smoker", "Ex smoker", "Smokers"))

#Set function to help ploting the boxplots
get_box_stats <- function(y, upper_limit = max(y) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "N =", length(y), "\n"
    )
  ))
}



a <- ggplot(data_macro_smk_never_smokers, aes(x=Smoking, y=value)) +
  geom_violin(aes(fill = Smoking), col = "black") +
  geom_boxplot(col = "black", outlier.shape = NA, notch = T, width = 0.25) +
  geom_jitter(col = "black", alpha = 0.1, size = 2) +
  scale_fill_manual(values=c(safe_colorblind_palette[c(1,3,2)])) +
  ylab("Macrophages proportions") + 
  stat_summary(aes(x = Smoking, y = value),
               fun.data = get_box_stats, geom = "text",
               hjust = 0.5, vjust = 0.9, size = 12) + 
  ggtitle("") +
  theme_minimal() +
  theme(axis.title.y = element_text(margin = margin(r = 20), size = 24),
        axis.text.y = element_text(size = 24, colour = "black"), 
        axis.text.x = element_text(size = 28, colour = "black"),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"),
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 24),
        strip.background = element_rect(fill="#88CCEE"), 
        strip.text = element_text(size = 20)) + 
  xlab("") + 
  ylab("Macrophages Proportions") + 
  scale_x_discrete(labels=c("Never smokers", "Ex-smoker","Smokers"))

pdf("figures/figure_7/figure_7d.pdf", w = 8, h = 8.5)
a
dev.off()


## Figure 7f ---


medianFollicleDiameter <- read.csv(file = "data/medianFollicleDiameter.csv")

b <- ggplot(data = medianFollicleDiameter, aes(x = smoking, y =MedianFollicleDiameter)) + 
  geom_violin(aes(fill = smoking),
              col = "black") +
  geom_boxplot(col = "black",
               outlier.shape = NA,
               notch = T,
               width = 0.25) +
  geom_jitter(col = "black",
              alpha = 0.1,
              size = 0.8) +
  # stat_summary( #to add median
  #   geom = "point",
  #   fun = "median",
  #   col = "black",
  #   size = 3,
  #   shape = 24,
  #   fill = "red"
  # ) +
  theme_minimal() +
  scale_fill_manual(values=c("#88CCEE", "#CC6677")) +
  ylab("median(Follicle diameter) (px)") +  ggtitle("") +
  stat_summary(aes(x = smoking, y = MedianFollicleDiameter),
               fun.data = get_box_stats, geom = "text",
               hjust = 0.5, vjust = 0.9, size = 12) + 
  theme(axis.title.y = element_text(margin = margin(r = 20), size = 24),
        axis.text.y = element_text(size = 24, colour = "black"), 
        axis.text.x = element_text(size = 28, colour = "black"),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"),
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 24),
        strip.background = element_rect(fill="#88CCEE"), 
        strip.text = element_text(size = 20)) +
  scale_x_discrete(labels=c("Never smokers", "Smokers")) + xlab("")

## Figure 7g ---

sdFollicleDiameter <- read.csv(file = "data/sdFollicleDiameter.csv")

c <- ggplot(data = sdFollicleDiameter, aes(x = smoking, y =SdFollicleDiameter)) + 
  geom_violin(aes(fill = smoking),
              col = "black") +
  geom_boxplot(col = "black",
               outlier.shape = NA,
               notch = T,
               width = 0.25) +
  geom_jitter(col = "black",
              alpha = 0.1,
              size = 0.8) +
  # stat_summary( #to add mean
  #   geom = "point",
  #   fun = "mean",
  #   col = "black",
  #   size = 3,
  #   shape = 24,
  #   fill = "red"
  # ) +
  theme_minimal() +
  scale_fill_manual(values=c("#88CCEE", "#CC6677")) +
  ylab("sd(Follicle diameter) (px)") +  ggtitle("") +
  stat_summary(aes(x = smoking, y = SdFollicleDiameter),
               fun.data = get_box_stats, geom = "text",
               hjust = 0.5, vjust = 0.9, size = 12) + 
  theme(axis.title.y = element_text(margin = margin(r = 20), size = 24),
        axis.text.y = element_text(size = 24, colour = "black"), 
        axis.text.x = element_text(size = 28, colour = "black"),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"),
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 24),
        strip.background = element_rect(fill="#88CCEE"), 
        strip.text = element_text(size = 20)) +
  scale_x_discrete(labels=c("Never smokers", "Smokers")) + xlab("")


pdf("figures/figure_7/figure_7fg.pdf", w = 12, h = 8.5)
ggpubr::ggarrange(b, c, ncol = 2)
dev.off()