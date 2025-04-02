# @Author: Rogerio Ribeiro and Jose Miguel Ramirez
# @E-mail: rogerio.e.ramos.ribeiro@gmail.com and jose.ramirez1@bsc.es
# @Date:   2023-03-13
# @Description: Generate publication Figure s8


#Set path 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Load libs
library(tidyverse)
library(ggrepel)
library(ggpubr)


# Load data 
source("safe_colourblind_pallete.R")
tissue.data <- read.csv(file = "tissue_abreviation.txt")


## Load data

res_smoking_overlap <- readRDS(file = "data/smoking_age_causal_cpg_overlaps.rds")
res_smoking_overlap_dir <- readRDS(file = "data/smoking_overlap_causal_direction.rds")

dir <- "figures/figure_s9//"
#dir.create(dir)




# Plot the overlap in smoking genes
size <- 15
res_smoking_overlap$trait <- str_trim(res_smoking_overlap$trait)
res_smoking_overlap <- merge(res_smoking_overlap, tissue.data, by = "tissue")
res_smoking_overlap$trait <-  gsub("_", " ", res_smoking_overlap$trait)
res_smoking_overlap$sig <- factor(ifelse(res_smoking_overlap$FDR > 0.05, "FDR > 0.05", "FDR < 0.05"), levels = c("FDR > 0.05", "FDR < 0.05"))



g <- ggplot(res_smoking_overlap, aes(x=log(OR), y=trait, alpha=sig)) + 
  geom_errorbar(aes(xmin=log(CL), xmax=log(CH)), width=.3, colour = "black") + 
  geom_vline(xintercept = 0) +
  geom_point(size=3, colour = "black") + ylab('') + theme_bw() +
  xlab("Log(Odds ratio)") +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour="black", size=size),
        axis.text.y = element_text(colour="black", size=size+1),
        axis.title.x = element_text(size=size+1),
        legend.text = element_text(size=size+1),
        legend.key.size = unit(0.9,"cm"), 
        strip.background = element_rect(fill = "#88CCEE"),
        strip.text = element_text(size = 16)) +
  guides(color = guide_legend(override.aes = list(size = 3)),
         alpha = guide_legend(override.aes = list(size = 3))) +
  facet_wrap(~ Name,  strip.position = "top") +
  xlim(-4, 4)

pdf(paste0(dir, "smoking_overlap.pdf"), w = 10)
g 
dev.off()



## Plot overlap of hyper and hypo positions with causal and dmg
wrangle_data <- function(data){
  data <- data %>% 
    filter(tissue == "Lung") %>%  #Keep only lung samples
    mutate(FDR = p.adjust(pvalue, method = "BH")) %>%
    merge(tissue.data, by = "tissue")
}


hyper_dmg_lung <- wrangle_data(res_smoking_overlap_dir[[1]])
hypo_dmg_lung <- wrangle_data(res_smoking_overlap_dir[[3]])
hyper_protective_lung <- wrangle_data(res_smoking_overlap_dir[[2]])
hypo_protective_lung <- wrangle_data(res_smoking_overlap_dir[[4]])


dmg_lung <- rbind(hyper_dmg_lung %>% mutate(dir = "Hypermethylation"), hypo_dmg_lung %>% mutate(dir = "Hypomethylation"))
dmg_lung$sig <- factor(ifelse(dmg_lung$FDR < 0.05, "FDR < 0.05", "FDR < 0.05"), levels = c("FDR < 0.05", "FDR > 0.05"))

g2 <- ggplot(dmg_lung, aes(x=log(OR), y=trait, colour=dir, alpha=sig)) + 
  geom_errorbar(aes(xmin=log(CL), xmax=log(CH)), width=.3) + 
  geom_vline(xintercept = 0) +
  geom_point(size=3) + ylab('') + theme_bw() +
  scale_colour_manual(values=c("#CC6677", "#88CCEE")) +
  xlab("Log(Odds ratio)") +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour="black", size=size),
        axis.text.y = element_text(colour="black", size=size+1),
        axis.title.x = element_text(size=size+1),
        legend.text = element_text(size=size+1),
        legend.key.size = unit(0.9,"cm")) +
  scale_alpha_manual(values=c("FDR < 0.05" = 1, "FDR > 0.05" = 1),
                     breaks=c("FDR < 0.05", "FDR > 0.05")) +
  guides(color = guide_legend(override.aes = list(size = 3)),
         alpha = guide_legend(override.aes = list(size = 3))) +
  xlim(-4, 4) +facet_wrap(~ dir,  strip.position = "top")



pdf(paste0(dir, "smoking_overlap_dmg.pdf"), w = 11.5)
g2
dev.off()


protc_lung <- rbind(hyper_protective_lung %>% mutate(dir = "Hypermethylation"), hypo_protective_lung %>% mutate(dir = "Hypomethylation"))
protc_lung$sig <- factor(ifelse(hyper_protective_lung$FDR < 0.05, "FDR < 0.05", "FDR < 0.05"), levels = c("FDR < 0.05", "FDR > 0.05"))

g3 <- ggplot(protc_lung, aes(x=log(OR), y=trait, colour=dir, alpha=sig)) + 
  geom_errorbar(aes(xmin=log(CL), xmax=log(CH)), width=.3) + 
  geom_vline(xintercept = 0) +
  geom_point(size=3) + ylab('') + theme_bw() +
  scale_colour_manual(values=c("#CC6677", "#88CCEE")) +
  xlab("Log(Odds ratio)") +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour="black", size=size),
        axis.text.y = element_text(colour="black", size=size+1),
        axis.title.x = element_text(size=size+1),
        legend.text = element_text(size=size+1),
        legend.key.size = unit(0.9,"cm")) +
  scale_alpha_manual(values=c("FDR < 0.05" = 1, "FDR > 0.05" = 1),
                     breaks=c("FDR < 0.05", "FDR > 0.05")) +
  guides(color = guide_legend(override.aes = list(size = 3)),
         alpha = guide_legend(override.aes = list(size = 3))) + 
  xlim(-4, 4) + 
  facet_wrap(~ dir,  strip.position = "top")


pdf(paste0(dir, "smoking_overlap_prot.pdf"), w = 11.5)
g3 
dev.off()


pdf(paste0(dir, "plots_no_labels.pdf"), w = 7, h = 12)
ggarrange(g+theme(legend.position = "None"), g2+theme(legend.position = "None"), g3+theme(legend.position = "None"), nrow = 3)
dev.off()

# Methylation plot
tissue_convert <- tissue.data$TISSUENAMEABREV
names(tissue_convert) <- tissue.data$SMTSDNOSEP

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")


# Boxplot aux function

get_box_stats <- function(y, upper_limit = max(y) * 1.25) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "N =", length(y), "\n"
    )
  ))
}

size <- 15

theme_agingclock <- theme(legend.title = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          axis.text.x = element_text(colour="black", size=size),
                          axis.text.y = element_text(colour="black", size=size+1),
                          axis.title.x = element_text(size=size+1),
                          legend.text = element_text(size=size+1),
                          legend.key.size = unit(0.9,"cm"), 
                          strip.background = element_rect(fill = "#88CCEE"),
                          strip.text = element_text(size = 13)) +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) 

## Load the clock data

## Hovarth clock 
data_ht <- readRDS("scripts/methylation//hovarth_molecular_clocks.rds")
#  DamAge and AdaptAge
vadimAge <- readRDS("scripts/methylation/vadims_molecular_clocks.rds")

# Figure S8a -----
# Hovarth Clocks
data_ht$comparison$sig <- factor(ifelse(data_ht$comparison$FDR > 0.05, "FDR > 0.05", "FDR < 0.05"), levels = c("FDR > 0.05", "FDR < 0.05"))
data_ht$comparison$estimate <- as.numeric(data_ht$comparison$estimate)
data_ht$comparison$lCI <- as.numeric(data_ht$comparison$lCI)
data_ht$comparison$hCI <- as.numeric(data_ht$comparison$hCI)

data_ht$comparison$tissue <- tissue_convert[data_ht$comparison$tissue]
data_ht$comparison$tissue <- factor(data_ht$comparison$tissue, levels = unique(data_ht$comparison$tissue)[9:1])
data_ht$comparison$comparison <- factor(data_ht$comparison$comparison, levels = c("Never vs Smoking","Never vs Ex-Smoker", "Ex-Smoker vs Smoker"))

g1 <- ggplot(data_ht$comparison, aes(x=estimate, y=tissue, alpha=sig)) + 
  geom_errorbar(aes(xmin=lCI, xmax=hCI), width=.3, colour = "black") + 
  geom_vline(xintercept = 0) +
  geom_point(size=3, colour = "black") + ylab('') + theme_bw() +
  xlab("Coefficient (Years)") +
  ggtitle("Horvath clock") + 
  theme_agingclock + 
  guides(color = guide_legend(override.aes = list(size = 3)),
         alpha = guide_legend(override.aes = list(size = 3))) +
  facet_wrap(~ comparison,  strip.position = "top")




# Vadim Clocks
vadimAge$adapt$comparison$sig <- factor(ifelse(vadimAge$adapt$comparison$FDR > 0.05, "FDR > 0.05", "FDR < 0.05"), levels = c("FDR > 0.05", "FDR < 0.05"))
vadimAge$adapt$comparison$estimate <- as.numeric(vadimAge$adapt$comparison$estimate)
vadimAge$adapt$comparison$lCI <- as.numeric(vadimAge$adapt$comparison$lCI)
vadimAge$adapt$comparison$hCI <- as.numeric(vadimAge$adapt$comparison$hCI)

vadimAge$adapt$comparison$tissue <- tissue_convert[vadimAge$adapt$comparison$tissue]
vadimAge$adapt$comparison$tissue <- factor(vadimAge$adapt$comparison$tissue, levels = unique(vadimAge$adapt$comparison$tissue)[9:1])
vadimAge$adapt$comparison$comparison <- factor(vadimAge$adapt$comparison$comparison, levels = c("Never vs Smoking","Never vs Ex-Smoker", "Ex-Smoker vs Smoker"))

g2 <- ggplot(vadimAge$adapt$comparison, aes(x=estimate, y=tissue, alpha=sig)) + 
  geom_errorbar(aes(xmin=lCI, xmax=hCI), width=.3, colour = "black") + 
  geom_vline(xintercept = 0) +
  geom_point(size=3, colour = "black") + ylab('') + theme_bw() +
  xlab("Coefficient (Years)") +
  ggtitle("AdaptAge clock") + 
  theme_agingclock + 
  guides(color = guide_legend(override.aes = list(size = 3)),
         alpha = guide_legend(override.aes = list(size = 3))) +
  facet_wrap(~ comparison,  strip.position = "top")


vadimAge$damn$comparison$sig <- factor(ifelse(vadimAge$damn$comparison$FDR > 0.05, "FDR > 0.05", "FDR < 0.05"), levels = c("FDR > 0.05", "FDR < 0.05"))
vadimAge$damn$comparison$estimate <- as.numeric(vadimAge$damn$comparison$estimate)
vadimAge$damn$comparison$lCI <- as.numeric(vadimAge$damn$comparison$lCI)
vadimAge$damn$comparison$hCI <- as.numeric(vadimAge$damn$comparison$hCI)

vadimAge$damn$comparison$tissue <- tissue_convert[vadimAge$damn$comparison$tissue]
vadimAge$damn$comparison$tissue <- factor(vadimAge$damn$comparison$tissue, levels = unique(vadimAge$damn$comparison$tissue)[9:1])
vadimAge$damn$comparison$comparison <- factor(vadimAge$damn$comparison$comparison, levels = c("Never vs Smoking","Never vs Ex-Smoker", "Ex-Smoker vs Smoker"))

g3 <- ggplot(vadimAge$damn$comparison, aes(x=estimate, y=tissue, alpha=sig)) + 
  geom_errorbar(aes(xmin=lCI, xmax=hCI), width=.3, colour = "black") + 
  geom_vline(xintercept = 0) +
  geom_point(size=3, colour = "black") + ylab('') + theme_bw() +
  xlab("Coefficient (Years)") +
  ggtitle("DamAge clock") + 
  theme_agingclock +
  guides(color = guide_legend(override.aes = list(size = 3)),
         alpha = guide_legend(override.aes = list(size = 3))) +
  facet_wrap(~ comparison,  strip.position = "top")

pdf("scripts/methylation/figure_s8.pdf", h = 9, w = 10)
ggarrange(g1,g2,g3, nrow = 3, labels = c("a", "b", "c"),
          font.label = list(size = 18, color = "black", face = "bold", family = NULL))
dev.off()