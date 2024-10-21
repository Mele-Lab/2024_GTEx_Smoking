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
