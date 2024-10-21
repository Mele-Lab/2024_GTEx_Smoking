# @Author: Rogerio Ribeiro and Jose Miguel Ramirez
# @E-mail: rogerio.e.ramos.ribeiro@gmail.com and jose.ramirez1@bsc.es
# @Date:   2024-09-06
# @Description: Generate publication Figure S8 (Bioclocks)

#Set working dir 
#setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load libraries
library(tidyverse)
library(ggpubr)

# Create output dir
dir <- "figures/figure_s8"
if (!dir.exists("figures/figure_s8")){
  dir.create("figures/figure_s8")
}

source("safe_colourblind_pallete.R")

# Boxplot aux function

get_box_stats <- function(y, upper_limit = max(y) * 1.25) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "N =", length(y), "\n"
    )
  ))
}


## Load the clock data

## Hovarth clock 
data_ht <- readRDS("data/figure_clock_hovath.rds")
#  DamAge and AdaptAge
vadimAge <- readRDS("data/figure_clock_vadim.rds")


# Figure S8a -----
# Get data 
data_lung <- data_ht$Lung_residuals_regression
data_lung$smoking_status <- factor(ifelse(data_lung$smoking_status == 0, "Never\nsmoker", "Smoker"))


a <- ggplot(data_lung, aes(x=smoking_status, y=error)) +
  geom_violin(aes(fill = smoking_status), col = "black") +
  geom_boxplot(col = "black", outlier.shape = NA, notch = T, width = 0.25) +
  geom_jitter(col = "black", alpha = 0.1, size = 1.2) +
  theme_minimal() +
  scale_fill_manual(values=c(safe_colorblind_palette[c(1,2)])) +
  scale_x_discrete(labels=c("Never\nsmokers", "Smokers")) + 
  ylab("Expression") + xlab("") + 
  facet_grid(. ~ "Horvath clock") + # Add the title strip
  stat_summary(aes(x = smoking_status, y = error),
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
  ylab("Adjusted Error")



# Figure S8b -----
# Get data 
data_lung_adapt <- vadimAge$Lung_residuals_adapt
data_lung_adapt$smoking_status <- factor(ifelse(data_lung_adapt$smoking_status == 0, "Never\nsmoker", "Smoker"))


b <- ggplot(data_lung_adapt, aes(x=smoking_status, y=error)) +
  geom_violin(aes(fill = smoking_status), col = "black") +
  geom_boxplot(col = "black", outlier.shape = NA, notch = T, width = 0.25) +
  geom_jitter(col = "black", alpha = 0.1, size = 1.2) +
  theme_minimal() +
  scale_fill_manual(values=c(safe_colorblind_palette[c(1,2)])) +
  scale_x_discrete(labels=c("Never\nsmokers", "Smokers")) + 
  ylab("Expression") + xlab("") + 
  facet_grid(. ~ "AdaptAge") + # Add the title strip
  stat_summary(aes(x = smoking_status, y = error),
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
  ylab("Adjusted Error")



# Figure S8c -----
# Get data 
data_lung_dam <- vadimAge$Lung_residuals_damge
data_lung_dam$smoking_status <- factor(ifelse(data_lung_dam$smoking_status == 0, "Never\nsmoker", "Smoker"))


c <- ggplot(data_lung_dam, aes(x=smoking_status, y=error)) +
  geom_violin(aes(fill = smoking_status), col = "black") +
  geom_boxplot(col = "black", outlier.shape = NA, notch = T, width = 0.25) +
  geom_jitter(col = "black", alpha = 0.1, size = 1.2) +
  theme_minimal() +
  scale_fill_manual(values=c(safe_colorblind_palette[c(1,2)])) +
  scale_x_discrete(labels=c("Never\nsmokers", "Smokers")) + 
  ylab("Expression") + xlab("") + 
  facet_grid(. ~ "DamAge") + # Add the title strip
  stat_summary(aes(x = smoking_status, y = error),
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
  ylab("Adjusted Error")


pdf(paste0(dir, "clocks.pdf"), w = 12, h = 4)
ggpubr::ggarrange(a,b,c, ncol = 3)
dev.off()
