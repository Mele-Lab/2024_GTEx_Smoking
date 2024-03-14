#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Generate publication Figure S8 (Reversibility on Methylation)


#Set path 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Loading libraries
library(ggplot2)
library(ggpubr)

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

ggarrange(g1, g2, ncol=2, common.legend = TRUE, legend="right", widths = c(1, 0.4))
