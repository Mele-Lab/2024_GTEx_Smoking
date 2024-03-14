# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Generate publication Figure 5 (DNA methylation)


#Set path 
# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio

library(ggplot2)

#C
read_data <- function(variables, data){ #Function to prepare data to plot and compute adjusted p value
  
  odds_ratio <- lapply(variables, function(trait) data[[trait]][['f']]$estimate)
  adj.P.Val <- p.adjust(sapply(variables, function(trait) data[[trait]][['f']]$p.value), method = "BH")
  CI_down <- lapply(variables, function(trait) data[[trait]][['f']]$conf.int[1])
  CI_up <- lapply(variables, function(trait) data[[trait]][['f']]$conf.int[2])
  sample_size <- lapply(variables, function(trait) data[[trait]][['m']])
  
  
  names(odds_ratio) <- variables
  names(adj.P.Val) <- variables
  names(CI_down) <- variables
  names(CI_up) <- variables
  names(sample_size) <- variables
  
  odds_ratio_df <- as.data.frame(unlist(odds_ratio))
  odds_ratio_df$label <- variables
  odds_ratio_df$type <- deparse(substitute(data)) #Either hypo or hyper
  colnames(odds_ratio_df) <- c('oddsRatio', 'region','type')
  
  adj.P.Val_df <- as.data.frame(unlist(adj.P.Val))
  adj.P.Val_df$label <- variables
  adj.P.Val_df$type <- deparse(substitute(data))
  colnames(adj.P.Val_df) <- c('adjPvalue','region','type')
  
  CI_down_df <- as.data.frame(unlist(CI_down))
  CI_down_df$label <- variables
  CI_down_df$type <- deparse(substitute(data))
  colnames(CI_down_df) <- c('CI_down','region','type')
  
  CI_up_df <- as.data.frame(unlist(CI_up))
  CI_up_df$label <- variables
  CI_up_df$type <- deparse(substitute(data))
  colnames(CI_up_df) <- c('CI_up','region','type')
  
  sample_size_df <- as.data.frame(unlist(sample_size))
  sample_size_df$label <- variables
  sample_size_df$type <- deparse(substitute(data))
  colnames(sample_size_df) <- c('sample_size','region','type')
  
  all <- Reduce(function(x, y) merge(x, y, all=TRUE), list(odds_ratio_df, adj.P.Val_df, CI_down_df, CI_up_df, sample_size_df))
  head(all)
  all$sig <- 'not Sig'
  all$sig[all$adjPvalue<0.05] <- 'Sig'
  all <- all[,c("region","oddsRatio","adjPvalue","CI_down","CI_up","sig","type", "sample_size")]
  return(all)
}

hypo <- readRDS(paste0('data/enrichment_chromhmm_hypo_age.rds'))
hyper <- readRDS(paste0('data/enrichment_chromhmm_hyper_age.rds'))

hypo_d <- read_data(names(hypo), hypo)
hyper_d <- read_data(names(hyper), hyper)
hyper_hypo <- rbind(hypo_d, hyper_d)
#Change order to plot
hyper_hypo$region <- factor(hyper_hypo$region, levels=rev(c("Active TSS", "Flanking TSS", "Bivalent TSS",
                                                            "ZNF genes & repeats", "Heterochromatin", "Quiescent",
                                                            "Weak repressed polycomb", "Repressed polycomb",
                                                            "Weak transcription", "Strong transcription",
                                                            "Weak enhancer", "Active enhancer", "Genic enhancer", "Bivalent enhancer")))
hyper_hypo$type[hyper_hypo$type =="hypo"] <- "Hypomethylation"
hyper_hypo$type[hyper_hypo$type =="hyper"] <- "Hypermethylation"

hyper_hypo$sig[hyper_hypo$sig =="Sig"] <- "FDR < 0.05"
hyper_hypo$sig[hyper_hypo$sig =="not Sig"] <- "FDR >= 0.05"
hyper_hypo$sig <- factor(hyper_hypo$sig, levels=c("FDR >= 0.05", "FDR < 0.05"))

g <- ggplot(hyper_hypo, aes(x=log(oddsRatio), y=region, colour=type, alpha=sig)) + 
  geom_errorbar(aes(xmin=log(CI_down), xmax=log(CI_up)), width=.3) + 
  geom_vline(xintercept = 0) +
  geom_point() + ylab('') + theme_bw() +
  scale_colour_manual(values=c("#CC6677", "#88CCEE")) +
  xlab("Log(Odds ratio)") +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour="black", size=10),
        axis.text.y = element_text(colour="black", size=10),
        axis.title.x = element_text(size=10)) +
  scale_alpha_discrete(range = c(0.3, 1))

g2 <- ggplot(hyper_hypo) + geom_col(aes(sample_size, region, fill=type), width = 0.7) +
  theme_bw() + xlab("Number of DMPs") + ylab("") +
  scale_fill_manual(values=c("#CC6677", "#88CCEE")) +
  theme(legend.position = "none",
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour="black", size=9),
        axis.title.x = element_text(size=10))

pdf(file = paste0("figures/figure_5/chromHMM_age_smoking.pdf"), w = 5.5, h = 3)
ggarrange(g, g2, ncol=2, common.legend = TRUE, legend="right", widths = c(1, 0.4))
dev.off()

#D is in figure_S4.R