# suppressMessages(suppressWarnings(library(bsseq)))

# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio

#Loading libraries
library(ggpubr)
library(valr)
library(dplyr)
library(ggplot2)
library(GenomicRanges)


#Option regions that are CpG Islands:
library(annotatr)
annotations <- build_annotations(genome = "hg38", annotations = "hg38_cpg_islands")



to_plot <- data.frame("oddsRatio"=1, "pvalue"=1, "CI_down"=1, "CI_up"=1, "type"=1, "sample_size"=1, "Tissue"=1)
tissues <- c("Hippocampus", "Frontal_Cortex", "Amygdala", "Hypothalamus", "BA24", "Caudate", "Nucleus", "Thyroid", "Lung") #Putamen has no DMRs
tissues <- c("Hippocampus", "Frontal_Cortex", "Amygdala", "Caudate", "Thyroid", "Lung") #Putamen has no DMRs
for(tissue in tissues){
  
  print(paste0("Starting code for tissue: ", tissue))
  Sys.time()
  
  dmrs <- readRDS(paste0(tissue,"_dmrs_small.all.CG.rds")) #DMRs with more than 2 CpGs
  dmrs <- dmrs[dmrs$n>5,]
  hypo <- dmrs[dmrs$meanDiff>0.1,] #The reference is not the one expected
  hyper <- dmrs[dmrs$meanDiff<(-0.1),] #The reference is not the one expected
  
  #Transform intro granges:
  hypo_gr <- GRanges(seqnames = hypo$chr, ranges = IRanges(hypo$start, hypo$end))
  hyper_gr <- GRanges(seqnames = hyper$chr, ranges = IRanges(hyper$start, hyper$end))
  
  annotated_hypo <- annotate_regions(regions = hypo_gr, annotations = annotations, ignore.strand = TRUE)
  annotated_hyper <- annotate_regions(regions = hyper_gr, annotations = annotations, ignore.strand = TRUE)
  
  hyper_CGI <- length(annotated_hyper$annot)
  hypo_CGI <- length(annotated_hypo$annot)
  
  hyper_no_CGI <- length(hyper_gr) - hyper_CGI
  hypo_no_CGI <- length(hypo_gr) - hypo_CGI
  

  # Fisher's exact test
  fisher_results <- fisher.test(matrix(c(hyper_CGI, hypo_CGI, hyper_no_CGI, hypo_no_CGI), ncol=2))
  to_plot <- rbind(to_plot, c("oddsRatio"=fisher_results$estimate, "pvalue"=fisher_results$p.value, "CI_down"=fisher_results$conf.int[1], "CI_up"=fisher_results$conf.int[2], "type"="hyper", "sample_size"=hyper_CGI, "Tissue"=tissue))
  
  fisher_results <- fisher.test(matrix(c(hypo_CGI, hyper_CGI, hypo_no_CGI, hyper_no_CGI), ncol=2))
  to_plot <- rbind(to_plot, c("oddsRatio"=fisher_results$estimate, "pvalue"=fisher_results$p.value, "CI_down"=fisher_results$conf.int[1], "CI_up"=fisher_results$conf.int[2], "type"="hypo", "sample_size"=hypo_CGI, "Tissue"=tissue))
}

colors_traits <- list('Smoking'=c("#CC6677", "#88CCEE"))
trait <- 'Smoking'

to_plot <- to_plot[-1,]
to_plot$type[to_plot$type =="hyper"] <- "Hypermethylation"
to_plot$type[to_plot$type =="hypo"] <- "Hypomethylation"

to_plot$oddsRatio <- as.numeric(to_plot$oddsRatio)
to_plot$CI_up <- as.numeric(to_plot$CI_up)
to_plot$CI_down <- as.numeric(to_plot$CI_down)
to_plot$pvalue <- as.numeric(to_plot$pvalue)
to_plot$sample_size <- as.numeric(to_plot$sample_size)

to_plot$padj <- p.adjust(to_plot$pvalue, method = "BH")

to_plot$sig <- to_plot$padj<0.05
to_plot$sig[to_plot$sig ==T] <- "FDR < 0.05"
to_plot$sig[to_plot$sig ==F] <- "FDR >= 0.05"
to_plot$sig <- factor(to_plot$sig, levels=c("FDR >= 0.05", "FDR < 0.05"))

to_plot$Tissue[to_plot$Tissue=="Caudate"] <- "Brain caudate"
to_plot$Tissue[to_plot$Tissue=="Frontal_Cortex"] <- "Brain frontal cortex"
to_plot$Tissue[to_plot$Tissue=="Amygdala"] <- "Brain amygdala"
to_plot$Tissue[to_plot$Tissue=="Hippocampus"] <- "Brain hippocampus"
to_plot$Tissue[to_plot$Tissue=="BA24"] <- "Brain cingulate cortex"
to_plot$Tissue[to_plot$Tissue=="Hypothalamus"] <- "Brain hypothalamus"
to_plot$Tissue[to_plot$Tissue=="Nucleus"] <- "Brain nucleus accumbens"

g <- ggplot(to_plot, aes(x=log2(oddsRatio), y=Tissue, colour=type, alpha=sig)) +
  geom_errorbar(aes(xmin=log2(CI_down), xmax=log2(CI_up)), width=.3) +
  geom_vline(xintercept = 0) +
  #xlim(0,20) + #Only for Lung to show the 0
  geom_point(size=3) + ylab('') + theme_bw() +
  scale_colour_manual(values=colors_traits[[trait]]) +
  xlab("log2(Odds ratio)") +
  scale_alpha_discrete(range = c(0.4, 1), drop = FALSE) +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(colour="black", size=13),
        axis.text.y = element_text(colour="black", size=14),
        legend.text = element_text(colour="black", size=13),
        axis.title.x = element_text(size=16),
        legend.spacing.y = unit(-0.05, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth=1))


# #Plot sample sizes:
g2 <- ggplot(to_plot) + geom_col(aes(sample_size, Tissue, fill=type), width = 0.6) +
  theme_classic() + xlab("Number of DMRs") + ylab("") +
  scale_fill_manual(values=colors_traits[[trait]]) +
  theme(legend.position = "none",
        axis.text.x = element_text(colour="black", size=13),
        axis.text.y=element_blank(),
        axis.title.x = element_text(size=16)) +
  scale_x_continuous(breaks=c(0, 750, 1500)) 

p <- ggarrange(g, g2, labels = c("A", "B"),
               common.legend = TRUE, legend = "right", widths = c(0.8,0.3))
# pdf(file = paste0("Plots/enrichment", tissue, "_0.1_Winona_all.pdf"), w = 8, h = 4)
print(p)
# dev.off()

pdf(file = paste0("Plots/Enrichment_CGI.pdf"), w = 7.5, h = 3)
print(p)
dev.off()



#Add a barplot with the number of samples per tissue
# annotation <- read.csv("../metadata.csv")
# annotation <- annotation[annotation$Smoking %in% c(0,2),]
# 
# # Aggregate number of samples per Tissue and Smoking Status
# smoking_counts <- annotation %>%
#   group_by(Tissue, Smoking) %>%
#   summarise(count = n()) %>%
#   mutate(Smoking = factor(Smoking, levels = c(0, 2), labels = c("Never", "Smoker")))
# 
# smoking_counts$Tissue[smoking_counts$Tissue=="Brain - Caudate (basal ganglia)"] <- "Brain caudate"
# smoking_counts$Tissue[smoking_counts$Tissue=="Brain - Frontal Cortex (BA9)"] <- "Brain frontal cortex"
# smoking_counts$Tissue[smoking_counts$Tissue=="Brain - Amygdala"] <- "Brain amygdala"
# smoking_counts$Tissue[smoking_counts$Tissue=="Brain - Hippocampus"] <- "Brain hippocampus"
# smoking_counts$Tissue[smoking_counts$Tissue=="Brain - Anterior cingulate cortex (BA24)"] <- "Brain cingulate cortex"
# smoking_counts$Tissue[smoking_counts$Tissue=="Brain - Hypothalamus"] <- "Brain hypothalamus"
# smoking_counts$Tissue[smoking_counts$Tissue=="Brain - Nucleus accumbens (basal ganglia)"] <- "Brain nucleus accumbens"
# smoking_counts$Tissue[smoking_counts$Tissue=="Brain - Putamen (basal ganglia)"] <- "Brain putamen"
# 
# 
# smoking_counts <- smoking_counts[smoking_counts$Tissue!="Brain putamen",] #No DMRs
# # Create the smoking status bar plot
# g3 <- ggplot(smoking_counts, aes(x = count, y = Tissue, fill = Smoking)) +
#   geom_col(width = 0.6) +
#   theme_classic() +
#   xlab("Number of Samples") +
#   ylab("") +
#   scale_fill_manual(values = c("Never" = "#88CCEE", "Smoker" = "#CC6677"), 
#                     labels = c("Never smokers", "Smokers")) +
#   theme(legend.title = element_blank(),
#         axis.text.x = element_text(colour="black", size=13),
#         axis.text.y = element_text(colour="black", size=14),
#         legend.text = element_text(colour="black", size=13),
#         axis.title.x = element_text(size=16))
# g3
# 
