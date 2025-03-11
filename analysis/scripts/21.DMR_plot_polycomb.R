# suppressMessages(suppressWarnings(library(bsseq)))

# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio

#Loading libraries
library(ggpubr)
library(valr)
library(dplyr)
library(ggplot2)

#ChromHMM IDs: https://personal.broadinstitute.org/cboix/epimap/metadata/Short_Metadata.html
# Caudate -> BSS00173 female
# Frontal_Cortex -> BSS00369 female
# Amygdala -> not exact match, so we take Hippocampus BSS01124
# BA24 -> not exact match, so we take cingulate gyrus BSS00219
# Hippocampus -> BSS01124 female
# Hypothalamus -> not exact match, so we take substantia nigra BSS01675
# Nucleus -> not exact match, so BSS00173 CAUDATE NUCLEUS
# Putamen -> BSS01469 male
# Thyroid -> BSS01832 female
# Lung -> BSS01190 female, BSS01835 male to compare if we wanted

#Functions:
#Option: >=(0.01), >(0), >0.1
my_fisher <- function(type, direction, threshold=0.1){
  #             DS      Not DS
  # type
  # other_types
  # print(type)
  chrom_tissue <- dmrs_regions
  
  type_df <- chrom_tissue[chrom_tissue$region_chromhmm_new == type,c("chr_dmr", "start_dmr", "end_dmr", "direction_dmr", "meanDiff_dmr", "n_dmr")] %>% distinct()
  other_type <- chrom_tissue[chrom_tissue$region_chromhmm_new != type,c("chr_dmr", "start_dmr", "end_dmr", "direction_dmr", "meanDiff_dmr", "n_dmr")] %>% distinct()
  if(direction=="hyper"){
    type_diff <- nrow(type_df[(type_df$meanDiff_dmr)>=(threshold),]) #At least 3 CpGs per region is already the background
    type_notdiff <- nrow(type_df[(type_df$meanDiff_dmr)<=(-threshold),]) #hypo
  }else if(direction=="hypo"){
    type_diff <- nrow(type_df[(type_df$meanDiff_dmr)<=(-threshold),]) #hypo
    type_notdiff <- nrow(type_df[(type_df$meanDiff_dmr)>=(threshold),]) #hyper
    
  }
  # type_notdiff <- nrow(type_df) - type_diff #One option was to use as background when testing hypermethylation, the hypomethylated and the DMRs that were slightly hypermethylated
  if(direction=="hyper"){
    other_type_diff <- nrow(other_type[(other_type$meanDiff_dmr)>=(threshold) & other_type$n_dmr>=3,])
    other_type_notdiff <- nrow(other_type[(other_type$meanDiff_dmr)<=(-threshold) & other_type$n_dmr>=3,])
  }else if(direction=="hypo"){
    other_type_diff <- nrow(other_type[(other_type$meanDiff_dmr)<=(-threshold) & other_type$n_dmr>=3,])
    other_type_notdiff <- nrow(other_type[(other_type$meanDiff_dmr)>=(threshold) & other_type$n_dmr>=3,])
  }
  # other_type_notdiff <- nrow(other_type) - other_type_diff

  ### test significance
  m <- matrix(c(type_diff, type_notdiff, other_type_diff, other_type_notdiff), 2, 2, byrow = T)
  # print(m)
  
  m[is.na(m)] <- 0
  #m <- m[c(type,paste0('No ',type)),]
  rownames(m) <- c(type, "Other")
  colnames(m) <- c("Hyper","Not Hyper")
  # print(m)
  f <- fisher.test(m)
  # print(f)
  return(list("f" = f, "m" = type_diff))
}


to_plot <- data.frame("oddsRatio"=1, "pvalue"=1, "CI_down"=1, "CI_up"=1, "type"=1, "sample_size"=1, "Tissue"=1)
tissues <- c("Hippocampus", "Frontal_Cortex", "Amygdala", "Hypothalamus", "BA24", "Caudate", "Nucleus", "Thyroid", "Lung") #Putamen has no DMRs
tissues <- c("Hippocampus", "Frontal_Cortex", "Amygdala", "Caudate", "Thyroid", "Lung") #Putamen has no DMRs
for(tissue in tissues){
  
  print(paste0("Starting code for tissue: ", tissue))
  Sys.time()
  
  #Reading data
  # BS_fit_subset_small <- readRDS(paste0('Data/',tissue,'_WGBS_CG/',tissue,'_small_smooth_all_samples_filt.CG.rds'))
  dmrs_small <- readRDS(paste0(tissue,"_dmrs_small.all.CG.rds")) #DMRs with more than 2 CpGs
  dmrs_small <- dmrs_small[dmrs_small$n>5,]
  # dmrs_small <- readRDS(paste0(tissue,"_dmrs_small.all.CG.nofilt.rds")) #DMRs with more than 2 CpGs
  dmrs0 <- dmrs_small
  #Putamen has 0 DMRs
  
  ### stats ####
  chr_number <- paste0('chr',c(1:22))
  dmrs_number <- dmrs_small
  
  table(dmrs_number$chr)
  
  sum(dmrs_number$meanDiff>0) 
  sum(dmrs_number$meanDiff<0) 
  table(dmrs_number$direction)
  
  # #### annotate with chromhmm #### 
  files <- list.files('ChromHMM/', pattern='.bed.gz',full.names=T)
  
  if(tissue=="Hippocampus"){
    id <- "BSS01124"
  } else if(tissue=="Frontal_Cortex"){
    id <- "BSS00369"
  } else if(tissue=="Hypothalamus"){ 
    id <- "BSS01675"
  } else if(tissue=="Amygdala"){ 
    id <- "BSS01124"
  } else if(tissue=="BA24"){
    id <- "BSS00219"
  } else if(tissue=="Caudate"){
    id <- "BSS00173"
  } else if(tissue=="Nucleus"){
    id <- "BSS00173"
  } else if(tissue=="Putamen"){
    id <- "BSS01469"
  } else if(tissue=="Thyroid"){
    id <- "BSS01832"
  } else if(tissue=="Lung"){
    id <- "BSS01190"
  }
  file <- files[grep(id, files)]
  chromhmm <- read.delim(file, sep='\t', header=F)
  
  
  chrom_df <- chromhmm[,c(1:4)]
  colnames(chrom_df) <- c('chrom','start','end','region')
  dmrs_number$chrom <- dmrs_number$chr
  dmrs0$chrom <- dmrs0$chr
  dmrs_regions <- bed_intersect(dmrs0, chrom_df, suffix = c("_dmr", "_chromhmm"))

  
  dmrs_regions$region_chromhmm_new <- dmrs_regions$region_chromhmm
  dmrs_regions$region_chromhmm_new[dmrs_regions$region_chromhmm %in% c("ReprPCWk","ReprPC", "EnhBiv", "TssBiv")] <- "ReprPC"
  dmrs_regions$region_chromhmm_new[!dmrs_regions$region_chromhmm %in% c("ReprPCWk","ReprPC", "EnhBiv", "TssBiv")] <- "Others"
  # dmrs_regions$region_chromhmm_new[dmrs_regions$region_chromhmm %in% c("ReprPC")] <- "ReprPC"
  # dmrs_regions$region_chromhmm_new[!dmrs_regions$region_chromhmm %in% c("ReprPC")] <- "Others"
  
  ### enrichment shared positions
  # Two-tailed Fisher test
  fisher_results <- my_fisher("ReprPC", "hypo", 0.1)
  to_plot <- rbind(to_plot, c("oddsRatio"=fisher_results$f$estimate, "pvalue"=fisher_results$f$p.value, "CI_down"=fisher_results$f$conf.int[1], "CI_up"=fisher_results$f$conf.int[2], "type"="hypo", "sample_size"=fisher_results$m, "Tissue"=tissue))
  
  fisher_results <- my_fisher("ReprPC", "hyper", 0.1)
  to_plot <- rbind(to_plot, c("oddsRatio"=fisher_results$f$estimate, "pvalue"=fisher_results$f$p.value, "CI_down"=fisher_results$f$conf.int[1], "CI_up"=fisher_results$f$conf.int[2], "type"="hyper", "sample_size"=fisher_results$m, "Tissue"=tissue))

}

colors_traits <- list('Smoking'=c("#CC6677", "#88CCEE"))
trait <- 'Smoking'

to_plot <- to_plot[-1,]
#Very important, the DMR software says something is hypo when group1 has lower methylation values than group2, so we need to switch what it considers as hypo and hyper to be consistent with the reference level we use
to_plot$type[to_plot$type =="hypo"] <- "Hypermethylation"
to_plot$type[to_plot$type =="hyper"] <- "Hypomethylation"

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
  scale_x_continuous(breaks=c(0, 2000, 4000)) 

p <- ggarrange(g, g2, labels = c("A", "B"),
               common.legend = TRUE, legend = "right", widths = c(0.8,0.3))
pdf(file = paste0("Plots/Enrichment.pdf"), w = 7.5, h = 3)
print(p)
dev.off()


#Add a barplot with the number of samples per tissue
annotation <- read.csv("../metadata.csv")
annotation <- annotation[annotation$Smoking %in% c(0,2),]

# Aggregate number of samples per Tissue and Smoking Status
smoking_counts <- annotation %>%
  group_by(Tissue, Smoking) %>%
  summarise(count = n()) %>%
  mutate(Smoking = factor(Smoking, levels = c(0, 2), labels = c("Never", "Smoker")))

smoking_counts$Tissue[smoking_counts$Tissue=="Brain - Caudate (basal ganglia)"] <- "Brain caudate"
smoking_counts$Tissue[smoking_counts$Tissue=="Brain - Frontal Cortex (BA9)"] <- "Brain frontal cortex"
smoking_counts$Tissue[smoking_counts$Tissue=="Brain - Amygdala"] <- "Brain amygdala"
smoking_counts$Tissue[smoking_counts$Tissue=="Brain - Hippocampus"] <- "Brain hippocampus"
smoking_counts$Tissue[smoking_counts$Tissue=="Brain - Anterior cingulate cortex (BA24)"] <- "Brain cingulate cortex"
smoking_counts$Tissue[smoking_counts$Tissue=="Brain - Hypothalamus"] <- "Brain hypothalamus"
smoking_counts$Tissue[smoking_counts$Tissue=="Brain - Nucleus accumbens (basal ganglia)"] <- "Brain nucleus accumbens"
smoking_counts$Tissue[smoking_counts$Tissue=="Brain - Putamen (basal ganglia)"] <- "Brain putamen"


smoking_counts <- smoking_counts[smoking_counts$Tissue!="Brain putamen",] #Not enough samples
smoking_counts <- smoking_counts[smoking_counts$Tissue!="Brain cingulate cortex",] #Not enough samples
smoking_counts <- smoking_counts[smoking_counts$Tissue!="Brain hypothalamus",] #Not enough samples
smoking_counts <- smoking_counts[smoking_counts$Tissue!="Brain nucleus accumbens",] #Not enough samples

# Create the smoking status bar plot
smoking_counts$Smoking <- factor(smoking_counts$Smoking, levels = c("Smoker", "Never"))
g3 <- ggplot(smoking_counts, aes(x = count, y = Tissue, fill = Smoking)) +
  geom_col(width = 0.6) +
  theme_classic() +
  xlab("Number of Samples") +
  ylab("") +
  scale_fill_manual(values = c("Smoker" = "#CC6677", "Never" = "#88CCEE"),
                    labels = c("Smokers", "Never smokers")) +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(colour="black", size=13),
        axis.text.y = element_text(colour="black", size=14),
        legend.text = element_text(colour="black", size=13),
        axis.title.x = element_text(size=16)) +
  guides(fill = guide_legend(reverse = TRUE))
g3

pdf(file = paste0("Plots/Samples.pdf"), w = 6, h = 3)
print(g3)
dev.off()
