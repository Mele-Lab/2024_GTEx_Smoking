#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Methylation results analysis on reversibility
# @software version: R=4.2.2

#Set path 
# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")


tissue_info <- read.csv("data/public/tissue_abreviation.txt")
tissues <- c("Lung", "ColonTransverse")

data <- matrix(nrow = length(tissues), ncol=3, dimnames = list(tissues, c("Reversible", "Partially reversible", "Non reversible")))
percentage_closer_to_never <- c()
logFC_never_ex <- c()
logFC_ex_smokers <- c()
for(tissue in tissues){ 
  print(tissue)
  model <- readRDS(paste0("tissues/", tissue, "/DML_results.rds"))
  signif <- model$Smoking2[model$Smoking2$adj.P.Val<0.05,]
  if(nrow(signif)==0){
    data[tissue,] <- rep(NA, 3)
    next
  }
  tissue_id <- tissue_info$Name[tissue_info$tissue==tissue]
  
  
  signif1 <- model$Smoking1[model$Smoking1$adj.P.Val<0.05,] #Never vs ex
  signif2 <- model$"Smoking1-Smoking2"[model$"Smoking1-Smoking2"$adj.P.Val<0.05,] #Ex vs smokers
  reversed <- rownames(signif)[rownames(signif) %in% rownames(signif2)]
  reversed <- reversed[!reversed %in% rownames(signif1)]
  
  data[tissue, 1] <- length(reversed) #reversed
  
  non_reversed <- rownames(signif)[rownames(signif) %in% rownames(signif1)]
  non_reversed <- non_reversed[!non_reversed %in% rownames(signif2)]
  data[tissue, 3] <- length(non_reversed) #non reversed
  
  data[tissue, 2] <- length(rownames(signif)[!rownames(signif) %in% reversed & !rownames(signif) %in% non_reversed])
  

  #Getting logFC between never smokers and ex-smokers and ex-smokers and smokers for the partially reversible positions
  to_exclude <- rownames(model$Smoking2)[model$Smoking2$adj.P.Val>=0.05] #these will not be partially reversible
  never_vs_ex <- model$Smoking1
  to_exclude <- c(to_exclude, rownames(never_vs_ex)[never_vs_ex$adj.P.Val<0.05])
  ex_vs_smokers <- model$`Smoking1-Smoking2`
  to_exclude <- c(to_exclude, rownames(ex_vs_smokers)[ex_vs_smokers$adj.P.Val<0.05])
  #We only want the logFC for the partialy reversible
  never_vs_ex <- never_vs_ex[!rownames(never_vs_ex) %in% to_exclude,]
  ex_vs_smokers <- ex_vs_smokers[!rownames(ex_vs_smokers) %in% to_exclude,]
  
  logFC_never_ex <- c(logFC_never_ex, abs(never_vs_ex$logFC))
  logFC_ex_smokers <- c(logFC_ex_smokers, abs(ex_vs_smokers$logFC))

  
  booleans <- abs(never_vs_ex$logFC) < abs(ex_vs_smokers$logFC)
  # booleans <- sum(booleans)/length(ex_vs_smokers$logFC)
  # percentage_closer_to_never <- c(percentage_closer_to_never, sum(booleans))
  percentage_closer_to_never <- c(percentage_closer_to_never, booleans)
}

#Percentage of partially reversed positions
round((sum(data[,2])/sum(data))*100, 2)

#Percentahe of ex-smokers in partially reversed positions closer to smokers:
1-sum(percentage_closer_to_never)/length(percentage_closer_to_never)


# test <- wilcox.test(logFC_never_ex, logFC_ex_smokers, paired=T, alternative = "greater")
test <- wilcox.test(logFC_never_ex, logFC_ex_smokers, paired=T)
p_val <- test$p.value

logFC <- cbind(logFC_never_ex, logFC_ex_smokers)
figure_S9 <- list(logFC, data, test)
names(figure_S9) <- c("logFC", "data", "test")
saveRDS(figure_S9, paste0("../figures/data/figure_S8_logFC_reversible.rds"))


#Are the reversible genes enriched in something?
lung <- readRDS(paste0("tissues/Lung/DML_results.rds"))
signif <- rownames(lung$Smoking2)[lung$Smoking2$adj.P.Val<0.05]
non_reversible <- signif[signif %in% rownames(lung$Smoking1)[lung$Smoking1$adj.P.Val<0.05]]
reversible <- signif[signif %in% rownames(lung$`Smoking1-Smoking2`)[lung$`Smoking1-Smoking2`$adj.P.Val<0.05]]
partially_reversible <- signif[!signif %in% c(reversible, non_reversible)]

#Are the rev or non_rev enriched in hypo or hyper?

#Are the rev or non_rev enriched in chromHMM?
annotation <- read.csv(paste0("~/Documents/mn4/scratch/bsc83/bsc83535/GTEx/v9/Oliva/", "GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv"))
chromhmm <- read.delim('~/Documents/mn4/projects/bsc83/Projects/GTEx_v8/Methylation/Data/EpiMap/Lung.BSS01190_18_CALLS_segments.bed.gz', sep='\t', header=F) #One donor of lung from EpiMap to run ChromHMM #for colon we should use another file
ann_bed <- annotation[!is.na(annotation$MAPINFO) & !is.na(annotation$CHR),] %>%
  dplyr::select(chrom=CHR, start=MAPINFO, end=MAPINFO, name=IlmnID) %>% distinct()
ann_bed$chrom <- paste0('chr',ann_bed$chrom)
head(ann_bed)
library(valr)
chrom_df <- chromhmm[,c(1:4)]
colnames(chrom_df) <- c('chrom','start','end','region')
chromhmm_cpgs <- bed_intersect(ann_bed, chrom_df, suffix = c("_ann", "_chromhmm"))
chromhmm_cpgs$region_chromhmm[chromhmm_cpgs$region_chromhmm %in% c("TssFlnkD", "TssFlnk", "TssFlnkU")] <- "Flanking TSS"
chromhmm_cpgs$region_chromhmm[chromhmm_cpgs$region_chromhmm %in% c("TssBiv")] <- "Bivalent TSS"
chromhmm_cpgs$region_chromhmm[chromhmm_cpgs$region_chromhmm %in% c("TssA")] <- "Active TSS"
chromhmm_cpgs$region_chromhmm[chromhmm_cpgs$region_chromhmm %in% c("EnhA2", "EnhA1")] <- "Active enhancer"
chromhmm_cpgs$region_chromhmm[chromhmm_cpgs$region_chromhmm %in% c("EnhWk")] <- "Weak enhancer"
chromhmm_cpgs$region_chromhmm[chromhmm_cpgs$region_chromhmm %in% c("EnhBiv")] <- "Bivalent enhancer"
chromhmm_cpgs$region_chromhmm[chromhmm_cpgs$region_chromhmm %in% c("EnhG1", "EnhG2")] <- "Genic enhancer"
chromhmm_cpgs$region_chromhmm[chromhmm_cpgs$region_chromhmm %in% c("ReprPCWk")] <- "Weak repressed polycomb"
chromhmm_cpgs$region_chromhmm[chromhmm_cpgs$region_chromhmm %in% c("ReprPC")] <- "Repressed polycomb"
chromhmm_cpgs$region_chromhmm[chromhmm_cpgs$region_chromhmm %in% c("Quies")] <- "Quiescent"
chromhmm_cpgs$region_chromhmm[chromhmm_cpgs$region_chromhmm %in% c("Het")] <- "Heterochromatin"
chromhmm_cpgs$region_chromhmm[chromhmm_cpgs$region_chromhmm %in% c("ZNF/Rpts")] <- "ZNF genes & repeats"
chromhmm_cpgs$region_chromhmm[chromhmm_cpgs$region_chromhmm %in% c("TxWk")] <- "Weak transcription"
chromhmm_cpgs$region_chromhmm[chromhmm_cpgs$region_chromhmm %in% c("Tx")] <- "Strong transcription"


results <- data.frame("subset"="1", "region"="1", "OR"="1", "p.value"="1")
for(subset in c("reversible", "partially_reversible", "non_reversible")){
  for(region in unique(chromhmm_cpgs$region_chromhmm)){
    # bg <- chromhmm_cpgs$name_ann #All cpgs studied
    bg <- signif #Smoking-DMPs
    print(paste(subset, region))
    region_cpgs <- chromhmm_cpgs$name_ann[chromhmm_cpgs$region_chromhmm==region]
    common <- region_cpgs[region_cpgs %in% get(subset)]
    only_region <- region_cpgs[!region_cpgs %in% get(subset)]
    only_subset <- get(subset)[!get(subset) %in% region_cpgs]
    only_region <- only_region[only_region %in%bg]
    only_subset <- only_subset[only_subset %in%bg]
    common <- common[common %in%bg]
    m <- matrix(c(length(common), length(only_correlated),
                  length(only_subset), length(bg)-length(only_subset)-length(only_region)-length(common) ), nrow=2)
    print(m)
    test<- fisher.test(m)
    print(test)
    results <- rbind(results, cbind("subset"=subset, "region"=region, "OR"=test$estimate, "p.value"=test$p.value))
  }
}
results <- results[-1,]
results$p.adjusted <- p.adjust(results$p.value, method = "BH")
results_subset <- results
results_subset <- results_subset[results_subset$p.adjusted<0.05,]
results_subset <- results_subset[results_subset$subset!="partially_reversible",]
#Are the rev or non_rev enriched in functions?
library(missMethyl)
for(subset in c("reversible", "partially_reversible", "non_reversible")){
  print(subset)
  go_results <- gometh(sig.cpg=get(subset), all.cpg = signif, collection="GO", array.type="EPIC") # Background: smoking-DMPs
  print(go_results[go_results$FDR<0.05,])
}

#Are the rev or non_rev enriched in correlated positions?
cor <- readRDS("tissues/Lung/Correlations_DMP_DEG_new.rds")
bg <- cor$probe #DMPs paired to a DEG
correlated <- cor$probe[cor$p.adj<0.05]
for(subset in c("reversible", "partially_reversible", "non_reversible")){
  print(subset)
  common <- correlated[correlated %in% get(subset)]
  only_correlated <- correlated[!correlated %in% get(subset)]
  only_subset <- get(subset)[!get(subset) %in% correlated]
  only_subset <- only_subset[only_subset %in% bg]
  m <- matrix(c(length(common), length(only_correlated),
                length(only_subset), length(bg)-length(only_subset)-length(only_correlated)-length(common) ), nrow=2)
  print(m)
  print(fisher.test(m))
}
p.adjust(c(1.054863e-16, 7.438e-07, 0.03965))

#Are the rev or non_rev enriched in age-DMPs?
age <- rownames(lung$AGE)[lung$AGE$adj.P.Val<0.05]
for(subset in c("reversible", "partially_reversible", "non_reversible")){
  bg <- signif
  print(subset)
  common <- age[age %in% get(subset)]
  only_age <- correlated[!correlated %in% get(subset)]
  only_subset <- get(subset)[!get(subset) %in% age]
  only_subset <- only_subset[only_subset %in% bg]
  only_age <- only_age[only_age %in% bg]
  m <- matrix(c(length(common), length(only_age),
                length(only_subset), length(bg)-length(only_subset)-length(only_age)-length(common) ), nrow=2)
  print(m)
  print(fisher.test(m))
}
  
#  no enrichment in functions, no enrichment in correlations, but enrichment in age-DMPs



#The chromHMM shows similar enrichment for both reversible and non-reversible, except quiescent in non-reversible, which makes sense
#Let's see if PR closer to smokers show more enrichments in chromHMM. Well, all CpGs where ex-smokers are closer to smokers vs all where ex-smokers are closer to never-smokers
model <- readRDS(paste0("tissues/Lung/DML_results.rds"))
signif <- model$Smoking2[model$Smoking2$adj.P.Val<0.05,]

#Getting logFC between never smokers and ex-smokers and ex-smokers and smokers for the partially reversible positions
never_vs_ex <- model$Smoking1
ex_vs_smokers <- model$`Smoking1-Smoking2`

#We only want the logFC for the partialy reversible
never_vs_ex <- never_vs_ex[rownames(never_vs_ex) %in% rownames(signif),]
ex_vs_smokers <- ex_vs_smokers[rownames(ex_vs_smokers) %in% rownames(signif),]

never_vs_ex <- cbind(rownames(never_vs_ex), abs(never_vs_ex$logFC))
colnames(never_vs_ex) <- c("CpG", "never_vs_ex")
ex_vs_smokers <- cbind(rownames(ex_vs_smokers), abs(ex_vs_smokers$logFC))
colnames(ex_vs_smokers) <- c("CpG", "ex_vs_smokers")

merged <- merge(never_vs_ex, ex_vs_smokers, by="CpG")
merged$ex_vs_smokers <- as.numeric(merged$ex_vs_smokers)
merged$never_vs_ex <- as.numeric(merged$never_vs_ex)

merged$closer_to_smoker <- merged$never_vs_ex > merged$ex_vs_smokers
to_never <- merged$CpG[!merged$closer_to_smoker]
to_smoker <- merged$CpG[merged$closer_to_smoker]

results <- data.frame("subset"="1", "region"="1", "OR"="1", "CI_lower"="1", "CI_upper"="1", "p.value"="1", "sample_size"=1)
for(subset in c("to_never", "to_smoker")){
  for(region in unique(chromhmm_cpgs$region_chromhmm)){
    # bg <- chromhmm_cpgs$name_ann #All cpgs studied
    bg <- rownames(signif) #Smoking-DMPs
    print(paste(subset, region))
    region_cpgs <- chromhmm_cpgs$name_ann[chromhmm_cpgs$region_chromhmm==region]
    common <- region_cpgs[region_cpgs %in% get(subset)]
    only_region <- region_cpgs[!region_cpgs %in% get(subset)]
    only_subset <- get(subset)[!get(subset) %in% region_cpgs]
    only_region <- only_region[only_region %in%bg]
    only_subset <- only_subset[only_subset %in%bg]
    common <- common[common %in%bg]
    m <- matrix(c(length(common), length(only_region),
                  length(only_subset), length(bg)-length(only_subset)-length(only_region)-length(common) ), nrow=2)
    print(m)
    test<- fisher.test(m)
    print(test)
    results <- rbind(results, cbind("subset"=subset, "region"=region, "OR"=test$estimate, "CI_lower"=test$conf.int[1], "CI_upper"=test$conf.int[2], "p.value"=test$p.value, "sample_size"=length(common)))
  }
}
results <- results[-1,]
results$p.adjusted <- p.adjust(results$p.value, method = "BH")
# results_subset <- results
# results_subset <- results_subset[results_subset$p.adjusted<0.05,]
# results_subset <- results_subset[results_subset$subset!="partially_reversible",]


results$sig <- "FDR >= 0.05"
results$sig[results$p.adjusted<0.05] <- "FDR < 0.05"
results$sig <- factor(results$sig, levels=c("FDR >= 0.05", "FDR < 0.05"))
results$OR <- as.numeric(results$OR)
results$CI_lower <- as.numeric(results$CI_lower)
results$CI_upper <- as.numeric(results$CI_upper)
results$sample_size <- as.numeric(results$sample_size)
results$region <- factor(results$region, levels=rev(c("Active TSS", "Flanking TSS", "Bivalent TSS",
                                                            "ZNF genes & repeats", "Heterochromatin", "Quiescent",
                                                            "Weak repressed polycomb", "Repressed polycomb",
                                                            "Weak transcription", "Strong transcription",
                                                            "Weak enhancer", "Active enhancer", "Genic enhancer", "Bivalent enhancer")))


library(ggplot2)
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

library(ggpubr)
ggarrange(g1, g2, ncol=2, common.legend = TRUE, legend="right", widths = c(1, 0.4))
