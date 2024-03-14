#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Methylation results analysis on aging and smoking
# @software version: R=4.2.2

#Set path 
# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")

tissue_info <- read.csv("data/public/tissues_sorted.csv")
source("../figures/safe_colourblind_pallete.R")

# Number of DML per tissue
results_DML <- readRDS("tissues/Lung/DML_results.rds")
# results_DML <- readRDS("tissues/ColonTransverse/DML_results.rds")
smoking <- results_DML$Smoking2[results_DML$Smoking2$adj.P.Val<0.05,]

fisher_function <- function(variable){
  common <- rownames(variable)[rownames(variable) %in% rownames(smoking)]
  only_smoking <- rownames(smoking)[!rownames(smoking) %in% rownames(variable)]
  only_variable <- rownames(variable)[!rownames(variable) %in% rownames(smoking)]
  background <- rownames(results_DML$Smoking2)
  bg <- background[!background %in% common]
  bg <- background[!background %in% only_variable]
  bg <- background[!background %in% only_smoking]
  m <- matrix(c(length(common), length(only_variable), length(only_smoking), length(bg)), nrow=2)
  f <- fisher.test(m, alternative = "greater") 
  return(f)
}

age <- results_DML$AGE[results_DML$AGE$adj.P.Val<0.05,]
fisher_function(age)

sex <- results_DML$SEX[results_DML$SEX$adj.P.Val<0.05,]
fisher_function(sex)

bmi <- results_DML$BMI[results_DML$BMI$adj.P.Val<0.05,]
fisher_function(bmi)

#Chi-square test -> directionality and plot
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)

bias.fun <- function(tissue){
  # bias.fun(my_names[i], variable, my_acronyms[i], tissues[i])
  print(tissue)
  results_DML <- readRDS(paste0("tissues/", tissue, "/DML_results.rds"))
  # results_DML <- readRDS(paste0(project_path, "ColonTransverse/DML_results.rds"))
  smoking <- results_DML$Smoking2[results_DML$Smoking2$adj.P.Val<0.05,]
  age <- results_DML$AGE[results_DML$AGE$adj.P.Val<0.05,]
  # Do we observe a higher than expected overlap of DMPs in a particular direction of change?
  genes.overlap <- intersect(rownames(smoking), rownames(age))
  trait1.up <- rownames(age[age$logFC>0,])
  trait1.down <- rownames(age[age$logFC<0,])
  
  trait2.up <- rownames(smoking[smoking$logFC>0,])
  trait2.down <- rownames(smoking[smoking$logFC<0,])
  
  # Observed counts
  counts <- c(sum(trait1.up %in% trait2.up), # upup
              sum(trait1.down %in% trait2.up), # downup
              sum(trait1.up %in% trait2.down), # updown
              sum(trait1.down %in% trait2.down) # downdown
  )
  # Expected proportions
  trait1.up.p <- length(trait1.up)/length(c(trait1.up, trait1.down)) # or / nrow(age)
  trait1.down.p <- length(trait1.down)/length(c(trait1.up, trait1.down))
  trait2.up.p <- length(trait2.up)/length(c(trait2.up, trait2.down))
  trait2.down.p <- length(trait2.down)/length(c(trait2.up, trait2.down))
  expected_prob <- c(trait1.up.p * trait2.up.p, trait1.down.p * trait2.up.p, trait1.up.p * trait2.down.p, trait1.down.p * trait2.down.p)
  expected_counts <- expected_prob*length(genes.overlap)
  # print(expected_counts)
  if(length(genes.overlap)==0){
    return(list("p-value" = NA, "which" = NA,
                "OvsE" = "NA:NA:NA:NA",
                "observed" = paste(as.character(counts), collapse = "/"),
                "enough" = F,
                "twenty" = F
    ))
  } 
  
  # Chi-square goodness of fit test
  if(min(expected_counts)<5){
    chi.test <- chisq.test(counts, p = expected_prob, simulate.p.value = T)
    enough <- F
  }else{
    chi.test <- chisq.test(counts, p = expected_prob)
    enough <- T
  }
  return(list("p-value" = chi.test$p.value,
              "which" = paste(c("up_up", "down_up", "up_down", "down_down")[chi.test$observed/round(chi.test$expected) > 1],collapse = ":"),
              "OvsE" = paste0(as.character(round(chi.test$observed/round(chi.test$expected),2)),collapse = ":"),
              "observed" = paste(as.character(counts), collapse = "/"),
              "enough" = enough
  ))
}

# Enrichment analysis  --
chi_square.results <- lapply(c("Lung", "ColonTransverse"), function(variable)
  bias.fun(variable))
names(chi_square.results) <- c("Lung", "ColonTransverse")


# col_names <- c("old-smoker", "young-smoker", "old-never smoker", "young-never smoker")
col_names <- c("up - up", "down - up", "up - down", "down - down")
data1 <- as.matrix(sapply(c("Lung", "ColonTransverse"), function(tissue) chi_square.results[[tissue]][["observed"]]))
data1 <- cbind(sapply(strsplit(data1,'/'), "[", 1), 
               sapply(strsplit(data1,'/'), "[", 2),
               sapply(strsplit(data1,'/'), "[", 3),
               sapply(strsplit(data1,'/'), "[", 4))
data1 <- as.data.frame(data1)
data1[c(1:4)] <- sapply(data1[c(1:4)],as.numeric) #I add c(1:4) to keep rownames (not necessary now)
data1 <- t(data1)
rownames(data1) <- col_names
# colnames(data1) <- "Lung"
colnames(data1) <- c("Lung", "ColonTransverse")
data1 <- as.data.frame(data1)
data1$type <- rownames(data1)
data3 <- as.data.frame(t(data1)[1:ncol(data1)-1,])
# data3[,1] <- as.numeric(data3[,1])
# data3 <- data3/colSums(data3)
data3[c(1:ncol(data3))] <- sapply(data3[c(1:ncol(data3))],as.numeric) #I add [c(1:4)] to keep rownames
data3 <- data3/rowSums(data3)
# colnames(data3) <- "Lung"
# colnames(data3) <- c("Lung", "ColonTransverse")
data3$type <- rownames(data3)

library(reshape2)
data <- melt(data1)
# data <- data$type
melted <- melt(data3)[c(1,3,5,7,2,4,6,8),]
data$y <- 100*melted[,3] #Sorted by variable instead of type
data$type <- factor(data$type, levels = rev(c("up - up", "down - down", "down - up", "up - down")), order = T)
library(colorspace) #To make lighter colors
saveRDS(data, "../figures/data/figure_6_bias_age.rds")


#The overlap between smoking and age is enriched in which locations?
#Check 13.4.location_
read_data <- function(variables, data){
  odds_ratio <- list()
  adj.P.Val <- list()
  CI_down <- list()
  CI_up <- list()
  
  odds_ratio <- lapply(variables, function(trait) data[[trait]][['f']]$estimate)
  adj.P.Val <- p.adjust(sapply(variables, function(trait) data[[trait]][['f']]$p.value), method = "BH")
  CI_down <- lapply(variables, function(trait) data[[trait]][['f']]$conf.int[1])
  CI_up <- lapply(variables, function(trait) data[[trait]][['f']]$conf.int[2])
  
  names(odds_ratio) <- variables
  names(adj.P.Val) <- variables
  names(CI_down) <- variables
  names(CI_up) <- variables
  
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
  
  all <- merge(odds_ratio_df, merge(adj.P.Val_df, merge(CI_down_df, CI_up_df, by=c('region', 'type')),  by=c('region', 'type')),  by=c('region', 'type'))
  head(all)
  all$sig <- 'not Sig'
  all$sig[all$adjPvalue<0.05] <- 'Sig'
  all <- all[,c("region","oddsRatio","adjPvalue","CI_down","CI_up","sig","type")]
  return(all)
}

# hypo <- readRDS('output/enrichment_chromhmm_hypo_age_CI.rds')
# hyper <- readRDS('output/enrichment_chromhmm_hyper_age_CI.rds')
# 
# hypo_d <- read_data(names(hypo), hypo)
# hyper_d <- read_data(names(hyper), hyper)
# hyper_hypo <- rbind(hypo_d, hyper_d)

# g <- ggplot(hyper_hypo[hyper_hypo$adjPvalue<0.05,], aes(x=(oddsRatio), y=region, colour=type)) + 
#   geom_errorbar(aes(xmin=(CI_down), xmax=(CI_up)), width=.15) + #xlim(c(0,3))+
#   geom_vline(xintercept = 1) + xlab("Odds ratio") +
#   geom_point() + ylab('') + theme_bw() +
#   theme(legend.title=element_blank())





# 
# 
# 
# #TFBS and Chromatin states
# hypo <- readRDS('output/enrichment_anno_age_hypo.rds')
# hyper <- readRDS('output/enrichment_anno_age_hyper.rds')
# 
# hypo_d <- read_data(c("tfbs", "open_chrm"), hypo)
# hyper_d <- read_data(c("tfbs", "open_chrm"), hyper)
# hyper_hypo <- rbind(hypo_d, hyper_d)
# g <- ggplot(hyper_hypo[hyper_hypo$adjPvalue<0.05,], aes(x=(oddsRatio), y=region, colour=type)) + 
#   geom_errorbar(aes(xmin=(CI_down), xmax=(CI_up)), width=.15) + #xlim(c(0,3))+
#   geom_vline(xintercept = 1) + xlab("Odds ratio") +
#   geom_point() + ylab('') + theme_bw() +
#   scale_y_discrete(breaks=c("tfbs", "open_chrm"),
#                    labels=c("TFBS", "Open Chromatin")) +
#   theme(legend.title=element_blank())
# 
# ggsave("../figures/figures/figure_6/OR_age_tfbs.png", g, device="png", width = 4, height = 1.5)
# 
# 
# 
# 
# #Promoters or enhancers
# names <- c("promoter", "enhancer", "gene_body")
# hypo_d <- read_data(names, hypo)
# hyper_d <- read_data(names, hyper)
# hyper_hypo <- rbind(hypo_d, hyper_d)
# hyper_hypo$region <- factor(hyper_hypo$region, levels=rev(c("promoter", "enhancer", "gene_body")))
# 
# g <- ggplot(hyper_hypo[hyper_hypo$adjPvalue<0.05,], aes(x=(oddsRatio), y=region, colour=type)) + 
#   geom_errorbar(aes(xmin=(CI_down), xmax=(CI_up)), width=.15) + #xlim(c(0,3))+
#   geom_vline(xintercept = 1) + xlab("Odds ratio") +
#   geom_point() + ylab('') + theme_bw() +
#   scale_y_discrete(breaks=names,
#                    labels=c("Promoter", "Enhancer", "Gene body")) +
#   theme(legend.title=element_blank())
# g
# 
# ggsave("../figures/figures/figure_6/OR_age_enhancer.png", g, device="png", width = 4, height = 2)
# 
# 
# #CpG Island
# names <- c("island", "shore", "shelf", "open_sea")
# hypo_d <- read_data(names, hypo)
# hyper_d <- read_data(names, hyper)
# hyper_hypo <- rbind(hypo_d, hyper_d)
# hyper_hypo$region <- factor(hyper_hypo$region, levels=rev(names))
# 
# g <- ggplot(hyper_hypo[hyper_hypo$adjPvalue<0.05,], aes(x=(oddsRatio), y=region, colour=type)) + 
#   geom_errorbar(aes(xmin=(CI_down), xmax=(CI_up)), width=.15) + #xlim(c(0,3))+
#   geom_vline(xintercept = 1) + xlab("Odds ratio") +
#   geom_point() + ylab('') + theme_bw() +
#   scale_y_discrete(breaks=names,
#                    labels=c("CpG island", "CpG shore", "CpG shelf", "Open sea")) +
#   theme(legend.title=element_blank())
# g
# 
# ggsave("../figures/figures/figure_6/OR_age_island.png", g, device="png", width = 4, height = 2)
# 
# 
# #genomic location
# names <- c("TSS200", "TSS1500", "3UTR", "5UTR", "1stExon", "Body", "intergenic")
# hypo_d <- read_data(names, hypo)
# hyper_d <- read_data(names, hyper)
# hyper_hypo <- rbind(hypo_d, hyper_d)
# hyper_hypo$region <- factor(hyper_hypo$region, levels=rev(names))
# 
# g <- ggplot(hyper_hypo[hyper_hypo$adjPvalue<0.05,], aes(x=(oddsRatio), y=region, colour=type)) + 
#   geom_errorbar(aes(xmin=(CI_down), xmax=(CI_up)), width=.15) + #xlim(c(0,3))+
#   geom_vline(xintercept = 1) + xlab("Odds ratio") +
#   geom_point() + ylab('') + theme_bw() +
#   scale_y_discrete(breaks=names,
#                    labels=c("TSS200", "TSS1500", "3'UTR", "5'UTR", "1st exon", "Gene body", "Intergenic")) +
#   theme(legend.title=element_blank())
# g
# 
# ggsave("../figures/figures/figure_6/OR_age_genomic.png", g, device="png", width = 4, height = 3)
# 


#Functional enrichments
library(missMethyl)
library(ggplot2)
plot_go <- function(go, direction){
  go <- go[go$ONTOLOGY=="BP",] #Only BP
  sig <- sum(go$FDR<0.05)
  print(sig)
  if(sig>20){
    max<-20
  } else{max<-sig}
  topgo <- topGSA(go, max)
  ggplot(data = topgo, aes(x=DE/N, y = factor(TERM, levels=rev(TERM)),
                           color = FDR, size = DE)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    theme_bw() + ylab("") + xlab("Gene Ratio") +
    ggtitle(paste0("GO ", direction, " BP"))
}

smoking <- rownames(results_DML$Smoking2[results_DML$Smoking2$adj.P.Val<0.05 & results_DML$Smoking2$logFC<0,])
age <- rownames(results_DML$AGE[results_DML$AGE$adj.P.Val<0.05 & results_DML$AGE$logFC<0,])
hypo_hypo <- smoking[smoking %in% age]

# go_hypo <- gometh(sig.cpg=hypo_hypo, all.cpg = rownames(signif), collection="GO", array.type="EPIC") #Output with genes takes a lot
go_hypo <- gometh(sig.cpg=hypo_hypo, all.cpg = annotation$IlmnID, collection="GO", array.type="EPIC")
sum(go_hypo$FDR<0.05)
plot_go(go_hypo, "hypo")
go_hypo <- go_hypo[go_hypo$FDR<0.05,]
go_hypo <- go_hypo[go_hypo$ONTOLOGY=="BP",] 

smoking <- rownames(results_DML$Smoking2[results_DML$Smoking2$adj.P.Val<0.05 & results_DML$Smoking2$logFC>0,])
age <- rownames(results_DML$AGE[results_DML$AGE$adj.P.Val<0.05 & results_DML$AGE$logFC>0,])
hyper_hyper <- smoking[smoking %in% age]
go_hyper <- gometh(sig.cpg=hyper_hyper, all.cpg = annotation$IlmnID, collection="GO", array.type="EPIC")
sum(go_hyper$FDR<0.05)
go_hyper <- go_hyper[go_hyper$FDR<0.05,]
go_hyper <- go_hyper[go_hyper$ONTOLOGY=="BP",] 
plot_go(go_hyper, "hyper")

go_hypo$direction <- "Hypomethylated"
go_hyper$direction <- "Hypermethylated"
to_plot_f <- rbind(go_hypo, go_hyper)
to_plot_f <- to_plot_f[order(to_plot_f$FDR),]

saveRDS(to_plot_f, "../figures/data/figure_6_enrichment_age.rds")
