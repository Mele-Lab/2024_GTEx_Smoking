#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Enrichment on TFBS for differentially methylated positions/loci
# @software version: R=4.2.2

library(dplyr)
library(parallel)
library(ggplot2)

#Set path 
# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")

Sys.time()
data_path <- "~/Documents/mn4/scratch/bsc83/bsc83535/GTEx/v9/Oliva/"
annotation <- read.csv(paste0(data_path, "chip_seq_processed.csv")) 
Sys.time()
colnames(annotation)[8] <- c("TF")
annotation <- annotation[,-1]

length(unique(annotation$TF)) #221

results <- annotation[,c(4,7)]
results <- results %>% distinct()

#Are DMPs enriched in particular TFs?
tissue <- "Lung"
results_DML <- readRDS(paste0("tissues/", tissue , "/DML_results.rds"))
results_DML <- results_DML$Smoking2
smoking_hyper <- rownames(results_DML[results_DML$adj.P.Val<0.05 & results_DML$logFC>0,])
smoking_hypo <- rownames(results_DML[results_DML$adj.P.Val<0.05 & results_DML$logFC<0,])

fisher_function <- function(target, dmps){
  target_positions <- unique(results$name_ann[results$TF==target])
  non_target_positions <- rownames(results_DML)[!rownames(results_DML) %in% target_positions]
  non_dmp <- rownames(results_DML)[!rownames(results_DML) %in% dmps]
  m <- matrix(c(sum(dmps %in% target_positions),
                sum(dmps %in% non_target_positions),
                sum(non_dmp %in% target_positions),
                sum(non_dmp %in% non_target_positions)), nrow=2)
  t <- fisher.test(m)
  return(list(t$p.value, t$estimate, t$conf.int[1], t$conf.int[2]))
}

print("Computing fisher's for hypermethylated positions")
Sys.time()
hyper <- mclapply(unique(results$TF), function(target) fisher_function(target, smoking_hyper))
names(hyper) <- unique(results$TF)
hyper <- do.call(rbind, hyper)
hyper <- as.data.frame(hyper)
hyper[,1] <- unlist(hyper[,1])
hyper[,2] <- unlist(hyper[,2])
hyper[,3] <- unlist(hyper[,3])
hyper[,4] <- unlist(hyper[,4])
colnames(hyper) <- c("p_value", "odds_ratio", "CI_low", "CI_high")
hyper$adj_p_val <- p.adjust(hyper$p_value, "BH")
Sys.time()

saveRDS(hyper, paste0("tissues/", tissue , "/TFBS_ChipAtlas_hyper.rds"))
# hyper <- readRDS(paste0("tissues/", tissue , "/TFBS_ChipAtlas_hyper.rds"))

print("Computing fisher's for hypomethylated positions")
Sys.time()
hypo <- mclapply(unique(results$TF), function(target) fisher_function(target, smoking_hypo)) #15 minutes
names(hypo) <- unique(results$TF)
hypo <- do.call(rbind, hypo)
hypo <- as.data.frame(hypo)
colnames(hypo) <- c("p_value", "odds_ratio", "CI_low", "CI_high")
hypo[,1] <- unlist(hypo[,1])
hypo[,2] <- unlist(hypo[,2])
hypo[,3] <- unlist(hypo[,3])
hypo[,4] <- unlist(hypo[,4])
hypo$adj_p_val <- p.adjust(hypo$p_value, "BH")

Sys.time()

saveRDS(hypo, paste0("tissues/", tissue , "/TFBS_ChipAtlas_hypo.rds"))
# hypo <- readRDS(paste0("tissues/", tissue , "/TFBS_ChipAtlas_hypo.rds"))


#Interpreting results
table(hyper$adj_p_val<0.05 & hyper$odds_ratio>1)
table(hyper$adj_p_val<0.05 & hyper$odds_ratio<1)

table(hypo$adj_p_val<0.05 & hypo$odds_ratio>1)
table(hypo$adj_p_val<0.05 & hypo$odds_ratio<1)

#Plot top 15 enriched in hypo and hyper
# test <- hyper[hyper$adj_p_val<0.05 & hyper$odds_ratio>1,]
# hyper_plot <- test[order(test$odds_ratio, decreasing = T)[1:6],]
# hyper_plot$TF <- rownames(hyper_plot)
# hyper_plot$TF <- factor(hyper_plot$TF, levels = rev(hyper_plot$TF))
# 
# g2 <- ggplot(hyper_plot, aes(x=odds_ratio, y=TF, col="#CC6677")) + 
#   geom_errorbar(aes(xmin=CI_low, xmax=CI_high), width=.4) + 
#   geom_vline(xintercept = 1) + 
#   geom_point(size=2.5) + ylab('') + theme_bw() +
#   scale_colour_manual(values=c("#CC6677")) +
#   xlab("Odds ratio") +
#   theme(legend.position = "none",
#         axis.text.x = element_text(colour="black", size=11),
#         axis.text.y = element_text(colour="black", size=13),
#         legend.text = element_text(colour="black", size=12),
#         axis.title.x = element_text(size=13),
#         legend.spacing.y = unit(-0.05, "cm"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(colour = "black", linewidth=1))
# g2
# 
# test <- hypo[hypo$adj_p_val<0.05 & hypo$odds_ratio>1,]
# hyper_plot <- test[order(test$odds_ratio, decreasing = T)[1:15],]
# hyper_plot$TF <- rownames(hyper_plot)
# hyper_plot$TF <- factor(hyper_plot$TF, levels = rev(hyper_plot$TF))
# g3 <- ggplot(hyper_plot, aes(x=odds_ratio, y=TF, col="#88CCEE")) + 
#   geom_errorbar(aes(xmin=CI_low, xmax=CI_high), width=.4) + 
#   geom_vline(xintercept = 0) + 
#   geom_point(size=2.5) + ylab('') + theme_bw() +
#   scale_colour_manual(values=c("#88CCEE")) +
#   xlab("Odds ratio") +
#   theme(legend.position = "none",
#         axis.text.x = element_text(colour="black", size=11),
#         axis.text.y = element_text(colour="black", size=13),
#         legend.text = element_text(colour="black", size=12),
#         axis.title.x = element_text(size=13),
#         legend.spacing.y = unit(-0.05, "cm"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(colour = "black", linewidth=1))
# g3




#What about aging?

tissue <- "Lung"
results_DML <- readRDS(paste0("tissues/", tissue , "/DML_results.rds"))
results_DML <- results_DML$AGE
smoking_hyper <- rownames(results_DML[results_DML$adj.P.Val<0.05 & results_DML$logFC>0,])
smoking_hypo <- rownames(results_DML[results_DML$adj.P.Val<0.05 & results_DML$logFC<0,])


fisher_function <- function(target, dmps){
  target_positions <- unique(results$name_ann[results$TF==target])
  non_target_positions <- rownames(results_DML)[!rownames(results_DML) %in% target_positions]
  non_dmp <- rownames(results_DML)[!rownames(results_DML) %in% dmps]
  m <- matrix(c(sum(dmps %in% target_positions),
                sum(dmps %in% non_target_positions),
                sum(non_dmp %in% target_positions),
                sum(non_dmp %in% non_target_positions)), nrow=2)
  t <- fisher.test(m)
  return(list(t$p.value, t$estimate, t$conf.int[1], t$conf.int[2]))
}

print("Computing fisher's for hypermethylated positions")
Sys.time()
hyper <- mclapply(unique(results$TF), function(target) fisher_function(target, smoking_hyper))
names(hyper) <- unique(results$TF)
hyper <- do.call(rbind, hyper)
hyper <- as.data.frame(hyper)
hyper[,1] <- unlist(hyper[,1])
hyper[,2] <- unlist(hyper[,2])
hyper[,3] <- unlist(hyper[,3])
hyper[,4] <- unlist(hyper[,4])
colnames(hyper) <- c("p_value", "odds_ratio", "CI_low", "CI_high")
hyper$adj_p_val <- p.adjust(hyper$p_value, "BH")
Sys.time()

saveRDS(hyper, paste0("tissues/", tissue , "/TFBS_ChipAtlas_hyper_age.rds"))

print("Computing fisher's for hypomethylated positions")
Sys.time()
hypo <- mclapply(unique(results$TF), function(target) fisher_function(target, smoking_hypo)) 
names(hypo) <- unique(results$TF)
hypo <- do.call(rbind, hypo)
hypo <- as.data.frame(hypo)
colnames(hypo) <- c("p_value", "odds_ratio", "CI_low", "CI_high")
hypo[,1] <- unlist(hypo[,1])
hypo[,2] <- unlist(hypo[,2])
hypo[,3] <- unlist(hypo[,3])
hypo[,4] <- unlist(hypo[,4])
hypo$adj_p_val <- p.adjust(hypo$p_value, "BH")

Sys.time()

saveRDS(hypo, paste0("tissues/", tissue , "/TFBS_ChipAtlas_hypo_age.rds"))

table(hyper$adj_p_val<0.05 & hyper$odds_ratio>1)
table(hyper$adj_p_val<0.05 & hyper$odds_ratio<1)

table(hypo$adj_p_val<0.05 & hypo$odds_ratio>1)
table(hypo$adj_p_val<0.05 & hypo$odds_ratio<1)


#Do these significantly overlap those of smoking?
smoking_hypo <- readRDS(paste0("tissues/", tissue , "/TFBS_ChipAtlas_hypo.rds"))
smoking_hyper <- readRDS(paste0("tissues/", tissue , "/TFBS_ChipAtlas_hyper.rds"))
age_hypo <- readRDS(paste0("tissues/", tissue , "/TFBS_ChipAtlas_hypo_age.rds"))
age_hyper <- readRDS(paste0("tissues/", tissue , "/TFBS_ChipAtlas_hyper_age.rds"))

s_hyper <- rownames(smoking_hyper)[smoking_hyper$adj_p_val<0.05 & smoking_hyper$odds_ratio>1]
a_hyper <- rownames(age_hyper)[age_hyper$adj_p_val<0.05 & age_hyper$odds_ratio>1]
common <- sum(s_hyper %in% a_hyper)
only_s <- sum(!s_hyper %in% a_hyper)
only_a <- sum(!a_hyper %in% s_hyper)
none <- nrow(smoking_hyper) - common - only_a - only_s
fisher.test(matrix(c(common, only_s, only_a, none), nrow=2))

s_hyper <- rownames(smoking_hypo)[smoking_hypo$adj_p_val<0.05 & smoking_hypo$odds_ratio>1]
a_hyper <- rownames(age_hypo)[age_hypo$adj_p_val<0.05 & age_hypo$odds_ratio>1]
common <- sum(s_hyper %in% a_hyper)
only_s <- sum(!s_hyper %in% a_hyper)
only_a <- sum(!a_hyper %in% s_hyper)
none <- nrow(smoking_hypo) - common - only_a - only_s
fisher.test(matrix(c(common, only_s, only_a, none), nrow=2))


#Get table with the information we want:
function_to_parse <- function(data, variable){
  colnames(data) <- c("p_value", "Odds ratio", "CI lower bound", "CI upper bound", "Adjusted p-value" )
  data$"Transcrition factor" <- rownames(data)
  data$Direction <- "Enriched"
  data$Direction[data$`Odds ratio`<1] <- "Depleted"
  data$Trait <- variable
  data <- data[,c(8, 6, 5, 7, 2, 3, 4)]
  data <- data[order(data$Direction, data$`Adjusted p-value`, method= "radix", decreasing = c(T, F)),]
  data <- data[data$`Adjusted p-value`<0.05,]
  return(data)
}

test <- rbind(function_to_parse(smoking_hyper, "Smoking"), 
              function_to_parse(smoking_hypo, "Smoking"),
              function_to_parse(age_hyper, "Age"),
              function_to_parse(age_hypo, "Age"))

library("xlsx")
write.xlsx(test, "output/Supplementary_table_15.xlsx",
           col.names = TRUE, row.names = F, append = FALSE)


#Plot aging
# test <- hyper[hyper$adj_p_val<0.05,]
# hyper_plot <- test[order(test$odds_ratio, decreasing = T)[1:15],]
# hyper_plot$TF <- rownames(hyper_plot)
# hyper_plot$TF <- factor(hyper_plot$TF, levels = rev(hyper_plot$TF))
# 
# g2 <- ggplot(hyper_plot, aes(x=odds_ratio, y=TF, col="#CC6677")) + 
#   geom_errorbar(aes(xmin=CI_low, xmax=CI_high), width=.4) + 
#   geom_vline(xintercept = 1) + 
#   geom_point(size=2.5) + ylab('') + theme_bw() +
#   scale_colour_manual(values=c("#CC6677")) +
#   xlab("Odds ratio") +
#   theme(legend.position = "none",
#         axis.text.x = element_text(colour="black", size=11),
#         axis.text.y = element_text(colour="black", size=13),
#         legend.text = element_text(colour="black", size=12),
#         axis.title.x = element_text(size=13),
#         legend.spacing.y = unit(-0.05, "cm"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(colour = "black", linewidth=1))
# g2
# 
# test <- hypo[hypo$adj_p_val<0.05,]
# hyper_plot <- test[order(test$odds_ratio, decreasing = T)[1:15],]
# hyper_plot$TF <- rownames(hyper_plot)
# hyper_plot$TF <- factor(hyper_plot$TF, levels = rev(hyper_plot$TF))
# g3 <- ggplot(hyper_plot, aes(x=odds_ratio, y=TF, col="#88CCEE")) + 
#   geom_errorbar(aes(xmin=CI_low, xmax=CI_high), width=.4) + 
#   geom_vline(xintercept = 0) + 
#   geom_point(size=2.5) + ylab('') + theme_bw() +
#   scale_colour_manual(values=c("#88CCEE")) +
#   xlab("Odds ratio") +
#   theme(legend.position = "none",
#         axis.text.x = element_text(colour="black", size=11),
#         axis.text.y = element_text(colour="black", size=13),
#         legend.text = element_text(colour="black", size=12),
#         axis.title.x = element_text(size=13),
#         legend.spacing.y = unit(-0.05, "cm"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(colour = "black", linewidth=1))
# g3
