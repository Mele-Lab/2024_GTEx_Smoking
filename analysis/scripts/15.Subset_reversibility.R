#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Methylation results analysis on reversibility using the same subset as expression
# @software version: R=4.2.2

#Set path 
# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")

# print("Reading results")
# results <- readRDS("tissues/Lung/DML_results_subset.rds")
# 
# smoking2 <- results$Smoking2
# signif <- smoking2[smoking2$adj.P.Val<0.05,]
# table(signif$logFC>0)
# smoking2 <- rownames(smoking2[smoking2$adj.P.Val<0.05,])
# smoking1 <- results$Smoking1
# smoking1 <- rownames(smoking1[smoking1$adj.P.Val<0.05,])
# smoking2_1 <- results$`Smoking1-Smoking2`
# smoking2_1 <- rownames(smoking2_1[smoking2_1$adj.P.Val<0.05,])
# 
# #How mamy reversible or nonr reversible positions we find?
# non_rev <- smoking1[smoking1 %in% smoking2]
# length(non_rev)
# rev <- smoking2_1[smoking2_1 %in% smoking2]
# length(rev)
# 
# #Are they hypo- or hyper-methylated?
# table(signif$logFC[rownames(signif) %in% non_rev]>0)/length(non_rev) #36% hypo
# 
# partially_reversible <- smoking2
# partially_reversible <- partially_reversible[!partially_reversible %in% smoking2_1]
# partially_reversible <- partially_reversible[!partially_reversible %in% smoking1]
# length(partially_reversible)/length(smoking2)



#Are ex-smoker betas of PR positions closer to smokers or ex smokers levels?
# Sys.time()
# betas <- readRDS("tissues/Lung/methylation_residuals.rds")
# Sys.time() 
# 
# #Subset of PR
# dim(betas)
# betas <- betas[rownames(betas) %in% partially_reversible,]
# dim(betas)
# 
# metadata <- readRDS("tissues/Lung/methylation_metadata_subset.rds")[,c(1,36)] 
# smokers <- betas[,colnames(betas) %in% metadata$SUBJID[metadata$Smoking==2]]
# never_smokers <- betas[,colnames(betas) %in% metadata$SUBJID[metadata$Smoking==0]]
# ex_smokers <- betas[,colnames(betas) %in% metadata$SUBJID[metadata$Smoking==1]]
# 
# #Are ex smokers closer to smoker or never smoker levels?
# table(abs(rowMeans(ex_smokers) - rowMeans(smokers)) > abs(rowMeans(ex_smokers) - rowMeans(never_smokers)))
# table(abs(rowMeans(ex_smokers) - rowMeans(smokers)) > abs(rowMeans(ex_smokers) - rowMeans(never_smokers)))/nrow(ex_smokers)
# #60% of ex-smokers are closer to smokers
# binom.test(18451, 30474, 0.5) #It is significant



#LogFC plot:
results <- readRDS("tissues/Lung/DML_results_subset.rds")
to_exclude <- rownames(results$Smoking2)[results$Smoking2$adj.P.Val>=0.05] #Smoking-DEGs.
never_vs_ex <- results$Smoking1
# to_exclude <- c(to_exclude, rownames(never_vs_ex)[never_vs_ex$adj.P.Val<0.05])
ex_vs_smokers <- results$`Smoking1-Smoking2`
# to_exclude <- c(to_exclude, rownames(ex_vs_smokers)[ex_vs_smokers$adj.P.Val<0.05])
never_vs_ex <- never_vs_ex[!rownames(never_vs_ex) %in% to_exclude,]
ex_vs_smokers <- ex_vs_smokers[!rownames(ex_vs_smokers) %in% to_exclude,]

boxplot <- data.frame(abs(never_vs_ex$logFC), "Never smokers\nvs\nex-smokers")
names(boxplot) <- c("logFC", "Variable")
second_boxplot <- data.frame(abs(ex_vs_smokers$logFC), "Ex-smokers\nvs\nsmokers")
names(second_boxplot) <- c("logFC", "Variable")

# combined <- cbind(boxplot, second_boxplot)
# table(combined[,1]> combined[,3])/nrow(combined)

boxplot <- rbind(boxplot, second_boxplot)
test <- wilcox.test(abs(never_vs_ex$logFC), abs(ex_vs_smokers$logFC), paired=T)
# test <- wilcox.test(combined[,1], combined[,3], paired=T)
#Never vs ex is less than ex vs smokers



#What about expression results?
results <- readRDS("tissues/Lung/voom_limma_results_subset.rds")
smoking <- rownames(results$Smoking2)[results$Smoking2$adj.P.Val<0.05]
sum(rownames(results$Smoking1)[results$Smoking1$adj.P.Val<0.05] %in% smoking)
sum(rownames(results$`Smoking1-Smoking2`)[results$`Smoking1-Smoking2`$adj.P.Val<0.05] %in% smoking)

#Are ex smokers closer to smoker or never smoker levels?
# tpm <- readRDS("tissues/Lung/tpm.rds")
# #Subset of PR
# partially_reversible <- smoking
# partially_reversible <- partially_reversible[!partially_reversible %in% rownames(results$Smoking1)[results$Smoking1$adj.P.Val<0.05]]
# partially_reversible <- partially_reversible[!partially_reversible %in% rownames(results$`Smoking1-Smoking2`)[results$`Smoking1-Smoking2`$adj.P.Val<0.05]]
# length(partially_reversible)/length(smoking)
# 
# dim(tpm)
# tpm <- tpm[rownames(tpm) %in% partially_reversible,]
# colnames(tpm) <- sapply(colnames(tpm), function(name) paste(strsplit(name, "-")[[1]][-3], collapse="-"))
# dim(tpm)
# 
# smokers <- tpm[,colnames(tpm) %in% metadata$SUBJID[metadata$Smoking==2]]
# never_smokers <- tpm[,colnames(tpm) %in% metadata$SUBJID[metadata$Smoking==0]]
# ex_smokers <- tpm[,colnames(tpm) %in% metadata$SUBJID[metadata$Smoking==1]]
# 
# table(abs(rowMeans(ex_smokers) - rowMeans(smokers)) > abs(rowMeans(ex_smokers) - rowMeans(never_smokers)))
# table(abs(rowMeans(ex_smokers) - rowMeans(smokers)) > abs(rowMeans(ex_smokers) - rowMeans(never_smokers)))/nrow(ex_smokers)
#70% of ex-smokers are closer to never smokers


to_exclude <- rownames(results$Smoking2)[results$Smoking2$adj.P.Val>=0.05] #these will not be partially reversible
never_vs_ex <- results$Smoking1
# to_exclude <- c(to_exclude, rownames(never_vs_ex)[never_vs_ex$adj.P.Val<0.05])
ex_vs_smokers <- results$`Smoking1-Smoking2`
# to_exclude <- c(to_exclude, rownames(ex_vs_smokers)[ex_vs_smokers$adj.P.Val<0.05])
never_vs_ex <- never_vs_ex[!rownames(never_vs_ex) %in% to_exclude,]
ex_vs_smokers <- ex_vs_smokers[!rownames(ex_vs_smokers) %in% to_exclude,]

boxplot_e <- data.frame(abs(never_vs_ex$logFC), "Never smokers\nvs\nex-smokers")
names(boxplot_e) <- c("logFC", "Variable")
second_boxplot <- data.frame(abs(ex_vs_smokers$logFC), "Ex-smokers\nvs\nsmokers")
names(second_boxplot) <- c("logFC", "Variable")
# 
# combined <- cbind(boxplot, second_boxplot)
# table(combined[,1]> combined[,3])/nrow(combined)
# 
boxplot_e <- rbind(boxplot_e, second_boxplot)
test_e <- wilcox.test(abs(never_vs_ex$logFC), abs(ex_vs_smokers$logFC), paired=T)


figure_S10 <- list(boxplot, test, boxplot_e, test_e)
names(figure_S10) <- c("boxplot_subset", "test_subset", "boxplot_expression", "test_expression")
saveRDS(figure_S10, paste0("../figures/data/figure_S8_logFC_subset.rds"))
