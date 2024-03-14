#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to downsample
# @software version: R=4.2.2


#Set path 
setwd(system("pwd", intern = T)) #If in linux
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
# setwd("..")


#Plot results

#This next part may take 2 hours
tissues <- list.dirs("Downsampling/11/", full.names = F)[-1]
medians <- matrix(NA, length(tissues), 6)
rownames(medians) <- tissues
colnames(medians) <- c(11, 25, 35, 50, 80, 100)

for(number in c(11, 25, 35, 50, 80, 100)){
# for(number in colnames(medians)){
  print(number)
  print(Sys.time())
  for(tissue in tissues){
    print(tissue)
    print(Sys.time())
    if(!dir.exists(paste0("Downsampling/", number, "/", tissue, "/"))){
      next
    }
    files <- list.files(paste0("Downsampling/", number, "/", tissue, "/"), pattern = "results", full.names = T)
    files <- files[!grepl("methylation", files)]
    vector <- c()
    for(file in files){
      data <- readRDS(file)
      vector <- c(vector, sum(data$Smoking2$adj.P.Val<0.05))
    }
    medians[tissue, as.character(number)] <- median(vector)

    print(Sys.time())
  }
  print("Finished with:")
  print(number)
  Sys.time()
}
print(medians)
save.image("Downsampling/downsampling_expression.Rdata")
# 
# 
# # 
# # 
# # to_plot <- matrix(1,46,9)
# # 
# # for(i in 1:length(metadata)){
# #   tissue <- names(metadata)[i]
# #   print(tissue)
# #   dea_res <- readRDS(paste0("tissues/", tissue, "/voom_limma_results_downsampling_2.rds"))
# #   to_plot[i,1] <- sum(dea_res$Smoking2$adj.P.Val<0.05)
# #   dea_res <- readRDS(paste0("tissues/", tissue, "/voom_limma_results_downsampling_3.rds"))
# #   to_plot[i,2] <- sum(dea_res$Smoking2$adj.P.Val<0.05)
# #   dea_res <- readRDS(paste0("tissues/", tissue, "/voom_limma_results_downsampling_4.rds"))
# #   to_plot[i,3] <- sum(dea_res$Smoking2$adj.P.Val<0.05)
# #   dea_res <- readRDS(paste0("tissues/", tissue, "/voom_limma_results_downsampling_5.rds"))
# #   to_plot[i,4] <- sum(dea_res$Smoking2$adj.P.Val<0.05)
# #   dea_res <- readRDS(paste0("tissues/", tissue, "/voom_limma_results_downsampling_6.rds"))
# #   to_plot[i,5] <- sum(dea_res$Smoking2$adj.P.Val<0.05)
# #   dea_res <- readRDS(paste0("tissues/", tissue, "/voom_limma_results_downsampling_7.rds"))
# #   to_plot[i,6] <- sum(dea_res$Smoking2$adj.P.Val<0.05)
# #   dea_res <- readRDS(paste0("tissues/", tissue, "/voom_limma_results_downsampling.rds"))
# #   to_plot[i,7] <- sum(dea_res$Smoking2$adj.P.Val<0.05)
# #   dea_res <- readRDS(paste0("tissues/", tissue, "/voom_limma_results_downsampling_50.rds"))
# #   to_plot[i,8] <- sum(dea_res$Smoking2$adj.P.Val<0.05)
# #   dea_res <- readRDS(paste0("tissues/", tissue, "/voom_limma_results_downsampling_80.rds"))
# #   to_plot[i,9] <- sum(dea_res$Smoking2$adj.P.Val<0.05)
# #   
# # }
# # rownames(to_plot) <- names(metadata)
# # 
# # 
# # #plot
# # library(ComplexHeatmap)
# # library(RColorBrewer)
# # Heatmap(apply(to_plot, 2,function(x) x/max(x,na.rm=T)), col = brewer.pal(9,"BuPu")[1:7],
# #         cell_fun = function(j, i, x, y, width, height, fill) {
# #           grid.text(prettyNum(to_plot[i, j], big.mark = ","), x, y, gp = gpar(fontsize = 10))},
# #         na_col = "white",
# #         cluster_rows = F,
# #         cluster_columns = F,
# #         name = "DE signal",
# #         row_names_side = "left",)
# 
# 
# #To do after job on this script has finished:
# load("Downsampling/downsampling.Rdata")
# 
# tissue_info <- read.csv("data/public/tissue_abreviation.txt")
# rownames(medians) <- sapply(rownames(medians), function(name) tissue_info$Name[tissue_info$X==name])
# 
# library(ComplexHeatmap)
# library(RColorBrewer)
# medians <- medians[rev(order(medians[,1])),]
# no_nas <- medians
# no_nas <- apply(apply(apply(no_nas, 2, as.numeric), 2, round), 2, prettyNum, big.mark=",")
# no_nas[no_nas=="NA"] <- ""
# 
# pdf("../figures/figures/figure_s1/downsampling.pdf", height = 6.5, width = 4.3)
# Heatmap(apply(medians, 2,function(x) x/max(x,na.rm=T)), col = brewer.pal(9,"BuPu")[c(1,4:7)],
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(no_nas[i, j], x, y, gp = gpar(fontsize = 9))},
#         na_col = "white",
#         cluster_rows = F,
#         column_names_side = "bottom",
#         cluster_columns = F,
#         name = "Differential\nexpression\nsignal",
#         column_names_rot = 45,
#         column_title_side = "bottom",
#         column_title = "Number of samples",
#         row_names_gp = gpar(fontsize = 10),
#         column_names_gp = gpar(fontsize = 10),
#         row_names_side = "left",)
# dev.off()
# #Sort tissues
# 
# 
# #Option B:
# library(reshape2)
# library(ggplot2)
# 
# df <- melt(medians)
# df$Var2 <- as.factor(df$Var2)
# cols <- tissue_info$color
# names(cols) <- tissue_info$Name
# 
# g <- ggplot(df, aes(Var2, value, col=Var1)) + geom_point() + 
#   scale_color_manual(values=cols) +
#   geom_line(aes(group = Var1)) + theme_bw() +
#   theme(legend.title = element_blank(),
#         axis.text = element_text(colour="black", size=14),
#         legend.text = element_text(colour="black", size=13),
#         axis.title = element_text(size=16),
#         # legend.spacing.y = unit(-0.05, "cm"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(colour = "black", linewidth=1)) +
#   ylab("Number of DEGs") + xlab("Number of samples")
# 
# pdf("../figures/figures/figure_s1/downsampling_pretty.pdf", height = 2.5, width = 9.5)
# g
# dev.off()
