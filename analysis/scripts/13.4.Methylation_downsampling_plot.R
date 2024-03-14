#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to downsample
# @software version: R=4.2.2


#Set path 
# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")

tissues <- c("BreastMammaryTissue", "ColonTransverse", "Lung", "MuscleSkeletal", "Ovary", "Prostate", "Testis", "WholeBlood")

#Parsing medians:

medians_hypo <- matrix(NA, length(tissues), 4)
rownames(medians_hypo) <- tissues
colnames(medians_hypo) <- c(11, 25, 35, 50)
medians_hyper <- matrix(NA, length(tissues), 4)
rownames(medians_hyper) <- tissues
colnames(medians_hyper) <- c(11, 25, 35, 50)

for(number in colnames(medians_hypo)){
# for(number in c(25, 35, 50)){
    print(number)
  print(Sys.time())
  for(tissue in tissues){
    print(tissue)
    print(Sys.time())
    if(!dir.exists(paste0("Downsampling/", number, "/", tissue, "/"))){
      next
    }
    files <- list.files(paste0("Downsampling/", number, "/", tissue, "/"), pattern = "results", full.names = T)
    files <- files[grep("methylation", files)]
    if(length(files)==0){next}
    hypo <- c()
    hyper <- c()
    for(file in files){
      # print(grep(file, files))
      data <- readRDS(file)
      hypo <- c(hypo, sum(data$Smoking2$logFC<0 & data$Smoking2$adj.P.Val<0.05))
      hyper <- c(hyper, sum(data$Smoking2$logFC>0 & data$Smoking2$adj.P.Val<0.05))
    }
    medians_hypo[tissue, as.character(number)] <- median(hypo)
    medians_hyper[tissue, as.character(number)] <- median(hyper)
    print(Sys.time())
  }
  print("Finished with:")
  print(number)
  Sys.time()
}
print(medians_hypo)
print(medians_hyper)

save.image("../figures/data/downsampling_methylation.Rdata")


#Plots are now in another script
# load("../figures/data/downsampling_methylation.Rdata")
# 
# tissue_info <- read.csv("data/public/tissue_abreviation.txt")
# rownames(medians_hypo) <- sapply(rownames(medians_hypo), function(name) tissue_info$Name[tissue_info$X==name])
# rownames(medians_hyper) <- sapply(rownames(medians_hyper), function(name) tissue_info$Name[tissue_info$X==name])
# 
# medians <- cbind(medians_hyper[,1], medians_hypo[,1], medians_hyper[,2], medians_hypo[,2], 
#                  medians_hyper[,3], medians_hypo[,3], medians_hyper[,4], medians_hypo[,4])
# medians <- round(medians)
# colnames(medians) <- c("22 hyper", "22 hypo", "50 hyper", "50 hypo", "70 hyper", "70 hypo", "100 hyper", "100 hypo")
# medians <- medians[c(3, 2, 1, 4, 5, 6, 7, 8),]
# 
# library(ComplexHeatmap)
# library(RColorBrewer)
# no_nas <- medians
# no_nas <- apply(apply(apply(no_nas, 2, as.numeric), 2, round), 2, prettyNum, big.mark=",")
# no_nas[no_nas=="NA"] <- ""
# 
# library(circlize)
# col_fun <- colorRamp2(c(0, 1, 5, 500, 7000, 20000, 25000), 
#                       c("white", brewer.pal(9,"BuPu")[c(1, 3:7)]))
# pdf("../figures/figures/figure_s9/downsampling_methylation.pdf", height = 2.4, width = 5.2)
# Heatmap(medians, col = col_fun,
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(no_nas[i, j], x, y, gp = gpar(fontsize = 9))},
#         na_col = "white",
#         cluster_rows = F,
#         column_names_side = "bottom",
#         cluster_columns = F,
#         name = "Number\nof DMPs",
#         column_title_side = "bottom",
#         column_title = "Number of samples",
#         row_names_gp = gpar(fontsize = 10),
#         column_names_gp = gpar(fontsize = 10),
#         row_names_side = "left",
#         heatmap_legend_param = list(
#           at=c(0, 12500, 25000)
#         ))
# dev.off()
#Sort tissues


#Option B:
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
# pdf("../figures/figures/figure_s9/downsampling_pretty.pdf", height = 2.5, width = 9.5)
# g
# dev.off()
