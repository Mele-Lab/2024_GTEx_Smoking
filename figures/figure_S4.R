#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Date:   2022-08-18
# @Description: Generate publication Figure s7 (Hierarchical partitioning on smoking and demographic traits)


#Set path 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Load libraries 
suppressPackageStartupMessages(library(yarrr))
library(reshape2)
library(scales)
library(RColorBrewer) 
suppressPackageStartupMessages(library(ComplexHeatmap))

#Load data
figure_data <- readRDS(file = "data/figure_S3.rds") #list of expression_barplot, expression_pirate, splicing_barplot and splicing_pirate
for(i in 1:length(figure_data)) assign(names(figure_data)[i], figure_data[[i]])

figure_data <- readRDS(file = "data/splicing_data.rds")
for(i in 1:length(figure_data)) assign(names(figure_data)[i], figure_data[[i]])
figure_data <- readRDS(file = "data/splicing_data2.rds")
for(i in 1:length(figure_data)) assign(names(figure_data)[i], figure_data[[i]])

## Figure 2A-C

# plot_bias <- function(variable, col_names){
plot_bias <- function(variable){
  data1 <- as.matrix(sapply(tissues, function(tissue) chi_square.results[[variable]][[tissue]][["observed"]]))
  names <- rownames(data1)
  data1 <- cbind(sapply(strsplit(data1,'/'), "[", 1), 
                 sapply(strsplit(data1,'/'), "[", 2),
                 sapply(strsplit(data1,'/'), "[", 3),
                 sapply(strsplit(data1,'/'), "[", 4))
  rownames(data1) <- names
  data1 <- as.data.frame(data1)
  data1[c(1:4)] <- sapply(data1[c(1:4)],as.numeric) #I add c(1:4) to keep rownames
  data1 <- t(data1)
  col_names <- c("up - up", "down - up", "up - down", "down - down")
  rownames(data1) <- col_names
  
  #Only plot bias when significant overlap
  tissues_to_plot <- to_plot_overlap[,variable]>1 #Only signif overlap
  data1 <- data1[,colnames(data1) %in% names(tissues_to_plot[tissues_to_plot])] #Only signif overlap
  if(is.null(dim(data1))){
    if(sum(data1)>10){print("correct")}
  }else{
    data1 <- data1[,colSums(data1)>10]
  }
  
  p_vals <- p_values[,colnames(p_values)==variable]
  if(is.null(dim(data1))){
    p_vals <- p_vals[names(tissues_to_plot[tissues_to_plot])]
  }else{
    p_vals <- p_vals[names(p_vals) %in% colnames(data1)] #Only signif overlap
  }
  
  p_vals <- p.adjust(p_vals, method = "BH")
  p_vals <- p_vals[p_vals<0.05]
  
  if(length(p_vals)==1){ #Only one tissue under study, only 1 column, so we need to transform the vector into data frame
    if(ncol(data1)>1){
      data1 <- data1[,colnames(data1) %in% names(p_vals)] #This would not work with only one column
    }
    data1 <- as.data.frame(data1)
    colnames(data1) <- names(p_vals)
  } else{
    my_col_names <- names(p_vals)
    data1 <- data1[,colnames(data1) %in% my_col_names] #This would not work with only one column
    data1 <- as.data.frame(data1)
    if(ncol(data1)==1){
      colnames(data1) <- my_col_names[!is.na(my_col_names)]
    }
  }
  
  
  data1$type <- rownames(data1)
  data3 <- as.data.frame(t(data1)[1:ncol(data1)-1,])
  data3[c(1:ncol(data3))] <- sapply(data3[c(1:ncol(data3))],as.numeric) #I add [c(1:4)] to keep rownames
  if(ncol(data1)==2){
    data3 <- data3/colSums(data3)
    colnames(data3) <- colnames(data1)[-2]
  }else{
    data3 <- data3/rowSums(data3)
  }
  data3 <- as.data.frame(t(data3))
  data3$type <- rownames(data3)
  
  data <- melt(data1)
  data$y <- 100*melt(data3)[,3]
  data$type <- factor(data$type, levels = rev(col_names), order = T)
  #Plot the tissue acronyms instead of the tissue whole names
  data$variable <- sapply(data$variable, function(tissue) tissue_info$Name[tissue_info$tissue==tissue])
  #Order by FDR
  # sorting <- rownames(fdr_bias)[order(fdr_bias[,1])]
  sorting <- names(p_vals)[order(p_vals)]
  sorting <- sapply(sorting, function(tissue) tissue_info$Name[tissue_info$tissue==tissue]) 
  data$variable <- factor(data$variable, levels=rev(sorting), order=T)
  data$variable <- droplevels(data$variable)
  data$type <- factor(data$type, levels=rev(c("up - up", "down - down", "down - up", "up - down")))
  if(variable=="Ancestry"){variable <- "African"}
  else if(variable=="Sex"){variable <- "Female"}
  else if(variable=="BMI"){variable <- "High BMI"}
  
  plot <- ggplot(data,
                 aes(x = variable, 
                     y = y, 
                     fill = type,
                     label = value) ) +
    geom_bar(stat = "identity") +
    coord_flip() +
    geom_text(size = 4.5, position = position_stack(vjust = 0.5)) +
    scale_fill_manual(values = rev(c(lighten("#88CCEE", amount = 0.5), "#88CCEE", 
                                     lighten("#CC6677", amount = 0.3), "#CC6677")), 
                      breaks = rev(levels(data$type))) +
    theme_bw() +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 14.5, colour="black"),
          axis.title = element_text(size = 14.5, colour="black"),
          legend.title = element_text(size = 14.5, colour="black"),
          axis.text.x = element_text(vjust = 0.5, size = 13),
          legend.key.size = unit(1, 'cm'),
          legend.text = element_text(size=13.5, colour="black")) +
    ylab("Differentially expressed genes (%)") + xlab("") + 
    labs(fill=paste(variable, "- Smoking")) 
  return(plot)
}



pdf("figures/figure_s3_4/Bias_ancestry.pdf", width = 7, height = 1.5)
# plot_bias("Ancestry", c("African - smokers", "European - smokers", "African - never smokers", "European - never smokers"))
plot_bias("Ancestry")
dev.off()

pdf("figures/figure_s3_4/Bias_BMI.pdf", width = 5.5, height = 1)
# plot_bias("BMI", c("High BMI - smokers", "Low BMI - smokers", "High BMI - never smokers", "Low BMI - never smokers"))
plot_bias("BMI")
dev.off()
pdf("figures/figure_s3_4/Bias_sex.pdf", width = 6.5, height = 1)
# plot_bias("Sex", c("Female - smokers", "Male - smokers", "Female - never smokers", "Male - never smokers"))
plot_bias("Sex")
dev.off()


## Figure 3D

#order:
figure_data <- readRDS(file = "data/degs_summary.rds")
order <- figure_data$DEGS.per.tissue
order <- sort(order[,3], decreasing = T)
order_vector <- match(names(order), rownames(interactions))
interactions <- interactions[order_vector[!is.na(order_vector)],]
#Figure 3D

tissue_info <- read.csv("tissue_abreviation.txt")

create_heatmap <- function(data, tissue_info, tissues_plot, diseases_plot, size=12, name="#DEGs", label_data=data, color_palette, row_names=T){ #Function for 7A and 7B
  data_plot <- data[tissues_plot, diseases_plot]
  without_NA <- replace(label_data, is.na(label_data), "")
  
  rownames(data_plot) <- sapply(rownames(data_plot), function(tissue) tissue_info$Name[tissue_info$tissue==tissue])

  if(row_names==F){
    rownames(data_plot) <- NULL
  }
  ht <- Heatmap(data_plot,
                heatmap_legend_param = list(legend_height = unit(5, "cm"),
                                            row_names_gp = gpar(fontface = "bold"),
                                            grid_width = unit(0.75, "cm"),
                                            labels_gp=gpar(fontsize=size),
                                            title_position="topcenter",
                                            title_gp=gpar(fontsize=size, fontface=2),
                                            direction="vertical"),
                col = color_palette, 
                na_col = "white",
                cluster_rows = F,
                cluster_columns = F,
                name = name,
                row_names_side = "left",
                column_names_side = "top",
                column_names_rot =  60,
                column_names_gp = gpar(fontsize = size),
                column_names_max_height= unit(9, "cm"),
                row_names_gp = gpar(fontsize = size),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(my_pretty_num_function(without_NA[i, j]), x, y, gp = gpar(fontsize = size))}
                
  )
  draw(ht, heatmap_legend_side="right")
}

my_pretty_num_function <- function(n){
  if(n==""){
    return(n)
  } else{
    prettyNum(n, big.mark = ",")
  }
}
colours <- brewer.pal(9,"BuPu")[c(1,2,3)]
pdf("figures/figure_s3/Interactions.pdf", width = 4, height = 8)
create_heatmap(interactions, tissue_info, tissues_plot = rownames(interactions), diseases_plot = colnames(interactions),
               label_data=interactions, color_palette = colours, name="Number\nof interactions")
dev.off()
