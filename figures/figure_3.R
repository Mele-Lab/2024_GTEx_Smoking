# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Generate publication Figure 2 (Alternative splicing and demographic traits)


#Set path 
# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio

#Load libraries 
library(ggplot2) 
library(RColorBrewer) 
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
library(reshape2)
library(ggrepel)
library(forcats)
library(colorspace) #To make lighter colors


#Load data
figure_data <- readRDS(file = "data/splicing_data.rds")
for(i in 1:length(figure_data)) assign(names(figure_data)[i], figure_data[[i]])
figure_data <- readRDS(file = "data/splicing_data2.rds")
for(i in 1:length(figure_data)) assign(names(figure_data)[i], figure_data[[i]])

source("safe_colourblind_pallete.R")
tissue_info <- read.csv("tissue_abreviation.txt")


#Additive effects plots

tissue_order <- readRDS(file = "data/degs_summary.rds")$DEGS
tissue_order <- sort(table(tissue_order$tissue))

#Plot heatmap of additive effects between smoking and the demographic traits
create_heatmap <- function(data, tissue_info, tissues_plot, diseases_plot, size=12, name="#DEGs", label_data=data, color_palette, row_names=T){ #Function for 7A and 7B
  data_plot <- data[tissues_plot, diseases_plot]
  without_NA <- replace(label_data, is.na(label_data), "")
 
  data_plot <- data_plot[match(rev(names(tissue_order)), rownames(data_plot)),]
  without_NA <- without_NA[match(rev(names(tissue_order)), rownames(without_NA)),]
  
  rownames(data_plot) <- sapply(rownames(data_plot), function(tissue) tissue_info$Name[tissue_info$tissue==tissue])
  traits_cols <- c("Age" = "#56B4E9", "Ancestry" = "#E69F00", "Sex" = "#009E73", "BMI" = "#CC79A7", "Disease" = "#000000")
  column_ha_top <- HeatmapAnnotation("Number of\nsignificant\noverlaps" = anno_barplot(colSums(to_plot_overlap>1),
                                                                                  border = F,
                                                                                  axis_param = list(gp=gpar(fontsize = 9)),
                                                                                  gp = gpar(fill = traits_cols,
                                                                                            col = traits_cols)),
                                     show_annotation_name = T,
                                     annotation_name_side = "left",
                                     annotation_name_rot = 90,
                                     annotation_name_gp = gpar(fontsize = 11),
                                     height = unit(1.8, "cm"))
  colnames(data_plot) <- c("Age + Smoking", "Ancestry + Smoking", "Sex + Smoking", "BMI + Smoking")
  ht <- Heatmap(data_plot,
          heatmap_legend_param = list(legend_height = unit(3, "cm"),
                                      row_names_gp = gpar(fontface = "bold"),
                                      grid_width = unit(0.5, "cm"),
                                      labels_gp=gpar(fontsize=size),
                                      title_position="topcenter",
                                      title_gp=gpar(fontsize=size, fontface=2)),#,
                                      #direction="horizontal"),
          col = color_palette, 
          top_annotation = column_ha_top,
          na_col = "white",
          cluster_rows = F,
          cluster_columns = F,
          name = name,
          row_names_side = "left",
          column_names_side = "bottom",
          column_names_rot =  60,
          column_names_gp = gpar(fontsize = size),
          column_names_max_height= unit(9, "cm"),
          row_names_gp = gpar(fontsize = size),
          # left_annotation = row_ha_left,
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(my_pretty_num_function(without_NA[i, j]), x, y, gp = gpar(fontsize = size))}
          
  )
  draw(ht)#, heatmap_legend_side="bottom")
}

my_pretty_num_function <- function(n){
  if(n==""){
    return(n)
  } else{
    prettyNum(n, big.mark = ",")
  }
}
# colours <- gray.colors(200, start = 0, end = 1)[c(200, 160:30)]

colours <- c("white", brewer.pal(9,"BuPu")[c(4,5,7)])
colours <- colorRamp2(c(1,2,3,4), colours)
to_plot_overlap <- to_plot_overlap[nrow(to_plot_overlap):1,]
overlap <- overlap[nrow(overlap):1,]

pdf("figures/figure_3/Overlaps.pdf", width = 4.2, height = 8) #if png -> units="in", res=200
create_heatmap(data = to_plot_overlap, tissue_info, tissues_plot = rownames(overlap), 
               diseases_plot = colnames(overlap), name="Odds ratio", label_data=overlap, 
               color_palette = colours)
dev.off()


#Bias images

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

pdf("figures/figure_3/Bias_age_overlap.pdf", width = 6, height = 2.5) #maybe 7?
# plot_bias("Age", c("old - smokers", "young - smokers", "old - never smokers", "young - never smokers"))
plot_bias("Age")
dev.off()



#Figure 2 B and D
get_box_stats <- function(y, upper_limit) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "N =", length(y), "\n"
    )
  ))
}
# p1 <- my.box_plot(data_PDGFRB, gene, "Smoking", "Age2")
my.box_plot <- function(data, gene, variable_col, x_variable){
  cols <- c("#88CCEE",  "#CC6677")
  facet_color <- "#805D93"
  if(x_variable == "Age_int"){ x_facet <- "Age"
  }else if(x_variable == "Age2"){ 
    x_facet <- "Age2"
    facet_color <- traits_cols["Age"]
  }else if(x_variable == "BMI_int"){ 
    x_facet <- "BMI"
    facet_color <- traits_cols["BMI"]
  }else if(x_variable == "BMI2"){ x_facet <- "BMI2"
  }else{x_facet <- x_variable}
  gene.name <- gene_annotation[gene_annotation$gene == gene, "symbol"]
  p <- ggplot(data = data,
              aes(x = x_dummy,
                  y = log2(TPM+1))
  ) +
    geom_violin(aes(fill = eval(parse(text=variable_col))),
                col = "black") +
    geom_boxplot(col = "black",
                 outlier.shape = NA,
                 notch = T,
                 width = 0.25) +
    geom_jitter(col = "black",
                alpha = 0.1,
                size = 0.8) +
    theme_minimal() +
    xlab("") +
    ylab("log2(TPM+1)") +
    scale_fill_manual(values = cols) +
    labs(title=paste0(gene.name,
                      " (",
                      tissue,
                      ")")) +
    theme(legend.position = "none",
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          plot.title = element_text(hjust = 0.5,
                                    size = 20)) +
    stat_summary(aes(x = x_dummy, y = log2(TPM)),
                 fun.data = get_box_stats, fun.args = list(upper_limit = max(log2(data$TPM + 1)) * 1.15), 
                 geom = "text", hjust = 0.5, vjust = 0.9, size = 5) +
    facet_grid(.~ eval(parse(text = x_facet)), drop = T, scales = "free") +
    theme(axis.title.y = element_text(margin = margin(r = 20), size = 18),
          axis.text.y = element_text(size = 16, colour = "black"), 
          axis.text.x = element_text(size = 16, colour = "black"),
          panel.background = element_rect(fill = "white"),
          axis.line = element_line(colour = "grey"),
          legend.position = "none", 
          plot.title = element_text(hjust = 0.5, size = 20),
          strip.background = element_rect(fill="#88CCEE"), 
          strip.text = element_text(size = 16)) +
    scale_x_discrete(labels=c('Never smokers', 'Smokers', 'Never smokers', 'Smokers'))
  return(p)
}


#Plot boxplots
traits_cols <- c("Age" = "#56B4E9", "Ancestry" = "#E69F00", 
                 "Sex" = "#009E73", "BMI" = "#CC79A7", "Smoking" = "#000000") #696969
my_traits <- names(traits_cols)


gene <- "ENSG00000113721.13" #PDGFRB
tissue <- "Thyroid"
p1 <- my.box_plot(data_PDGFRB, gene, "Smoking", "Age2")
pdf(paste0("figures/figure_3/", tissue, "_", gene_annotation[gene_annotation$gene == gene, "symbol"], ".pdf"),
    width = 6.5, height = 3.5)
p1
dev.off()

#I may want to find another lung example with age, maybe a known one
gene <- "ENSG00000244694.7" #PTCHD4
tissue <- "Lung"
p2 <- my.box_plot(data_PTCHD4, gene, "Smoking", "Age2")
pdf(paste0("figures/figure_3/", tissue, "_", gene_annotation[gene_annotation$gene == gene, "symbol"], ".pdf"),
    width = 6.7, height = 3.5)
p2
dev.off()
