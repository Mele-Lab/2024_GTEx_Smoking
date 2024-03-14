#Code to explore the DMPs associated with AHRR


#Set path
# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")

pairs <- readRDS("output/gene_probe_pairs.rds")
pairs <- pairs[!is.na(pairs$UCSC_RefGene_Name),]
ahrr <- pairs[pairs$UCSC_RefGene_Name=="AHRR",] #174 cpg associated with AHRR
 
cor <- readRDS("tissues/Lung/Correlations_DMP_DEG_new.rds")
ahrr_cor <- cor[cor$gene=="AHRR",] #DMP in AHRR
# 
# signif <- cor[cor$p.adj<0.05,]
# sum(signif$gene=="AHRR")

res <- readRDS("tissues/Lung/DML_results.rds")$Smoking2


#Plot location in chromHMM of DMPs for AHRR based on correlation, I would do volcano plot where y is -log10(p-value), x axis is correlation and color is chromHMM summarized
library(ggplot2)
to_plot_correlations <- merge(res, ahrr_cor, by.x="row.names", by.y = "probe")
to_plot_correlations <- merge(to_plot_correlations, ahrr, by.x="Row.names", by.y = "Name")

#Add chromHMM info
chromhmm_cpgs <- read.csv("data/public/lung_chromhmm.csv")
to_plot_correlations <- merge(to_plot_correlations, chromhmm_cpgs, by.x = "Row.names", by.y = "name_ann")

saveRDS(to_plot_correlations, "../figures/data/ahrr_correlations.rds")

#This code will be in figure_S5.R
# #Improve annotation based on ChromHMM
# table(to_plot_correlations$category) #Only 4 annotated as enhancer and 1 as promoter 
# to_plot_correlations$category[to_plot_correlations$region_chromhmm %in% c("Weak transcription", "Strong transcription")] <- "gene_body"
# to_plot_correlations$category[to_plot_correlations$region_chromhmm %in% c("Active enhancer", "Weak enhancer", "Genic enhancer")] <- "enhancer"
# to_plot_correlations$category[to_plot_correlations$region_chromhmm %in% c("Active TSS", "Flanking TSS")] <- "promoter"
# table(to_plot_correlations$category) #18 annotated as enhancer and 4 as promoter 
# to_test <- to_plot_correlations[to_plot_correlations$p.adj<0.05,]
# table(to_test$category, to_test$cor>0)
# 39/49
# 7/10
# # to_plot_correlations$category[to_plot_correlations$region_chromhmm %in% c("Genic enhancer")] <- "genic enhancer" #3 of the positively correlated with enhancers are genic enhancers
# 
# g <- ggplot(to_plot_correlations) + geom_point(aes(cor, -log10(p.adj), color=category)) + theme_bw() + #p value of correlation
#   ylab("-log10(p-value)") + geom_hline(yintercept=-log10(0.05)) +
#   # scale_color_manual(values=c("#CC6677", "#88CCEE", "#DDCC77", "#117733")) +
#   scale_color_manual(values=c("#CC6677", "#88CCEE", "#DDCC77"), labels=c("Enhancer", "Gene body", "Promoter")) +
#   theme(legend.title = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) + xlab("Correlation (rho)")
# ggsave("../figures/figure_s5/ahrr_correlations.pdf", g, width = 3)




#Everything below is previous analysis


# res <- res[res$adj.P.Val<0.05,]
# table(rownames(res) %in% ahrr$Name) #72 out of the 174 are DMP

#Are all of these 72 previously reported? Plot volcano and color by previously reported
#Check from EWAS the reported ones in AHRR
# library(readr)
# library(ggplot2)
# ewas <- read_tsv("data/AHRR_EWAS.tsv") #http://www.ewascatalog.org/?trait=smoking
# ewas <- ewas[ewas$gene=="AHRR" & !is.na(ewas$beta),]
# # test <- aggregate(beta~cpg, data=ewas, mean)
# # test$same <- sapply(test$cpg, function(probe) identical(sign(test$beta[test$cpg==probe]), sign(res$logFC[rownames(res)==probe])))
# 
# ahrr$signif <- ahrr$Name %in% rownames(res)[res$adj.P.Val<0.05]
# probes <- ahrr$Name[ahrr$signif]
# probes <- as.data.frame(cbind(probes, sapply(probes, function(probe) any(res$logFC[rownames(res) == probe] * ewas$beta[ewas$cpg == probe] > 0))))
# 
# to_plot <- res[rownames(res) %in% ahrr$Name,]
# to_plot$color <- "#CC6677"
# to_plot$color[to_plot$adj.P.Val>0.05] <- "grey"
# to_plot$color[rownames(to_plot) %in% probes$probes[probes$V2=="TRUE"]] <- "#88CCEE"
# g <- ggplot(to_plot) + geom_point(aes(logFC, -log10(adj.P.Val), color=color)) + theme_bw() +
#   ylab("-log10(p-value)") + xlim(-1.5, 1.5) +
#   scale_color_manual(values=c("#CC6677", "#88CCEE", "grey"), labels=c("Previously reported", "New", "Not significant")) +
#   theme(legend.title = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
# 
# ggsave("../figures/figures/others/AHRR_new.pdf", g, device = "pdf", width=3.5)

# #Color by correlated with expression and by location in AHRR:
# to_plot$color <- "grey"
# # to_plot$color[to_plot$adj.P.Val>0.05] <- "grey"
# to_plot$color[rownames(to_plot) %in% signif$probe[signif$cor>0]] <- "#88CCEE"
# to_plot$color[rownames(to_plot) %in% signif$probe[signif$cor<0]] <- "#CC6677"
# 
# #Color whether it is positive or negative
# g2 <- ggplot(to_plot) + geom_point(aes(logFC, -log10(adj.P.Val), color=color)) + theme_bw() +
#   ylab("-log10(p-value)") + xlim(-1.5, 1.5) +
#   scale_color_manual(values=c("#88CCEE", "#CC6677", "grey"), labels=c("Positive", "Negative", "Uncorrelated")) +
#   labs(color="Correlated with expression") + 
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
# 
# ggsave("../figures/figures/others/AHRR_correlation.pdf", g2, device = "pdf", width=4)

#Color by correlations with expression
# to_plot_2 <- merge(to_plot, signif, by.x="row.names", by.y="probe", all.x = T) 
# 
# #Color whether it is positive or negative
# library(RColorBrewer)
# suppressPackageStartupMessages(library(circlize))
# colorRamp2(c(0,1, 50000), brewer.pal(9,"RdBu")[c(1,2,7)])
# g2 <- ggplot(to_plot_2) + geom_point(aes(logFC, -log10(adj.P.Val), color=cor)) + theme_bw() +
#   ylab("-log10(p-value)") + xlim(-1.5, 1.5) +
#   labs(color="Correlation (rho)") + 
#   scale_color_gradient2(midpoint=0, low="#CC6677",
#                         high="#88CCEE", na.value="grey", limits=c(-0.65, 0.65) )  +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
# 
# ggsave("../figures/figures/others/AHRR_correlation_rho.pdf", g2, device = "pdf", width=3.2, height = 1.9)



# chromhmm_cpgs <- read_csv("data/public/lung_chromhmm.csv")
# 
# 
# to_plot$color <- sapply(rownames(to_plot), function(cpg) chromhmm_cpgs$region_chromhmm[chromhmm_cpgs$name_ann==cpg])
# to_plot$color[to_plot$adj.P.Val>0.05] <- "Not significant"
# to_plot$color[to_plot$color=="Active enhancer"] <- "Enhancer"
# to_plot$color[to_plot$color=="Genic enhancer"] <- "Enhancer"
# to_plot$color[to_plot$color=="Weak enhancer"] <- "Enhancer"
# to_plot$color[to_plot$color=="Active TSS"] <- "TSS"
# to_plot$color[to_plot$color=="Flanking TSS"] <- "TSS"
# to_plot$color[to_plot$color=="Strong transcription"] <- "Transcription"
# to_plot$color[to_plot$color=="Weak transcription"] <- "Transcription"
# to_plot$color[to_plot$color=="Weak repressed polycomb"] <- "Other"
# to_plot$color[to_plot$color=="ZNF genes & repeats"] <- "Other"
# to_plot$color <- factor(to_plot$color, levels=c("Enhancer", "TSS", "Transcription", "Other", "Not significant"))
# g3 <- ggplot(to_plot) + geom_point(aes(logFC, -log10(adj.P.Val), color=color)) + theme_bw() +
#   ylab("-log10(p-value)") + xlim(-1.5, 1.5) +
#   labs(color="Chromatin state") + 
#   scale_color_manual(values=c("#CC6677", "#88CCEE", "#117733", "#DDCC77", "grey")) +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
# 
# g3
# ggsave("../figures/figures/others/AHRR_chromHMM_colored.pdf", g3, device = "pdf", width=4, height = 2.5)

# numbers <- to_plot[to_plot$color!="Not significant",]
# pos <- numbers[numbers$logFC>0,]
# neg <- numbers[numbers$logFC<0,]
# table(pos$color)
# table(neg$color)
# 
# #Were these positions in the previous array?
# data_path <- "~/Documents/mn4/scratch/bsc83/bsc83535/GTEx/v9/Oliva/"
# annotation <- read.csv(paste0(data_path, "GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv"))
# 
# to_plot$color <- "#CC6677"
# to_plot$color[rownames(to_plot) %in% annotation$Name[annotation$Methyl450_Loci==TRUE]] <- "#88CCEE"
# 
# g4 <- ggplot(to_plot) + geom_point(aes(logFC, -log10(adj.P.Val), color=color)) + theme_bw() +
#   ylab("-log10(p-value)") + xlim(-1.5, 1.5) +
#   scale_color_manual(values=c("#88CCEE", "#CC6677"), labels=c("Yes", "No")) +
#   labs(color="Was it in the previous array?") + 
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
# 
# ggsave("../figures/figures/others/AHRR_450k.pdf", g4, device = "pdf", width=4)


#Plot number of CpGs / DMPs per gene (and normalized by kb or by number of probes)
# tab_1 <- as.data.frame(table(pairs$UCSC_RefGene_Name))
# tab_1$V4 <- "All"
# # plot(density(tab_1[,2]))
# # ggplot(tab_1, aes(V4, Freq)) + geom_boxplot() + xlab("") + ylab("Number of CpGs per gene") +
# #   geom_text(aes(V4, 174), label="AHRR") +
# #   geom_text(aes(V4, 37), label="CYP1A1") +
# #   geom_text(aes(V4, 7), label="GPR15")
# 
# tab_2 <- as.data.frame(table(cor$gene))
# tab_2$V4 <- "All"
# ggplot(tab_2, aes(V4, Freq)) + geom_violin() + xlab("") + ylab("Number of DMPs per gene") +
#   geom_text(data=tab_2[tab_2$Freq>40,], aes(V4, Freq, label=Var1)) +
#   geom_boxplot(col = "black",
#                outlier.shape = NA,
#                notch = T,
#                width = 0.25) +
#   # geom_text(aes(V4, 72), label="AHRR") +
#   # geom_text(aes(V4, 11), label="CYP1A1") +
#   # geom_text(aes(V4, 0), label="GPR15") +
#   theme_bw() + xlab("") +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
# 
# #Proportion of DMPs/CpGs
# names(tab_2) <- c("Var1", "Freq_DMP")
# tab <- merge(tab_1, tab_2, all=T) #repeted genes
# tab$Freq_DMP[is.na(tab$Freq_DMP)] <- 0
# tab$V5 <- tab$Freq_DMP/tab$Freq
# library(ggrepel)
# ggplot(tab, aes(V4, V5)) + geom_violin() + xlab("") + ylab("Number DMPs / Number CpGs") +
#   geom_text(aes(V4, 0.4137931), label="AHRR") +
#   geom_text_repel(data=tab[tab$V5>0.6,], aes(V4, V5, label=Var1)) +
#   # geom_text(aes(V4, 0.2972973), label="CYP1A1") +
#   # geom_text(aes(V4, 0), label="GPR15") +
#   theme_bw() + xlab("") +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())




# tab_3 <- as.data.frame(table(signif$gene))
# tab_3$V4 <- "All"
# ggplot(tab_3, aes(V4, Freq)) + geom_violin() + xlab("") + ylab("Number of Correlated DMPs per gene") +
#   geom_text(data=tab_3[tab_3$Freq>17,], aes(V4, Freq, label=Var1)) +
#   geom_boxplot(col = "black",
#                outlier.shape = NA,
#                notch = F,
#                width = 0.25) +
#   theme_bw() + xlab("") +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())



# #Proportion of DMPs/CpGs
# names(tab_3) <- c("Var1", "Freq_DMP")
# tab <- merge(tab_1, tab_3, all=T) #repeted genes
# tab$Freq_DMP[is.na(tab$Freq_DMP)] <- 0
# tab$V5 <- tab$Freq_DMP/tab$Freq
# library(ggrepel)
# ggplot(tab, aes(V4, V5)) + geom_violin() + xlab("") + ylab("Number correlated DMPs / Number CpGs") +
#   # geom_text(aes(V4, 0.4137931), label="AHRR") +
#   geom_text_repel(data=tab[tab$V5>0.3,], aes(V4, V5, label=Var1)) +
#   # geom_text(aes(V4, 0.2972973), label="CYP1A1") +
#   # geom_text(aes(V4, 0), label="GPR15") +
#   theme_bw() + xlab("") +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())



# #Other genes: ROBO4, CYP1A1
# robo4 <- pairs[pairs$UCSC_RefGene_Name=="ROBO4",]
# robo4 <- pairs[pairs$UCSC_RefGene_Name=="CYP1A1",]
# robo4 <- pairs[pairs$UCSC_RefGene_Name=="CYP1B1",]
# robo4 <- pairs[pairs$UCSC_RefGene_Name=="F2RL3",]
# to_plot <- res[rownames(res) %in% robo4$Name,]
# 
# to_plot$color <- sapply(rownames(to_plot), function(cpg) chromhmm_cpgs$region_chromhmm[chromhmm_cpgs$name_ann==cpg])
# g3 <- ggplot(to_plot) + geom_point(aes(logFC, -log10(adj.P.Val), color=color)) + theme_bw() +
#   ylab("-log10(p-value)") + 
#   labs(color="ChromHMM region in F2RL3") + 
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
# g3
