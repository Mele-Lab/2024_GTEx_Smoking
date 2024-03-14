#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to study the common effects of smoking and other demographic traits
# @software version: R=4.2.2

Sys.time()
#Set path 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")

cor <- readRDS("tissues/Lung/Correlations_DMP_DEG_new.rds")
library(readr)
chromhmm_cpgs <- read_csv("data/public/lung_chromhmm.csv")

cor_m <- merge(cor, chromhmm_cpgs, by.x = "probe", "name_ann")
cor_m_int <- cor_m[cor_m$region_chromhmm %in% c("Repressed polycomb", "Bivalent enhancer", "Bivalent TSS"),]

#keep hypermethylated in bivalent
dmps <- readRDS("tissues/Lung/DML_results.rds")$Smoking2
cor_m_int <- cor_m_int[cor_m_int$probe %in% rownames(dmps[dmps$adj.P.Val<0.05 & dmps$logFC>0,]),] #257 out of 279

#Percentage of unique(correlated)
length(unique(cor_m_int$probe[cor_m_int$p.adj<0.05]))/length(unique(cor_m_int$probe)) #From the hypermethylated PCR2 regions, only 4.2% correlate with expression

table(cor_m_int$p.adj<0.05)/nrow(cor_m_int) #4 % correlate

cor_m <- cor_m[!cor_m$probe %in% cor_m_int$probe,]
table(cor_m$p.adj<0.05)/nrow(cor_m) #The bivalent correlate less than the others, as they correlate 12 %

binom.test(11, 257, p=0.1275, alternative = "less") # p-value = 3.578e-06


#What if we account for logFC, because the bivalent DMPs had lower effect sizes
#Add logFC to cor_m_int
test <- merge(cor_m_int, dmps, by.x="probe", by.y="row.names")
test$variable <- 1
test_others <- merge(cor_m, dmps, by.x="probe", by.y="row.names")
test_others$variable <- 0
test <- rbind(test, test_others)

#get a subset of cor_m with a similar logFC distribution to cor_m_int
library(MatchIt)
matching_model <- matchit(variable ~ logFC, data=test, method = "optimal") #Optimal better than nearest neighbour. The matching is optimal in the sense that that sum of the absolute pairwise distances in the matched sample is as small as possible. Advantages: it is less likely that extreme within-pair distances will be large, unlike with nearest neighbor matching
summary(matching_model)

index_model <- as.numeric(matching_model$match.matrix[,1])
test_to_compare <- test[index_model,]

table(test_to_compare$p.adj<0.05)/nrow(test_to_compare) #The bivalent correlate less than the others, as they correlate 15 %
t <- binom.test(11, 257, p=0.1478599, alternative = "less") # p-value = 5.155e-08
p_val <- signif(t$p.value, digits=3)

# wilcox.test(abs(test_to_compare$cor), abs(cor_m_int$cor))

# #Plot comparing correlations. No, we just care about significant correlations
# to_plot <-  test[c(index_model, matching_model$match.matrix[!is.na(matching_model$match.matrix)]),]
# to_plot$variable <- as.factor(to_plot$variable)
# library(ggplot2)
# get_box_stats <- function(y, upper_limit) {
#   return(data.frame(
#     y = 0.95 * upper_limit,
#     label = paste(
#       "N =", length(y), "\n"
#     )
#   ))
# }
# 
# g <- ggplot(data = to_plot, aes(variable, abs(cor))) +
#   geom_violin(aes(fill = variable),
#               col = "black") +
#   geom_boxplot(col = "black",
#                outlier.shape = NA,
#                notch = T,
#                width = 0.25) +
#   theme_bw() + ggtitle(paste0("p-value = ", p_val)) +
#   xlab("") +
#   ylab("abs(cor)") +
#   scale_fill_manual(values = c("#88CCEE",  "#CC6677")) +
#   stat_summary(fun.data = get_box_stats, fun.args = list(upper_limit = max(to_plot$V1) * 1.15),
#                geom = "text", hjust = 0.5, vjust = 0.9, size = 4) +
#   theme(axis.title = element_text(size = 15),
#         axis.text = element_text(size = 13, colour = "black"),
#         panel.background = element_rect(fill = "white"),
#         axis.line = element_line(colour = "grey"),
#         legend.position = "none",
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         plot.title = element_text(hjust = 0.5, size = 13))

#Increase ratio, but until when, I have to check some variable about similarity of mean


#What about the genes for enriched hypermethylated genes?
#Maybe only get the hypermethylated in those areas and check if they are less correlated, because the enrichments can get me the genes but they are repetitive

#Do the enrichments again and check the genes associated with development?
smoking2 <- readRDS("tissues/Lung/DML_results.rds")$Smoking2
signif <- smoking2[smoking2$adj.P.Val<0.05,]
smoking_hyper <- rownames(signif[signif$logFC>0,])

library(missMethyl)
results <- data.frame("ONTOLOGY"="1", "TERM"="1", "N"="1", "DE"="1", "P.DE"="1", "FDR"="1", "region"="1", "direction"="1", "genes"="1")

for(region in c("Active enhancer", "Active TSS", "Flanking TSS", "Quiescent")){
  print(region)
  subset <- chromhmm_cpgs$name_ann[chromhmm_cpgs$region_chromhmm==region]
  subset_hyper <- subset[subset %in% smoking_hyper]
  go_hyper <- gometh(sig.cpg=subset_hyper, all.cpg = subset, collection="GO", array.type="EPIC", sig.genes = T) #Output with genes takes a lot
  go_hyper <- go_hyper[go_hyper$ONTOLOGY=="BP",]
  go_hyper <- go_hyper[go_hyper$FDR<0.05,]
  go_hyper <- go_hyper[grep("development", go_hyper$TERM),]
  genes <- unlist(strsplit(go_hyper$SigGenesInSet, ","))
}


testing <- cor_m[cor_m$gene %in% genes,]
smoking_hyper_int <- smoking_hyper[smoking_hyper %in% chromhmm_cpgs$name_ann[chromhmm_cpgs$region_chromhmm=="Active TSS"]]
testing_2 <- cor_m[cor_m$probe %in% smoking_hyper_int,]
testing_3 <- testing_2[testing_2$probe %in% testing$probe,]

#Are hypermethylated in Active enhancer less correlated?
cor_m_int <-cor_m[cor_m$probe %in% testing_2$probe,]
table(cor_m_int$p.adj<0.05)/nrow(cor_m_int) #14.9 active enhancer, 12 % active TSS

#Are hypermethylated in Active enhancer associated to developmental gene less correlated?
cor_m_int <- cor_m[cor_m$probe %in% testing_3$probe,]
table(cor_m_int$p.adj<0.05)/nrow(cor_m_int) #15 active enhancer, 11 % active TSS
