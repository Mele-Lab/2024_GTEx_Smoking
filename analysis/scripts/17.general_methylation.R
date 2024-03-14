
#Set path
# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")

library(readr)
chromhmm_cpgs <- read_csv("data/public/lung_chromhmm.csv")

beta <- readRDS("tissues/Lung/methylation_data.rds")
probes <- rownames(beta)
beta <- sapply(beta, as.numeric)
rownames(beta) <- probes

to_plot <- rowMeans(beta)
to_plot <- as.data.frame(to_plot)
to_plot$cpg <- rownames(to_plot)
to_plot <- merge(to_plot, chromhmm_cpgs, by.x="cpg", by.y="name_ann")

library(ggplot2)

library(dplyr)
to_plot %>%  group_by(region_chromhmm) %>% summarise(n=n()) ->Summary.data

to_plot$region_chromhmm <- factor(to_plot$region_chromhmm, levels=c("Active enhancer",
                                                                    "Genic enhancer",
                                                                    "Weak enhancer",
                                                                    "Active TSS", "Flanking TSS",
                                                                    "Strong transcription",
                                                                    "Weak transcription",
                                                                    "Heterochromatin",
                                                                    "Quiescent",
                                                                    "Bivalent enhancer",
                                                                    "Bivalent TSS",
                                                                    "Repressed polycomb",
                                                                    "Weak repressed polycomb",
                                                                    "ZNF genes & repeats"))

# g <- ggplot(to_plot) + geom_violin(aes(region_chromhmm, to_plot), scale="width", fill="lightgray") + xlab("") + ylab("Beta values") +
#   geom_text(data=Summary.data ,aes(x = region_chromhmm, y = 1.1, label=n), fontface =2, size = 3) +
#   theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#                      panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(),)
# g
# ggsave("../figures/figures/others/methylation.pdf", g, device = "pdf", width = 8, height = 4)

g <- ggplot(to_plot) + geom_violin(aes(region_chromhmm, to_plot), scale="width", fill="gray90") + xlab("") + ylab("Beta values") +
  geom_text(data=Summary.data ,aes(x = region_chromhmm, y = 1.1, label=n), fontface =2, size = 3) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title = element_text(size = 14),
                     axis.text = element_text(size = 12.5, colour = "black"))

ggsave("figures/figure_s9/methylation.pdf", g, device = "pdf", width = 6.5, height = 4.5)

# 
# #What about the DMPs in those regions?
# dmps <- readRDS("tissues/Lung/DML_results.rds")
# 
# 
# # hypo <- rownames(dmps$Smoking2)[dmps$Smoking2$adj.P.Val<0.05 & dmps$Smoking2$logFC<0]
# # hyper <- rownames(dmps$Smoking2)[dmps$Smoking2$adj.P.Val<0.05 & dmps$Smoking2$logFC>0]
# # 
# # hypo <- to_plot[to_plot$cpg %in% hypo,]
# # hypo %>%  group_by(region_chromhmm) %>% summarise(n=n()) ->Summary.data
# # 
# # g <- ggplot(hypo) + geom_violin(aes(region_chromhmm, to_plot), trim=T, scale="width") + xlab("Hypo") + ylab("Beta values") +
# #   geom_text(data=Summary.data ,aes(x = region_chromhmm, y = 1.1, label=n), fontface =2, size = 3) +
# #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# # 
# # ggsave("../figures/figures/others/methylation_hypo.pdf", g, device = "pdf", width = 8, height = 4)
# # 
# # hyper <- to_plot[to_plot$cpg %in% hyper,]
# # hyper %>%  group_by(region_chromhmm) %>% summarise(n=n()) ->Summary.data
# # g <- ggplot(hyper) + geom_violin(aes(region_chromhmm, to_plot), trim=T, scale="width") + xlab("Hyper") + ylab("Beta values") +
# #   geom_text(data=Summary.data ,aes(x = region_chromhmm, y = 1.1, label=n), fontface =2, size = 3) +
# #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# # 
# # ggsave("../figures/figures/others/methylation_hyper.pdf", g, device = "pdf", width = 8, height = 4)
# 
# # #What if only never smokers?
# metadata <- readRDS("tissues/Lung/methylation_metadata.rds")
# smokers <- beta[,colnames(beta) %in% metadata$SUBJID[metadata$Smoking==0]]
# to_plot <- rowMeans(smokers)
# to_plot <- as.data.frame(to_plot)
# to_plot$cpg <- rownames(to_plot)
# to_plot <- merge(to_plot, chromhmm_cpgs, by.x="cpg", by.y="name_ann")
# 
# to_plot %>%  group_by(region_chromhmm) %>% summarise(n=n()) ->Summary.data
# 
# to_plot$region_chromhmm <- factor(to_plot$region_chromhmm, levels=c("Active enhancer",
#                                                                     "Genic enhancer",
#                                                                     "Weak enhancer",
#                                                                     "Active TSS", "Flanking TSS",
#                                                                     "Strong transcription",
#                                                                     "Weak transcription",
#                                                                     "Heterochromatin",
#                                                                     "Quiescent",
#                                                                     "Bivalent enhancer",
#                                                                     "Bivalent TSS",
#                                                                     "Repressed polycomb",
#                                                                     "Weak repressed polycomb",
#                                                                     "ZNF genes & repeats"))
# g <- ggplot(to_plot) + geom_violin(aes(region_chromhmm, to_plot), trim=T, scale="width") + xlab("") + ylab("Beta values") +
#   geom_text(data=Summary.data ,aes(x = region_chromhmm, y = 1.1, label=n), fontface =2, size = 3) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# ggsave("../figures/figures/others/methylation_never_smokers.pdf", g, device = "pdf", width = 8, height = 4)
# 
# threshold <- 0
# threshold <- 0.5
# threshold <- 1
# # hypo <- rownames(dmps$Smoking2)[dmps$Smoking2$adj.P.Val<0.05 & dmps$Smoking2$logFC<0]
# hypo <- rownames(dmps$Smoking2)[dmps$Smoking2$adj.P.Val<0.05 & dmps$Smoking2$logFC< -threshold]
# # hyper <- rownames(dmps$Smoking2)[dmps$Smoking2$adj.P.Val<0.05 & dmps$Smoking2$logFC>0]
# hyper <- rownames(dmps$Smoking2)[dmps$Smoking2$adj.P.Val<0.05 & dmps$Smoking2$logFC>threshold]
# 
# hypo <- to_plot[to_plot$cpg %in% hypo,]
# hypo %>%  group_by(region_chromhmm) %>% summarise(n=n()) ->Summary.data
# 
# g <- ggplot(hypo) + geom_violin(aes(region_chromhmm, to_plot), trim=T, scale="width") + xlab("Hypo") + ylab("Beta values") +
#   geom_text(data=Summary.data ,aes(x = region_chromhmm, y = 1.1, label=n), fontface =2, size = 3) +
#   geom_jitter(aes(region_chromhmm, to_plot), alpha=0.6) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# ggsave(paste0("../figures/figures/others/methylation_never_smokers_hypo_", threshold,".pdf"), g, device = "pdf", width = 8, height = 4)
# 
# hyper <- to_plot[to_plot$cpg %in% hyper,]
# hyper %>%  group_by(region_chromhmm) %>% summarise(n=n()) ->Summary.data
# g <- ggplot(hyper) + geom_violin(aes(region_chromhmm, to_plot), trim=T, scale="width") + xlab("Hyper") + ylab("Beta values") +
#   geom_text(data=Summary.data ,aes(x = region_chromhmm, y = 1.1, label=n), fontface =2, size = 3) +
#   geom_jitter(aes(region_chromhmm, to_plot), alpha=0.6) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# ggsave(paste0("../figures/figures/others/methylation_never_smokers_hyper_", threshold, ".pdf"), g, device = "pdf", width = 8, height = 4)
# 
# 
# 
# 
# 
# #get an example of a cpg in never and in smokers:
# hypo <- dmps$Smoking2[dmps$Smoking2$adj.P.Val<0.05 & dmps$Smoking2$logFC<0,]
# hyper <- dmps$Smoking2[dmps$Smoking2$adj.P.Val<0.05 & dmps$Smoking2$logFC>0,]
# hyper_merged <- rownames(hyper)[rownames(hyper) %in% chromhmm_cpgs$name_ann[chromhmm_cpgs$region_chromhmm=="Strong transcription"]]
# # hyper_merged <- rownames(hyper)[rownames(hyper) %in% chromhmm_cpgs$name_ann[chromhmm_cpgs$region_chromhmm=="Active enhancer"]]
# hyper <- hyper[rownames(hyper) %in% hyper_merged,]
# hypo_example <- "cg00605306" #highest hypo
# hypo_example <- "cg08193650" #highest hyper
# hypo_example <- "cg07943658" #hyper example in transcription highest logFC
# hypo_example <- "cg25954028" #hyper example in active enhancer
# hypo_example <- "cg10761372" #hyper example in transcription lowest logFC
# hypo_example <- "cg04306817" #active tss logFC > 0.5 with very low beta value
# hypo_example <- "cg14965300"  #strong transcription hyper with beta 0.9 in never and logFC > 1
# 
# beta_example <- beta[rownames(beta)==hypo_example, ]
# beta_example <- cbind(beta_example, names(beta_example))
# beta_example <- as.data.frame(beta_example)
# beta_example$Smoking <- sapply(beta_example[,2], function(donor) metadata$Smoking[metadata$SUBJID==donor])
# beta_example$beta_example <- as.numeric(beta_example$beta_example)
# beta_example$Smoking <- as.factor(beta_example$Smoking)
# beta_example %>%  group_by(Smoking) %>% summarise(n=n()) ->Summary.data
# g <- ggplot(beta_example) + geom_violin(aes(Smoking, beta_example), trim=T, scale="width") + 
#   geom_jitter(aes(Smoking, beta_example)) + xlab("Smokers") + ylab("Beta values") +
#   geom_text(data=Summary.data ,aes(x = Smoking, y = 0.8, label=n), fontface =2, size = 3)
# g
# 
# 
