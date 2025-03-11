#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to create Table S1 and figure 1
# @software version: R=4.2.2

#Set path
# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")
# getwd()

# Reading tissue information: names, abbreviations and colors
tissue_info <- read.csv("data/public/tissue_abreviation.txt")
tissues <- tissue_info$SMTSD #4 tissues that we do not use
tissues <- tissues[!tissues %in% c("Cells - EBV-transformed lymphocytes", "Cells - Cultured fibroblasts")] #Non-tissues are excluded

#Reading protected metadata
donor_metadata <- read.delim("data/protected/GTEx_Subject_Phenotypes.GRU.txt.gz")
donor_metadata <- donor_metadata[, colnames(donor_metadata) %in% c("SUBJID", "DTHHRDY", "AGE", "SEX", "BMI")] #Variables of interest

sample_metadata <- read.delim("data/protected/GTEx_Sample_Attributes.GRU.txt.gz")
sample_metadata <- sample_metadata[sample_metadata$SMAFRZE=="RNASEQ",]
ancestry_metadata <- read.delim("data/protected/inferred_ancestry_838donors.txt")
  

#Reading smoking annotation:
library("readxl")
smoking <- read_excel("data/protected/exSmoker_annotation.xlsx") #814 Manual curation of the protected smoking annotation
smoking <- smoking[!smoking$MHSMKTP %in% c("Cigar", "Pipe", "Other"),] #792  Excluding non-cigarette smokers
smoking <- smoking[!smoking$SmokerStatus=="unknown",] #718
smoking$Smoking <- 0  #never-smokers
smoking$Smoking[smoking$SmokerStatus=="smoker"] <- 2
smoking$Smoking[smoking$SmokerStatus=="ex-smoker"] <- 1
smoking <- smoking[,c("SUBJID", "SmokerStatus", "Smoking")]
smoking <- smoking[!is.na(smoking$SUBJID),]
names(smoking)[1] <- c("Donor")

smoking <- smoking[,-2]
smoking$Smoking <- as.factor(smoking$Smoking) #We need smoking as a factor for downstream analysis
#Getting age distribution:
age_dist <- merge(smoking, donor_metadata, by.x="Donor", by.y="SUBJID")

#We exclude the Asian ancestries:

age_dist <- merge(age_dist, ancestry_metadata, by.x="Donor", by.y="ID")
age_dist <- age_dist[age_dist$inferred_ancestry!="ASN",]

#Plot distribution of all donors we use:
library(ggplot2)
library(showtext)
library(sysfonts)
font_add_google("Roboto", "roboto")
g <- ggplot(age_dist) + geom_density(aes(AGE), fill = "#88CCEE", alpha = 0.5) + theme_classic() + xlab("Age") + ylab(NULL) +
  theme(axis.text.y = element_text(size = 8, color = "black"),
        text = element_text(family = "roboto", color = "black"),
        axis.text = element_text(family = "roboto", color = "black")) 
g
ggsave("Plots/age_distribution.png", g, device = "png", width = 2.3, height = 1.25)

smoking_data <- data.frame(table(age_dist$Smoking))

# Create the bar plot
g2 <- ggplot(smoking_data, aes(x = Freq, y = Var1, fill = Var1)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Freq), hjust = 1.2) + 
  theme_classic() +
  xlab("Number of donors") +
  ylab(NULL) +
  ggtitle(NULL) +
  theme(legend.position = "none",
        text = element_text(family = "roboto", color = "black"),
        axis.text = element_text(family = "roboto", color = "black"))  +
  scale_fill_manual(values = c("#88CCEE", "#DDCC77", "#CC6677")) +
  scale_y_discrete(labels = c("0" = "Never smokers", "1" = "Ex-smokers", "2" = "Smokers"))
g2

ggsave("Plots/smoking_distribution.png", g2, device = "png", width = 2.5, height = 1.25)


#Plot the age distributions per tissue


#Keep samples with RNAseq as some samples are WGS:

dim(sample_metadata)
sample_metadata$SUBJID <- sapply(sample_metadata$SAMPID, function(i) paste(unlist(strsplit(i,split = "-"))[1:2],collapse="-" ) ) #Getting donor ID based on the sample ID
sample_metadata <- sample_metadata[sample_metadata$SUBJID %in% age_dist$Donor,]
dim(sample_metadata)
sample_metadata <- sample_metadata[,c("SAMPID", "SUBJID", "SMTSD")]

#Adding age, bmi, ancestry and sex to sample metadata:
sample_metadata <- merge(sample_metadata, age_dist, by.x="SUBJID", by.y="Donor")

#Removing tissues with less than 80 samples:
library(dplyr)

filtered_metadata <- sample_metadata %>%
  group_by(SMTSD) %>%
  filter(n() >= 80) %>%
  ungroup()

#Get a table with the mean and sd for all variables

library(dplyr)

# Assuming your data is stored in a data frame called df
summary_table <- filtered_metadata %>%
  group_by(SMTSD) %>%
  summarise(
    total_samples = n(),  # Count total samples per tissue
    mean_age = mean(AGE, na.rm = TRUE),
    sd_age = sd(AGE, na.rm = TRUE),
    mean_BMI = mean(BMI, na.rm = TRUE),
    sd_BMI = sd(BMI, na.rm = TRUE),
    AFR_count = sum(inferred_ancestry == "AFR", na.rm = TRUE),
    AMR_count = sum(inferred_ancestry == "AMR", na.rm = TRUE),
    EUR_count = sum(inferred_ancestry == "EUR", na.rm = TRUE),
    Smoking_0 = sum(Smoking == 0, na.rm = TRUE),
    Smoking_1 = sum(Smoking == 1, na.rm = TRUE),
    Smoking_2 = sum(Smoking == 2, na.rm = TRUE),
    Male_count = sum(SEX == 1, na.rm = TRUE),
    Female_count = sum(SEX == 2, na.rm = TRUE)
  )

summary_table <- summary_table[!summary_table$SMTSD %in% c("Cells - Cultured fibroblasts", "Cells - EBV-transformed lymphocytes"),]

colnames(summary_table) <- c("Tissue", "Sample size", "Age (mean years)", "Age (sd)", "BMI (mean)", "BMI (sd)", "Number of African-descent individuals", 
                             "Number of admixed individuals", "Number of European-descent individuals", "Number of never smokers", "Number of ex-smokers", "Number of smokers",
                             "Number of males", "Number of females")
summary_table <- summary_table[,c(1,2,10,11,12,3,4,5,6,7,8,9,13,14)]

# Print the summary table
print(summary_table)

library("xlsx")
summary_table <- as.data.frame(summary_table)  # Convert tibble to dataframe
write.xlsx(summary_table, "output/Supplementary_table_data_distribution.xlsx", 
           col.names = TRUE, row.names = F, append = FALSE)

#Old, plotting only age distribution

# 
# 
# library(ggplot2)
# library(dplyr)
# library(gridExtra)
# 
# #Adding tissue name:
# tissue_info <- read.csv("data/public/tissues_sorted.csv")
# filtered_metadata <- merge(filtered_metadata, tissue_info, by="SMTSD")
# 
# # Bar plot
# filtered_metadata$Smoking <- factor(filtered_metadata$Smoking, levels = c("2", "1", "0"))
# p1 <- ggplot(filtered_metadata) + 
#   geom_bar(data = filtered_metadata, aes(x = reorder(Name, table(Name)[Name]), fill = Smoking)) +  
#   coord_flip() + 
#   scale_fill_manual(values = c("0" = "#88CCEE", "1" = "#DDCC77", "2" = "#CC6677"), 
#                     labels = c("0" = "Never smoker", "1" = "Ex-smoker", "2" = "Smoker")) + 
#   theme_bw() +
#   xlab(NULL) +
#   ylab("Count") +
#   theme(panel.background = element_blank(),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         text = element_text(color = "black"),
#         axis.text = element_text(color = "black"),
#         legend.title = element_blank(), 
#         legend.position = "right"
#   )
# p1
# # Violin plot
# 
# p2 <- ggplot(filtered_metadata) + 
#   geom_violin(aes(x = reorder(Name, table(Name)[Name]), y = AGE, fill = color), 
#               trim = TRUE, alpha = 0.5) +
#   geom_boxplot(aes(x = reorder(Name, table(Name)[Name]), y = AGE, fill = color), 
#                width = 0.5, outlier.shape = NA, alpha = 0.9) +  coord_flip() + 
#   scale_fill_identity() +
#   xlab(NULL) + 
#   ylab("Age") +
#   theme_bw() +
#   theme(panel.background = element_blank(),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         text = element_text(color = "black"),
#         axis.text = element_text(color = "black"),
#         axis.text.y = element_blank()  
#   )
# p2
# 
# library(ggpubr)  # Make sure to load ggpubr
# 
# # Arrange the plots with a shared legend on the right
# p_combined <- ggarrange(p1, p2, 
#                         ncol = 2, 
#                         widths = c(2, 1.2), 
#                         common.legend = TRUE, 
#                         legend = "right")
# 
# # Save as PDF
# pdf("Plots/combined_plots.pdf", width = 8, height = 6)
# print(p_combined)
# dev.off()

