#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to prepare DSA metadata
# @software version: R=4.2.2

#Set path 
setwd(system("pwd", intern = T)) #If in linux
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
# setwd("..")

data_dir <- "/gpfs/scratch/bsc83/bsc83535/GTEx/v8/expression_data/" #My path to the transcript TPMs, which can be downloaded from the GTEx portal
# data_dir <- "~/Documents/mn4/scratch/bsc83/bsc83535/GTEx/v8/expression_data/" #My path to the transcript TPMs, which can be downloaded from the GTEx portal
transcript_tpm <- read.delim(paste0(data_dir, "GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz"), skip=2)
colnames(transcript_tpm) <- gsub("\\.","-", colnames(transcript_tpm))
print("finished reading")
print(transcript_tpm[1:5,1:5])

#Keeping only transcripts of genes of interest
print(dim(transcript_tpm))
gene_annotation <- read.csv("data/public/gene_annotation.csv")
transcript_tpm <- transcript_tpm[transcript_tpm$gene_id %in% gene_annotation$ensembl.id,]
rownames(transcript_tpm) <- transcript_tpm$transcript_id
print(dim(transcript_tpm))

#Reading sample_metadata because brain regions have the same short sample name, so we cannot distinguish with our subset metadata with short name format
sample_metadata <- read.delim("data/protected/GTEx_Sample_Attributes.GRU.txt.gz")
sample_metadata <- sample_metadata[sample_metadata$SMAFRZE=="RNASEQ",]

tissue_info <- read.csv("data/public/tissues_sorted.csv")
tissues <- tissue_info$tissue
  
for(tissue in tissues){
  print(tissue)
  tissue_name <- tissue_info$SMTSD[tissue_info$tissue==tissue]
  transcript_tpm_subset <- transcript_tpm[,colnames(transcript_tpm) %in% sample_metadata$SAMPID[sample_metadata$SMTSD==tissue_name]]
  colnames(transcript_tpm_subset) <- sapply(colnames(transcript_tpm_subset), function(i) paste(unlist(strsplit(i,split = "-"))[1:3],collapse="-" ) ) #This is the shorter id 
  print(transcript_tpm_subset[1:2,1:2])
  metadata <- readRDS(paste0("tissues/", tissue, "/metadata.rds"))
  transcript_tpm_subset <- transcript_tpm_subset[,colnames(transcript_tpm_subset) %in% metadata$Sample] #Samples that we analyze (i.e., samples with smoking annotation)
  print(transcript_tpm_subset[1:2,1:2])
  write.table(transcript_tpm_subset, paste0("SUPPA/TranscripExpressionFiles/", tissue, ".transcript_TPM.txt"),
              col.names = T,
              row.names = T,
              quote = F,
              sep = "\t")
}
print("Everything worked")

# 
# # Functions ####
# # source(paste0(first_dir, "Raquel/R_functions/DEA_and_DSA.R_functions.R"))
# 
# # To save:
# # -- PSI & TPM of alternatively spliced events (ASE)
# # -- hier.part of ASE
# # -- PSI residuals of ASE
# # -- Differential splicing analysis (DSA) results
# 
# print(tissue)
# 
# # Number of CPU (cores) to parallelize mclapply
# n_cores <- 32
# # n_cores <- 2
# 
# 
# # Data ####
# # Gene annotation ----
# # PCG and lincRNA genes with mathched biotype annotation in gencode v38
# # PAR genes excluded
# gene_annotation <- read.delim("data/public/gencode.v26.GRCh38.genes.bed")
# colnames(gene_annotation) <- c("chr","start","end","strand","feature","ensembl.id","gene.name", "biotype","source") #Renaming variables
# gene_annotation <- gene_annotation[gene_annotation$biotype %in% c("protein_coding","lincRNA"),] 
# PAR_genes <- sapply(grep(".Y", gene_annotation$ensembl.id, value = T), function(gene)
#   unlist(strsplit(gene, split = "_"))[[1]])
# gene_annotation <- gene_annotation[-unlist(lapply(PAR_genes, function(gene)
#   grep(gene, gene_annotation$ensembl.id))),]
# 
# #length(unique(gene_annotation$ensembl.id)) # 26,196 genes
# sex.biased.genes <- gene_annotation[gene_annotation$chr=="chrY" |
#                                     gene_annotation$gene.name=="XIST", "ensembl.id"]
# 
# # Transcript annotation ----
# transcript_annotation <- read.delim(paste0(first_dir, "Jose/00_Data/gencode.v26.GRCh38.transcripts.bed"), header = F)
# 
# colnames(transcript_annotation) <- c("chr","start", "end","strand", "feature","ensembl.id", "transcript.id","gene.name", "g.biotype", "transcript.name","t.biotype")
# transcript_annotation <- transcript_annotation[transcript_annotation$ensembl.id %in% gene_annotation$ensembl.id,]
# 
# 
# # Exon annotation ----
# exon_annotation <- read.delim(paste0(first_dir, "Jose/00_Data/gencode.v26.GRCh38.exons.bed"), header = F)
# colnames(exon_annotation) <- c("chr","start", "end","strand", "feature","ensembl.id", "transcript.id","gene.name", "transcript.name","exon.id","exon.numer","g.biotype", "t.biotype")
# exon_annotation <- exon_annotation[exon_annotation$ensembl.id %in% gene_annotation$ensembl.id,]
# 
# # Event annotation ----
# events.info <- readRDS(paste0(first_dir, "Raquel/Draft/SUPPA/gencode.v26.PC_lincRNA.biotype_matched_v38.splicing_events_coordinates.rds"))
# 
# # Metadata ----
# file <- list.files(outpath, pattern = "SampleMetadata", full.names = T)
# 
# metadata <- readRDS(file)
# # metadata$Ancestry <- relevel(metadata$Ancestry, ref="EUR") #No, in smoking, the reference is AFR
# for(disease in diseases){
#   metadata[[disease]] <- as.factor(metadata[[disease]])
# }
# 
# #Testing: including cohort and nucacisonbatch too:
# # metadata_og <- read.delim(paste0(first_dir, "Jose/04_Smoking/Data/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt"))
# # metadata_og <- metadata_og[,c("SAMPID", "SMNABTCHT")]
# # names(metadata_og) <- c("Sample", "NucAcIsoBatch")
# # 
# # metadata_donor_og <- read.delim(paste0(first_dir, "Jose/04_Smoking/Data/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt"))
# # metadata_donor_og <- metadata_donor_og[,c("SUBJID", "COHORT")]
# # names(metadata_donor_og) <- c("Donor", "Cohort")
# # 
# # metadata <- merge(metadata, metadata_og, by="Sample")
# # metadata <- merge(metadata, metadata_donor_og, by="Donor")
# 
# extended <- T  #In this metadata we don't have these variables, but just to double check
# # if(extended==TRUE){
# #   #Extended Model (EM)
# #   metadata <- metadata[, !colnames(metadata) %in% c("Cohort", "NucAcIsoBatch")]
# # } else {
# #   #Base Model (BM)
# #   metadata <- metadata[, !colnames(metadata) %in% c("Cohort", "NucAcIsoBatch", "PEER1")]
# # }
# 
# #if sex tissue, exclude sex
# if(tissue %in% c("Prostate", "Testis", "Vagina", "Uterus", "Ovary")){
#   metadata <- metadata[,names(metadata)!="Sex"]
# }
# 
# if(length(unique(metadata$HardyScale))<5){metadata$HardyScale <- droplevels(metadata$HardyScale)}
# if(length(unique(metadata$Ancestry))<3){metadata$Ancestry <- droplevels(metadata$Ancestry)}
# if("Sex" %in% colnames(metadata)){
#   if(length(unique(metadata$Sex))<2){metadata$Sex <- droplevels(metadata$Sex)}
# }
# 
# # Transcript TPM ----
# transcript.tpm <- read.delim(paste0(first_dir, "Raquel/01_SUPPA/TranscripExpressionFiles/", tissue, ".transcript_TPM.txt"))
# # transcript.tpm <- read.delim(paste0(first_dir, "Raquel/Draft/SUPPA/TranscriptExpressionFiles/", tissue, ".transcript_TPM.txt"))
# 
# colnames(transcript.tpm) <- gsub("\\.", "-", colnames(transcript.tpm))
# colnames(transcript.tpm) <- sapply(colnames(transcript.tpm), function(i) paste(unlist(strsplit(i,split = "-"))[1:3],collapse="-" ) )
# 
# metadata$Sample <- sapply(metadata$Sample, function(i) paste(unlist(strsplit(i,split = "-"))[1:3],collapse="-" ) )
# 
# # Subset to tissue samples --
# transcript.tpm <- transcript.tpm[, metadata$Sample]
# 
# # Reading in PSI and TPM values for alternative splicing events annotated in PC and lincRNA genes ----
# # SUPPA report values like 1.0000000000000002  & 0.9999999999999998 that when read into R appear as 1 but are not recognized internally as 1
# # These types of events would not count as 1 for instance when you try to count the number of samples with PSI == 1
# # Included round(psi,2) in previous step to prevent this issue.
# # psi <- readRDS(paste0(first_dir, "Raquel/Draft/Data/Tissues/",tissue,"/",tissue,".SelectedSamples.PSI_values.Splicing_events.PC_lincRNA.rds"))
# psi <- readRDS(paste0(first_dir, "Raquel/00_Data/Tissues/",tissue,"/",tissue,".PSI_values.Splicing_events.PC_lincRNA.rds"))
# 
# # tpm <- readRDS(paste0(first_dir, "Raquel/Draft/Data/Tissues/",tissue,"/",tissue,".SelectedSamples.TPM.Splicing_events.PC_lincRNA.rds"))
# tpm <- readRDS(paste0(first_dir, "Raquel/00_Data/Tissues/",tissue,"/",tissue,".TPM.Splicing_events.PC_lincRNA.rds"))
# 
# # Expressed genes ----
# exprs_genes <- readRDS(paste0(outpath,tissue,".SelectedSamples.Expressed_genes.rds"))
# 
# 
# # ----Code ---- #####
# is.event.exprs <- function(event.id, tpm.threshold = 0.5){
#   if(which(rownames(psi) == event.id) %% 1000 == 0){ system(paste("echo 'Processed: ", which(rownames(psi) == event.id)," out of ", nrow(psi), "events'"))}
#   #event.id <- rownames(psi)[1]
#   # The two most abundant isforms in numeratos and denominator[!numerator] median TPM >= 1
#   # All anotated isoforms
#   #transcripts.id <- transcript_annotation[transcript_annotation$ensembl.id == events.info[event.id, "ensembl.id"], "transcript.id"]
# 
#   # Isoforms that include the event
#   isoforms.in <- unlist(strsplit(events.info[event.id, "isoforms.spliced_in"],split = ","))
#   # Isoforms that excluded the event
#   isoforms.out <- unlist(strsplit(events.info[event.id, "isoforms.spliced_out"],split = ","))
#   if(is.na(isoforms.in) | is.na(isoforms.out)){return(FALSE)}
# 
#   # Isoforms exprs value for each tissue sample
#   isoforms.tpm <- lapply(metadata$Sample, function(sample)
#     sapply(c(isoforms.in, isoforms.out), function(isoform)
#       transcript.tpm[isoform,sample]
#     ))
#   names(isoforms.tpm) <- metadata$Sample
# 
#   # ---- #
#   # Option 1: Most abundant isoform median TPM >= 1 ----
#   # Most abundant isoform expression (median per isoform across samples)
#   # isoforms.median.tpm <- sapply(c(isoforms.in, isoforms.out), function(isoform)
#   #     median(sapply(metadata$Sample, function(sample)
#   #         isoforms.tpm[[sample]][isoform]
#   #     ))
#   # )
#   # isoforms.median.tpm[names(isoforms.median.tpm) %in% isoforms.in]
#   # isoforms.median.tpm[names(isoforms.median.tpm) %in% isoforms.out]
#   #
#   # condition <- max(isoforms.median.tpm[names(isoforms.median.tpm) %in% isoforms.in]) >  tpm.threshold &
#   #     max(isoforms.median.tpm[names(isoforms.median.tpm) %in% isoforms.out]) > tpm.threshold
#   # return(condition)
#   # ---- #
# 
#   # ---- #
#   # Option 2: At least 20% of the samples express the most abundant isofoorm above 0.5/1 TPM ----
#   # Most abundant isoform expression (median per isoform across samples)
#   iso.in.most_abundant <- names(which.max(sapply(isoforms.in, function(isoform)
#     median(sapply(metadata$Sample, function(sample)
#       isoforms.tpm[[sample]][isoform]
#     ))
#   )))
#   iso.out.most_abundant <- names(which.max(sapply(isoforms.out, function(isoform)
#     median(sapply(metadata$Sample, function(sample)
#       isoforms.tpm[[sample]][isoform]
#     ))
#   )))
#   # 20% of the samples express the most abundant above 1 TPM
#   condition1 <- sum(sapply(metadata$Sample, function(sample)
#     isoforms.tpm[[sample]][iso.in.most_abundant] >= tpm.threshold
#   )) >= round(0.2*nrow(metadata))
#   condition2 <- sum(sapply(metadata$Sample, function(sample)
#     isoforms.tpm[[sample]][iso.out.most_abundant] >= tpm.threshold
#   )) >= round(0.2*nrow(metadata))
# 
#   return(condition1 & condition2)
# 
# }
# 
# 
# 
# # 1. Filter AS events ----
# print("# ---- Selecting alternatively spliced events ---- #")
# print(paste0("Number of ASE in PC and lincRNA genes: ", nrow(psi)))
# 
# # * Track number of alternative splicing events ----
# no.ase <- vector()
# no.ase <- nrow(psi)
# 
# # 1.1 Subset events in PCG and lincRNA expressed in tissue ----
# print("# ---- Subset events in PCG and lincRNA expressed in tissue ---- #")
# 
# # Keep events in expressed PC and lincRNA --
# psi$ensembl_id <- sapply(rownames(psi), function(gene) unlist(strsplit(gene,split = ";"))[[1]])
# tpm$ensembl_id <- sapply(rownames(tpm), function(gene) unlist(strsplit(gene,split = ";"))[[1]])
# 
# psi <- psi[psi$ensembl_id %in% exprs_genes,]
# tpm <- tpm[tpm$ensembl_id %in% exprs_genes,]
# psi <- psi[!psi$ensembl_id %in% sex.biased.genes,]
# tpm <- tpm[!tpm$ensembl_id %in% sex.biased.genes,]
# 
# psi <- psi[,-ncol(psi)]
# tpm <- tpm[,-ncol(tpm)]
# 
# 
# #Subset for testing
# # psi <- psi[1:50,]
# 
# #Jose's addition, keep samples in psi for which we have metadata
# colnames(psi) <- sapply(colnames(psi), function(i) paste(unlist(strsplit(i,split = "-"))[1:3],collapse="-" ) )
# psi <- psi[,colnames(psi) %in% metadata$Sample]  #Not identical, why?
# # metadata <- metadata[metadata$Sample %in% colnames(psi),] #I just added this line
# colnames(tpm) <- sapply(colnames(tpm), function(i) paste(unlist(strsplit(i,split = "-"))[1:3],collapse="-" ) )
# tpm <- tpm[,colnames(tpm) %in% metadata$Sample]
# # tpm <- tpm[rownames(tpm) %in% rownames(psi),] #For testing only
# 
# # * Track number of alternative splicing events ----
# no.ase <- c(no.ase, nrow(psi))
# 
# # 1.2 Retain events with no samples having an NA ----
# print("# ---- Retain events with no samples having an NA ---- #")
# 
# # Retain events with no samples having an NA
# number_na <- rowSums(is.na(psi),na.rm=F) # vector with number of samples with NA per event
# noNA_events.psi <- names(number_na)[number_na==0] # events with 0 NAs
# psi <- psi[noNA_events.psi, ]
# 
# # * Track number of alternative splicing events ----
# no.ase <- c(no.ase, nrow(psi))
# 
# # 1.3 Calculate event residuals ####
# print("# ---- Calculating residuals ---- #")
# 
# # 1.3.1 Variables ----  I changed it becaus if I have more than 1 disease, I want them to be individual traits (e.g., Smoking and Pneumonia)
# covariates <- c("HardyScale","IschemicTime", "RIN", "Cohort", "NucAcIsoBatch", "ExonicRate", "PEER1", "PEER2")
# covariates <- covariates[covariates %in% names(metadata)]
# individual_traits <- names(metadata)[!names(metadata) %in% covariates]
# individual_traits <- individual_traits[!individual_traits %in% c("Donor", "Sample")]
# 
# contrast_names <- c()
# for(trait in individual_traits){
#   if(trait=="Age" | trait == "BMI"){
#     contrast_names <- c(contrast_names, trait)
#   } else if(trait=="Ancestry"){
#     contrast_names <- c(contrast_names, "AncestryAMR", "AncestryEUR")
#   } else if(trait=="Sex"){
#     contrast_names <- c(contrast_names, "Sex2")
#   } else{ #Assuming the rest are diseases
#     metadata[[trait]] <- as.factor(metadata[[trait]])
#     contrast_names <- c(contrast_names, paste0(trait, 1))
#   }
# }
# contrast_names <- c(contrast_names, "Smoking2")
# 
# 
# # 1.3.2 Compute residuals ----
# # -------------- #
# print(Sys.time())
# # -------------- #
# # Residuals    20 min
# fr <- mclapply(rownames(psi), function(event_id) get_residuals(event_id, metadata), mc.cores = n_cores )
# names(fr) <- rownames(psi)
# # -------------- #
# print(Sys.time())
# # -------------- #
# 
# # 1.3.3 Create dataframe ----
# psi_residuals <- do.call(rbind.data.frame,
#                          fr)
# colnames(psi_residuals) <- colnames(psi)
# rownames(psi_residuals) <- rownames(psi)
# psi_residuals <- round(psi_residuals, 2)
# 
# # 1.4 Exclude events with low complexity ----
# print("# ---- Exclude event with low complexity ---- #")
# # exclude events with fewer than max(10, 0.1n) unique values, where n is the sample size
# psi.complexity <- apply(psi, 1, function(x) length(unique(x)))
# # exclude events with fewer than max(10, 0.1n) unique values, where n is the sample size
# psi_residuals.complexity <- apply(psi_residuals, 1, function(x) length(unique(x)))
# 
# threshold_complexity <- 15
# # threshold_complexity <- 10
# 
# psi <- psi[intersect(names(psi.complexity[psi.complexity >= threshold_complexity]),
#                               names(psi_residuals.complexity[psi_residuals.complexity >= threshold_complexity])
# ),]
# psi_residuals <- psi_residuals[intersect(names(psi.complexity[psi.complexity >= threshold_complexity]),
#                                          names(psi_residuals.complexity[psi_residuals.complexity >= threshold_complexity])
#                                          ),]
# 
# 
# # * Track number of alternative splicing events ----
# no.ase <- c(no.ase, nrow(psi))
# 
# # 1.5 Exclude events with insufficient variability ----
# print("# ---- Exclude events with with insufficient variability ---- #")
# event_freq <- apply(psi, 1, function(x) sort(table(x),decreasing = T)[1]/sort(table(x),decreasing = T)[2] < 80/20)
# event_freq_residuals <- apply(psi_residuals, 1, function(x) sort(table(x),decreasing = T)[1]/sort(table(x),decreasing = T)[2] < 80/20)
# psi <- psi[event_freq & event_freq_residuals,]
# psi_residuals <- psi_residuals[event_freq & event_freq_residuals,]
# 
# # * Track number of alternative splicing events ----
# no.ase <- c(no.ase, nrow(psi))
# 
# # 1.6 Exclude events not sufficiently expressed ----
# print("# ---- Exclude events not sufficiently expressed ---- #")
# # -------------- #
# print(Sys.time())
# # -------------- #
# events_exprs <- unlist(mclapply(rownames(psi), function(i) is.event.exprs(i, 0.5),  mc.cores = n_cores  ))
# # events_exprs <- unlist(mclapply(rownames(psi), function(i) is.event.exprs(i, 1),  mc.cores = n_cores  ))
# # -------------- #
# print(Sys.time())
# # -------------- #
# save.image(file = paste0(first_dir, "/Jose/04_Smoking/Splicing/", disease, "_", tissue, "_three_levels_filter_data.RData"))
# 
# 
# 
# 
# # load(paste0(first_dir, "/Jose/04_Smoking/Splicing/", disease, "_", tissue,"_three_levels_data.RData"))
# 
# # Subset ASE sufficiently exprs ----
# psi <- psi[events_exprs,]
# psi_residuals <- psi_residuals[rownames(psi),]
# tpm <- tpm[rownames(psi),]
# 
# # * Track number of alternative splicing events ----
# no.ase <- c(no.ase, nrow(psi))
# 
# # 2. hier.part ----
# print("# ---- Running hier.part ---- #")
# 
# # 2.1 Calculate explained variance ----
# # -------------- #
# print(Sys.time())
# # -------------- #
# hier.part.results <- mclapply(rownames(psi_residuals), function(event)
#   hier.part.mod(y=as.numeric(psi_residuals[event,]), x=metadata[,individual_traits],
#                 fam = "quasibinomial", link = "logit", gof = "Rsqu",control = list(maxit = 100)), mc.cores = n_cores)
# names(hier.part.results) <- rownames(psi_residuals)
# # -------------- #
# print(Sys.time())
# # -------------- #
# 
# # 2.2 Exclude events with negative estimates ----
# print(paste0("ASE events with unestimable contributions: ", sum(is.na(hier.part.results))))
# if(length(which(is.na(hier.part.results)))>0){
#   psi <- psi[-which(is.na(hier.part.results)),]
#   psi_residuals <- psi_residuals[-which(is.na(hier.part.results)),]
#   hier.part.results <- hier.part.results[-which(is.na(hier.part.results))]
# }
# 
# # * Track number of alternative splicing events ----
# no.ase <- c(no.ase, nrow(psi))
# 
# # 2.3 Parse results ----
# rsq <- sapply(rownames(psi_residuals), function(event)
#   sum(hier.part.results[[event]]$IJ[,1]))
# names(rsq) <- rownames(psi_residuals)
# rel_perc <- do.call(rbind.data.frame,
#                     lapply(names(hier.part.results), function(event)
#                       as.numeric(unlist(hier.part.results[[event]]$I.perc))))
# rownames(rel_perc) <- rownames(psi_residuals)
# colnames(rel_perc) <- individual_traits
# abs_perc <-  do.call(rbind.data.frame,
#                      lapply(names(hier.part.results), function(event)
#                        hier.part.results[[event]]$IJ[,1])
# )
# rownames(abs_perc) <- rownames(psi_residuals)
# colnames(abs_perc) <- individual_traits
# hier_data <- cbind.data.frame(rsq,rel_perc,abs_perc)
# colnames(hier_data) <- c("R2",paste0(individual_traits,"_rel"),paste0(individual_traits,"_abs"))
# 
# save.image(file = paste0(first_dir, "/Jose/04_Smoking/Splicing/", disease, "_", tissue, "_hier_part_three_levels_filter_data.RData"))
# 
# # load(paste0(first_dir, "/Jose/04_Smoking/Splicing/", disease, "_", tissue,"_hier_part_three_levels_data.RData"))
# 
# 
# # 3. Differential splicing analysis ----
# print("# ---- Running differential splicing analysis ---- #")
# 
# 
# get_tables <- function(event_id, glm_model, mult, mult_ancestry){
#   df <-  as.data.frame(coef(summary(glm_model)))
#   df <- df[contrast_names,]
#   
#   # df$Trait <- c(individual_traits[-length(individual_traits)], "Smoking1", "Smoking2")
#   df$Trait <- contrast_names
#   df$Ensembl_id <- unlist(strsplit(event_id, split = ";"))[[1]]
#   df$Event_id <- event_id
#   df$Event <- unlist(strsplit(event_id, split = ";"))[[2]]
#   df$Type <- unlist(strsplit(df$Event, split = ":"))[[1]]
#   df$`t value`<- NULL
#   df$`Std. Error` = NULL
#   df$glm.p.value = df$`Pr(>|t|)` 
#   # ----- robust sandwich error #
#   cft <- coeftest(glm_model, vcov.=vcovHC(glm_model, type="HC0"))
#   df$Std.Error <- cft[contrast_names, 2]
#   df$Z.Value <- cft[contrast_names, 3]
#   df$P.Value <- cft[contrast_names, 4]
#   df$Iter <- glm_model$iter
#   #Jose's update: 
#   mult_cft <- coeftest(mult) #To add the last comparison needed Smoking1-2. P-values of Smoking 1 and Smoking2 match with the ones we get here. Some variables could be excluded though
#   df <- rbind(df, c(mult$test$coefficients[3], mult$test$pvalues[3], "Smoking1-2", 
#                     unlist(strsplit(event_id, split = ";"))[[1]], event_id, 
#                     unlist(strsplit(event_id, split = ";"))[[2]], unlist(strsplit(df$Event, split = ":"))[[1]],
#                     mult$test$pvalues[3], mult_cft[3,2:4], glm_model$iter))
#   rownames(df)[nrow(df)] <- "Smoking1-2"
#   mult_cft_ancestry <- coeftest(mult_ancestry) #To add the last comparison needed AncestryAMR-EUR. 
#   df <- rbind(df, c(mult_ancestry$test$coefficients[3], mult_ancestry$test$pvalues[3], "AncestryAMR-EUR", 
#                     unlist(strsplit(event_id, split = ";"))[[1]], event_id, 
#                     unlist(strsplit(event_id, split = ";"))[[2]], unlist(strsplit(df$Event, split = ":"))[[1]],
#                     mult_ancestry$test$pvalues[3], mult_cft_ancestry[3,2:4], glm_model$iter))
#   rownames(df)[nrow(df)] <- "AncestryAMR-EUR"
#   
#   # -------------------------------------- #
#   # If warning, event might be modelled properly and function might return an NA (NaN)
#   # Do this to keep track of warnigns. 
#   # Note how events with warnings are associated with extremely large betas
#   cft <- tryCatch({
#     coeftest(glm_model, vcov.=vcovHC(glm_model, type="HC0"))
#   }, warning = function(w) {
#     return(NA)
#     # "Warning message:
#     #   In sqrt(diag(se)) : NaNs produced"
#   }  
#   )
#   # Subset data used in downstream analysis
#   if(!is.na(cft[1])){
#     df$coeftest_Warning <- rep("0",nrow(df)) # 0 are events with no warnings
#   }else{
#     df$coeftest_Warning <- rep("1", nrow(df)) # 1 are events with  warnings
#   }
#   df <- df[,c("Event_id",
#               "Trait",
#               "Ensembl_id",
#               "Type",  
#               "Event",
#               "Estimate",
#               "glm.p.value",
#               "Iter",
#               "Std.Error",
#               "Z.Value", 
#               "P.Value",
#               "coeftest_Warning")]
#   return(df)
# }
# 
# # Function to fit a glm model per event
# model_psi <- function(event_id, mdata){
#   # Model one event at a time
#   psi_values <- pmin(pmax(as.numeric(psi[event_id,]),0),1) # forzar los límites de PSI a los valores teóricos de mínimo y máximo
#   glm_data <- cbind.data.frame(mdata, psi_values)
#   
#   # Model formula: psi ~ covariates + traits
#   mod_formula <- as.formula(paste("psi_values ~ ", paste(c(covariates, individual_traits), collapse = " + "), collapse = " ") )
#   myglm <- glm(mod_formula, 
#                data = glm_data, 
#                family = quasibinomial('logit'),
#                control = list(maxit = 100))
# 
#   #Compare with ex-smokers:
#   mult <- summary(glht(myglm, mcp(Smoking="Tukey"), vcov. = vcovHC(myglm, type="HC0")))
#   mult_ancestry <- summary(glht(myglm, mcp(Ancestry="Tukey"), vcov. = vcovHC(myglm, type="HC0")))
#   # If there is a glm warning cause algorithm did not converge 
#   my_warning <- tryCatch({
#     glm(mod_formula, 
#         data = glm_data, 
#         family = quasibinomial('logit'),
#         control = list( maxit = 100) )
#   }, warning = function(w) {
#     return(NA)
#     # Warning message:
#     #  glm.fit: algorithm did not converge 
#   }  
#   )
#   
#   # Create table
#   event_results <- get_tables(event_id, myglm, mult, mult_ancestry)
#   
#   if(!is.na(my_warning[1])){
#     event_results$glm_Warning <- rep("0", nrow(event_results)) # 0 are events with no warnings
#   }else{
#     event_results$glm_Warning <- rep("1", nrow(event_results)) # 1 are events with  warnings
#   }
#   
#   return(list('Event_id' = event_id,
#               'glm' = myglm,
#               'res' = event_results))  # table with estimates and P.Values per event
# }
# 
# # 3.1 Run PSI models ----
# # -------------- #
# print(Sys.time())
# # -------------- #
# # One model per event
# fr <- mclapply(rownames(psi), function(event_id) model_psi(event_id, metadata), mc.cores = n_cores )
# names(fr) <- rownames(psi)
# # -------------- #
# print(Sys.time())
# # -------------- #
# 
# # 3.2 Parse results table ----
# # individual_traits <- c(individual_traits[-length(individual_traits)], "Smoking1", "Smoking2", "Smoking1-2")
# individual_traits <- c(contrast_names, "Smoking1-2", "AncestryAMR-EUR")
# 
# results <- do.call(rbind.data.frame,
#                    lapply(rownames(psi), function(event_id) fr[[event_id]][['res']]))
# dsa_res <- lapply(individual_traits, function(trait)
#   results[results$Trait==trait,])
# names(dsa_res) <- individual_traits
# 
# for(trait in individual_traits){
#   # Set event ID as as rownames
#   rownames(dsa_res[[trait]]) <- dsa_res[[trait]]$Event_id
#   dsa_res[[trait]] <- dsa_res[[trait]][,-1]
# }
# 
# # 3.3 Exclude events with warnings ----
# # -- Beware of warnings -- #
# # The coeftest function raises warning at particular events when the variance-covariance matrix cannot be computed
# # In extreme instances it returns NA, for example, instances of events completely stratified that we excluded before-hand
# # Yet, some events might remain that cannot be modelled
# # If coef.test raises a warning we consider those events are not properly modelled and assign a P.Value of NA
# # Multiple testing pvalue correction ->
# # if in the p-values vector (not corrected) there is an NA, the p.adjust() does not consider it in the n (number of observations)
# # Remove events with warnings in coeftest function
# # These warnings appear in tissues with low sample size for events with low variance
# # We consider these events cannot be modelled
# print(paste0("ASE that raised glm warning: ", sum(dsa_res$Age$glm_Warning == "1") ))
# print(paste0("ASE that raised coef.test warning: ", sum(dsa_res$Age$coeftest_Warning == "1") ))
# 
# # * Track number of alternative splicing events ----
# no.ase <- c(no.ase, nrow(psi) - sum(dsa_res$Age$glm_Warning == "1"))
# no.ase <- c(no.ase, nrow(psi) - sum(dsa_res$Age$glm_Warning == "1") - sum(dsa_res$Age$coeftest_Warning == "1"))
# no.ase <- c(no.ase, nrow(dsa_res[["Age"]]))
# 
# if(sum(dsa_res$Age$glm_Warning == "1") > 0 ){
#   for(trait in individual_traits){
#     dsa_res[[trait]]$P.Value[which(dsa_res[[trait]]$glm_Warning=="1")] <- NA
#   }
# }
# if(sum(dsa_res$Age$coeftest_Warning == "1") > 0 ){
#   for(trait in individual_traits){
#     dsa_res[[trait]]$P.Value[which(dsa_res[[trait]]$coeftest_Warning=="1")] <- NA
#   }
# }
# for(trait in individual_traits){
#   dsa_res[[trait]] <- dsa_res[[trait]][!is.na(dsa_res[[trait]]$P.Value),]
# }
# 
# # 3.4 FDR correction ---
# for(trait in individual_traits){
#   dsa_res[[trait]]$adj.P.Val <- p.adjust(dsa_res[[trait]]$P.Value, method = "BH")
# }
# print(paste0("ASE modelled: ", nrow(dsa_res$Age)))
# 
# # 3.5 Save results ----
# trait <- "Age"
# psi <- psi[rownames(dsa_res[[trait]]),]
# psi_residuals <- psi_residuals[rownames(dsa_res[[trait]]),]
# tpm <- tpm[rownames(dsa_res[[trait]]),]
# hier_data <- hier_data[rownames(dsa_res[[trait]]),]
# 
# # PSI values
# saveRDS(psi, 
#         paste0(outpath,tissue,".Alternatively_spliced_events.PSI_values.rds"))
# # TPM values
# saveRDS(tpm, 
#         paste0(outpath,tissue,".Alternatively_spliced_events.TPM.rds"))
# # Save residuals
# saveRDS(psi_residuals,paste0(outpath, tissue,".Alternatively_spliced_events.psi_residuals.rds"))
# # Save hier.part
# saveRDS(hier_data,paste0(outpath,tissue, ".Alternatively_spliced_events.hier_part.rds"))
# # Save filtering
# saveRDS(no.ase,
#         paste0(outpath,tissue,".Alternatively_spliced_events.Filtering.rds"))
# saveRDS(dsa_res,
#         paste0(outpath,
#                tissue,".fractional_regression.covariates_and_traits.results.tmp.rds")	
# )
# 
# # 4. DeltaPSI and average TPM ----
# print("# ---- Calculating deltaPSI ---- #")
# 
# # 5.2 Compute average TPM expression and variance for each event
# avrg_TPM <- apply(tpm, 1, function(x) mean(log2(x+1)))
# var_TPM <- apply(tpm, 1, function(x) var(x))
# 
# AFR_samples <- metadata[metadata$Ancestry=="AFR", "Sample"]
# EUR_samples <- metadata[metadata$Ancestry=="EUR", "Sample"]
# AMR_samples <- metadata[metadata$Ancestry=="AMR", "Sample"]
# Female_samples <- metadata[metadata$Sex=="2", "Sample"]
# Male_samples <- metadata[metadata$Sex=="1", "Sample"]
# obese_samples <- metadata[metadata$BMI>=30, "Sample"]
# normal_samples <- metadata[metadata$BMI<25, "Sample"]
# 
# # 5.2 Compute deltaPSI ----
# get.deltaPSI <- function(trait){ #Add ancestries in the correct way, and double check the continous variables
#   print(trait)
#   # -------------- #
#   print(Sys.time())
#   # -------------- #)
#   
#   # if(trait == "Ancestry"){
#   #   delta.psi <- sapply(rownames(dsa_res[[trait]]), function(event)
#   #     mean(as.numeric(psi[event,EUR_samples]),na.rm=T) - mean(as.numeric(psi[event,AFR_samples]),na.rm=T))
#   # }else 
#   if(trait == "Sex" | trait == "Sex2"){
#     if(tissue %in% c("Ovary","Uterus","Vagina","Testis","Prostate","BreastMammaryTissue_Female","BreastMammaryTissue_Male")){
#       delta.psi <- NA
#     }else{
#       delta.psi <- sapply(rownames(dsa_res[[trait]]), function(event)
#         mean(as.numeric(psi[event,Female_samples]),na.rm=T) - mean(as.numeric(psi[event,Male_samples]),na.rm=T) )
#     }
#   }else if(trait == "Age"){
#     younger_samples <- metadata[metadata$Age<45,"Sample"] # [20-45)
#     older_samples <- metadata[metadata$Age>=45,"Sample"] # [45-70]
#     delta.psi <- sapply(rownames(dsa_res[[trait]]), function(event)
#       mean(as.numeric(psi[event,older_samples]),na.rm=T) - mean(as.numeric(psi[event,younger_samples]),na.rm=T))
#   }else if(trait == "BMI"){
#     delta.psi <-  sapply(rownames(dsa_res[[trait]]), function(event)
#       mean(as.numeric(psi[event,obese_samples]),na.rm=T) - mean(as.numeric(psi[event,normal_samples]),na.rm=T))
#   } else if(trait == "Smoking1"){
#     healthy_samples <- metadata[metadata[["Smoking"]]=="0", "Sample"]
#     disease_samples <- metadata[metadata[["Smoking"]]=="1", "Sample"]
#     delta.psi <- sapply(rownames(dsa_res[[trait]]), function(event)
#       mean(as.numeric(psi[event,disease_samples]),na.rm=T) - mean(as.numeric(psi[event,healthy_samples]),na.rm=T))
#   }else if(trait == "Smoking2"){
#     healthy_samples <- metadata[metadata[["Smoking"]]=="0", "Sample"]
#     disease_samples <- metadata[metadata[["Smoking"]]=="2", "Sample"]
#     delta.psi <- sapply(rownames(dsa_res[[trait]]), function(event)
#       mean(as.numeric(psi[event,disease_samples]),na.rm=T) - mean(as.numeric(psi[event,healthy_samples]),na.rm=T))
#   }else if(trait == "Smoking1-2"){
#     healthy_samples <- metadata[metadata[["Smoking"]]=="1", "Sample"]
#     disease_samples <- metadata[metadata[["Smoking"]]=="2", "Sample"]
#     delta.psi <- sapply(rownames(dsa_res[[trait]]), function(event)
#       mean(as.numeric(psi[event,disease_samples]),na.rm=T) - mean(as.numeric(psi[event,healthy_samples]),na.rm=T))
#   }else if(trait == "AncestryAMR"){
#     delta.psi <- sapply(rownames(dsa_res[[trait]]), function(event)
#       mean(as.numeric(psi[event,AMR_samples]),na.rm=T) - mean(as.numeric(psi[event,AFR_samples]),na.rm=T))
#   }else if(trait == "AncestryEUR"){
#     delta.psi <- sapply(rownames(dsa_res[[trait]]), function(event)
#       mean(as.numeric(psi[event,EUR_samples]),na.rm=T) - mean(as.numeric(psi[event,AFR_samples]),na.rm=T))
#   }else if(trait == "AncestryAMR-EUR"){
#     delta.psi <- sapply(rownames(dsa_res[[trait]]), function(event)
#       mean(as.numeric(psi[event,EUR_samples]),na.rm=T) - mean(as.numeric(psi[event,AMR_samples]),na.rm=T))
#   }else { #I am assuming all other traits are correct diseases:
#     # delta.psi <- NA
#     healthy_samples <- metadata[metadata[[trait]]=="0", "Sample"]
#     disease_samples <- metadata[metadata[[trait]]=="1", "Sample"]
#       delta.psi <- sapply(rownames(dsa_res[[trait]]), function(event)
#         mean(as.numeric(psi[event,disease_samples]),na.rm=T) - mean(as.numeric(psi[event,healthy_samples]),na.rm=T))
#     
#   }
#   # -------------- #
#   print(Sys.time())
#   # -------------- #
#   return(delta.psi)
# }
# 
# # Add deltaPSI, AvgExprs, ExprsVar and order data.frame
# for(trait in individual_traits){
#   if(trait == "Sex"){
#     if(tissue %in% c("Ovary","Uterus","Vagina","Testis","Prostate","BreastMammaryTissue_Female","BreastMammaryTissue_Male")){
#       dsa_res[[trait]] <- NA
#     }else{
#       # Add deltaPSI
#       dsa_res[[trait]]$deltaPSI <- get.deltaPSI(trait)
#       # Add average expression TPM
#       dsa_res[[trait]]$AvgExprs <- sapply(rownames(dsa_res[[trait]]), function(event) avrg_TPM[event])
#       dsa_res[[trait]]$ExprsVar <- sapply(rownames(dsa_res[[trait]]), function(event) var_TPM[event])
#       # Order by FDR
#       dsa_res[[trait]] <- dsa_res[[trait]][order(dsa_res[[trait]]$adj.P.Val),]
#     }
#   }else{
#     # Add deltaPSI
#     dsa_res[[trait]]$deltaPSI <- get.deltaPSI(trait)
#     # Add average expression TPM and variance
#     dsa_res[[trait]]$AvgExprs <- sapply(rownames(dsa_res[[trait]]), function(event) avrg_TPM[event])
#     dsa_res[[trait]]$ExprsVar <- sapply(rownames(dsa_res[[trait]]), function(event) var_TPM[event])
#     # Order by FDR
#     dsa_res[[trait]] <- dsa_res[[trait]][order(dsa_res[[trait]]$adj.P.Val),]
#   }
# }
# 
# print(paste0(outpath, tissue,".fractional_regression.covariates_and_traits.results_complexity_15.rds")	)
# # 5.3 Save data ----
# saveRDS(dsa_res,
#         paste0(outpath,
#                tissue,".fractional_regression.covariates_and_traits.results_complexity_15.rds")	
# )
# 
# # ---------------------- #
# end_time <- Sys.time()
# print("")
# print("# ---- Elapsed time ---- #")
# end_time - start_time
# print("# ---------------------- #\n")
# 
# print("")
# 
