#suppressPackageStartupMessages(library(DSS))
suppressMessages(suppressWarnings(library(bsseq)))
suppressMessages(suppressWarnings(library(HDF5Array)))
suppressMessages(suppressWarnings(library(BiocParallel)))
suppressMessages(suppressWarnings(library(rtracklayer)))
suppressMessages(suppressWarnings(library(here)))

#require(bsseq)
library(optparse)
options(warn=-1)
parser <- OptionParser()
parser <- add_option(parser, opt_str=c("-t", "--tissue"), type="character",
                     dest="tissue",
                     help="Tissue")
options=parse_args(parser)
tissue=options$tissue

#Reading input data
# workdir <- "~/Documents/mn4/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/DNA_methylation/"
# workdir <- "~/Documents/mn5/Projects/GTEx_v8/Jose/04_Smoking/old/DNA_methylation/"
# savedir <- "~/Documents/mn5/Projects/GTEx_v8/Jose/04_Smoking/old/DNA_methylation/DMRs_new/"
workdir <- "/gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/old/DNA_methylation/"
savedir <- "/gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/old/DNA_methylation/DMRs_new/"

# tissue <- "Lung"
# tissue <- "caudate"
# tissue <- "thyroid"
# tissue <- "frontal_Cortex"
print(paste0("Starting code for tissue: ", tissue))
Sys.time()

parsing_file <- function(path){
  print(path)
  data <- read.table(paste0(workdir, R.utils::capitalize.default(tissue), "/", path))
  data$sum <- data$V5 + data$V6
  data <- data[,c(1, 2, 7, 5)]
  colnames(data) <- c("chr", "pos", "N", "X") #N=number of total reads, X=number of methylated reads
  return(data)
}

folder_names <- read.table(paste0(workdir, tissue, "_samples.txt"))[,1]
file_names <- read.table(paste0(workdir, tissue, "_samples.txt"))[,2]
file_names <- gsub("\\..*","",file_names)
file_names <- paste0(file_names, "_val_1_bismark_bt2_pe.bismark.cov.gz") #In thyroid now we have the files with .gz, will they work?

#Subset of files that are either never-smokers or smokers for now, and create a vector with their new names, smoker_1...
annotation <- read.csv(paste0(workdir, "metadata.csv"))
if(tissue=="Caudate"){
  name <- "Brain - Caudate (basal ganglia)"
} else if(tissue=="Frontal_Cortex"){
  name <- "Brain - Frontal Cortex (BA9)"
} else if(tissue=="Amygdala"){
  name <- "Brain - Amygdala"
} else if(tissue=="BA24"){
  name <- "Brain - Anterior cingulate cortex (BA24)"
} else if(tissue=="Hippocampus"){
  name <- "Brain - Hippocampus"
} else if(tissue=="Hypothalamus"){
  name <- "Brain - Hypothalamus"
} else if(tissue=="Nucleus"){
  name <- "Brain - Nucleus accumbens (basal ganglia)"
} else if(tissue=="Putamen"){
  name <- "Brain - Putamen (basal ganglia)"
} else{
  # name <- R.utils::capitalize.default(tissue)
  name <- tissue
}
annotation <- annotation[annotation$Tissue==name,]
annotation <- annotation[annotation$Smoking %in% c(0,2),]

real_names <- c()
s <- 0
n <- 0
for(i in 1:nrow(annotation)){
  if(annotation$Smoking[i]==0){
    n <- n + 1
    real_names <- c(real_names, paste0("never_", n))
  } else{
    s <- s + 1
    real_names <- c(real_names, paste0("smoker_", s))
  }
}


metadata_to_share <- annotation
metadata_to_share$dss_name <- real_names
# write.table(metadata_to_share, paste0(savedir, "metadata_to_share_", tissue, "_v1.txt"))

final_file_names <- file_names[Reduce("|", lapply(annotation$SUBJID, grepl, file_names, fixed = TRUE))] #Subset of files for smokers and never smokers only
folder_names <- folder_names[folder_names %in% annotation$SUBJID]
final_names <- paste0(workdir,tissue,'/',folder_names, "/", final_file_names)

### get CpGs

loci_cg <-  findLoci(
  pattern = "CG",
  subject = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, include = paste0('chr',c(1:22)))

print(paste0("Created subset of regions with only GC for tissue: ", tissue))
Sys.time()

# list <- lapply(final_names, parsing_file)
rownames(metadata_to_share) <- metadata_to_share$dss_name
library(parallel)
bsseq <- read.bismark(
  files = final_names,
  loci = loci_cg,
  colData = metadata_to_share,
  rmZeroCov = FALSE,
  BPPARAM = MulticoreParam(workers = 10, progressbar = TRUE),
  BACKEND = "HDF5Array",
  dir = paste0(savedir,'/Data/',tissue,'_WGBS_CG/'),
  verbose = TRUE, replace=TRUE)
#
chr_interest <- unique(bsseq@rowRanges@seqnames)[c(1:24)]
chr_interest <- droplevels(chr_interest)
loci_interest <- which(bsseq@rowRanges@seqnames %in% c(chr_interest))
bsseq <- bsseq[loci_interest,]

print(paste0("Created bsseq object for tissue: ", tissue))
Sys.time()


# list <- mclapply(final_names, parsing_file, mc.cores=40)
# class(list)
# sapply(list, class)
# 
# print("bismark files finished reading")
# save.image(file = paste0(savedir, tissue,"_first_image.RData"))
# # load(paste0(workdir, "first_image.RData"))
# 
# BSobj = makeBSseqData( list, real_names )
# 
saveRDS(bsseq, paste0(savedir, "BSobj_", tissue,"CGs.rds")) #more than 10 GiB

dir.create(paste0(savedir,'Data/',tissue,'_WGBS_CG', ".small_smooth"))
file.copy(
  from = list.files(paste0(savedir,'Data/',tissue,'_WGBS_CG/'), full.names = TRUE),
  to = paste0(savedir,'Data/',tissue,'_WGBS_CG', ".small_smooth"))
bsseq <- loadHDF5SummarizedExperiment(paste0(savedir,'Data/',tissue,'_WGBS_CG', ".small_smooth"))

print(paste0("About to start smoothing for tissue: ", tissue))
Sys.time()

#Small smooth
# BPPARAM <- MulticoreParam(workers =10)
BS_fit_small <- BSmooth(
  BSseq = bsseq,
  ns = 70, #In the paper they use 70
  h = 1000,
  maxGap = 10^8,
  keep.se = FALSE,
  BPPARAM = SerialParam(progressbar = TRUE),
  verbose = TRUE)

#Large smooth
# BS_fit_large <- BSmooth(
#   BSseq = bsseq,
#   ns = 500,
#   h = 20000,
#   maxGap = 10^8,
#   keep.se = FALSE,
#   BPPARAM = MulticoreParam(workers=2,progressbar = TRUE),
#   verbose = TRUE)
# saveRDS(bsseq, '~/marenostrum/Projects/GTEx_v8/Methylation/Data/Lung_WGBS/Lung_large_smooth_six_samples.rds')
# #In the paper they keep the CpGs when at least 1x coverage in all samples, 

print(paste0("Finished smoothing for tissue: ", tissue))
Sys.time()

cov <- getCoverage(BS_fit_small)
keepLoci.ex <- which(rowSums(cov>0)==ncol(cov))
BS_fit_subset_small <- BS_fit_small[keepLoci.ex,]
# summary(getMeth(BS_fit_subset)) #I need to remove all NAs
meth <- getMeth(BS_fit_subset_small)
keepLoci.ex <- which(rowSums(is.na(meth))==0)
length(keepLoci.ex)
BS_fit_subset_small <- BS_fit_subset_small[keepLoci.ex,]
saveRDS(BS_fit_subset_small, paste0(savedir,'Data/',tissue,'_WGBS_CG/',tissue,'_small_smooth_all_samples_filt.CG.rds'))

print(paste0(savedir,'Data/',tissue,'_WGBS_CG/',tissue,'_small_smooth_all_samples_filt.CG.rds', " created for tissue: ", tissue))

# 
# #Same with large smooth
# cov <- getCoverage(BS_fit_large)
# keepLoci.ex <- which(rowSums(cov>0)==ncol(cov))
# BS_fit_subset_large <- BS_fit_large[keepLoci.ex,]
# summary(getMeth(BS_fit_subset_large)) #I need to remove all NAs
# meth <- getMeth(BS_fit_subset_large)
# keepLoci.ex <- which(rowSums(is.na(meth))==0)
# length(keepLoci.ex)
# BS_fit_subset_large <- BS_fit_subset_large[keepLoci.ex,]
# saveRDS(BS_fit_subset_large, '~/marenostrum/Projects/GTEx_v8/Methylation/Data/Lung_WGBS/Lung_large_smooth_six_samples_filt.rds')
# save.image(file = paste0(savedir, tissue,"_second_image.RData"))


#Back to correct pipeline (not multivariate)

BS.small.tstat <- BSmooth.tstat(BS_fit_subset_small,
                                group1 = real_names[grep("never", real_names)],
                                group2 = real_names[grep("smoker", real_names)],
                                local.correct = TRUE, verbose = 2)

pdf(file = paste0(savedir,tissue,"_tstat_small_all.CG.pdf"), width = 7, height = 5)
plot(BS.small.tstat) #we see some hypermethylation
dev.off()


# BS.large.tstat <- BSmooth.tstat(BS_fit_subset_large,
#                                 group1 = real_names[grep("Male", real_names)],
#                                 group2 = real_names[grep("Female", real_names)],
#                                 local.correct = TRUE, verbose = 2)
# saveRDS(BS.large.tstat, '~/marenostrum/Projects/GTEx_v8/Methylation/Data/Lung_WGBS/Lung_large_smooth_six_samples_filt_tstat.rds')
# 
# pdf(file = paste0(savedir, "Plots/",tissue,"_tstat_large_6samples.pdf"), width = 7, height = 5)
# plot(BS.large.tstat) #we see some hypermethylation
# dev.off()


dmrs0 <- dmrFinder(BS.small.tstat)
# dmrs_small <- subset(dmrs0, n >= 2 & abs(meanDiff) >= 0.1) #we filter out DMRs that do not have at least 2 CpGs in them and at least a mean difference (across the DMR) in methylation between never and smokers of at least 0.1
dmrs_small <- subset(dmrs0, n >= 3)  #we filter out DMRs that do not have at least 2 CpGs in them

# #Plotting
# pData <- pData(BS_fit_subset_small)
# pData$col <- rep(c("darkblue","darkred"),times=c(length(real_names[grep("Male", real_names)]),
#                                                  length(real_names[grep("Female", real_names)]))) #male are in blue
# pData(BS_fit_subset_small) <- pData
# # plotRegion(BS_fit, dmrs[1,], extend = 500, addRegions = dmrs)
# 
# pdf(file = paste0(savedir, "Plots/",tissue,"_many_small.CG.all.pdf"), width = 7, height = 5)
# plotManyRegions(BS_fit_subset_small, dmrs_small[1:50,], extend = 3000,
#                 addRegions = dmrs_small)
# dev.off()

saveRDS(dmrs_small, paste0(savedir, tissue,"_dmrs_small.all.CG.rds"))
saveRDS(dmrs0, paste0(savedir, tissue,"_dmrs_small.all.CG.nofilt.rds"))


### stats ####
# chr_number <- paste0('chr',c(1:22))
# dmrs_number <- dmrs_large[dmrs_large$chr %in% chr_number,]
# dmrs_number <- dmrs_small
# 
# table(dmrs_number$chr)
# 
# sum(dmrs_number$meanDiff>0) ## males 39829
# sum(dmrs_number$meanDiff<0) ## females 28242
# table(dmrs_number$direction)
# 
# dmrs_large_f <- subset(dmrs_number, n >= 3 & abs(meanDiff) >= 0.1)
# dmrs_small_f <- subset(dmrs_number, n >= 3 & abs(meanDiff) >= 0.1)
# sum(dmrs_large_f$meanDiff>0) ## males 294
# sum(dmrs_large_f$meanDiff<0) ## females 50
# 
# #### annotate with chromhmm #### 
# files <- list.files('~/marenostrum/Projects/GTEx_v8/Methylation/Data/EpiMap/hg38/', pattern='.bed.gz',full.names=T)
# 
# names_chrom <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "Breast", "MuscleSkeletal", "KidneyCortex", "Testis", "PBMC")
# names_chrom <- c('Lung')
# chromhmm <- lapply(names_chrom, function(tis)
#   read.delim(files[grep(tis, files)], sep='\t', header=F))
# names(chromhmm) <- c('Lung')
# 
# # ann_bed <- dmrs_number[!is.na(annotation$MAPINFO) & !is.na(annotation$CHR),] %>%
# #   dplyr::select(chrom=CHR, start=MAPINFO, end=MAPINFO, name=IlmnID) %>% distinct()
# # ann_bed$chrom <- paste0('chr',ann_bed$chrom)
# # head(ann_bed)
# # ann_bed$start <- ann_bed$start-1
# 
# dmrs0 <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Lung_dmrs_small.all.CG.nofilt.rds')
# 
# library(valr)
# chrom_df <- chromhmm[['Lung']][,c(1:4)]
# colnames(chrom_df) <- c('chrom','start','end','region')
# dmrs_number$chrom <- dmrs_number$chr
# dmrs0$chrom <- dmrs0$chr
# dmrs_regions <- bed_intersect(dmrs0, chrom_df, suffix = c("_dmr", "_chromhmm"))
# 
# dmrs_polycomb <- dmrs_regions[dmrs_regions$region_chromhmm %in% c("ReprPCWk","ReprPC"),]
# dmrs_polycomb <- dmrs_polycomb[,c("chr_dmr", "start_dmr", "end_dmr", "direction_dmr", "meanDiff_dmr")] %>% distinct()
# table(dmrs_polycomb$direction_dmr[abs(dmrs_polycomb$meanDiff_dmr)>0.15])
# 
# dmrs_regions$region_chromhmm_new <- dmrs_regions$region_chromhmm
# dmrs_regions$region_chromhmm_new[dmrs_regions$region_chromhmm %in% c("TssFlnkD", "TssFlnk", "TssFlnkU","TssA")] <- "TSS"
# dmrs_regions$region_chromhmm_new[dmrs_regions$region_chromhmm %in% c("EnhA2", "EnhA1","EnhWk","EnhG1", "EnhG2")] <- "Enh"
# dmrs_regions$region_chromhmm_new[dmrs_regions$region_chromhmm %in% c("ReprPCWk","ReprPC")] <- "ReprPC"
# dmrs_regions$region_chromhmm_new[dmrs_regions$region_chromhmm %in% c("TxWk","Tx")] <- "Tx"
# 
# ### enrichment shared positions
# my_fisher <- function(type){
#   #             DS      Not DS
#   # type
#   # other_types
#   print(type)
#   res <- dmrs_regions
#   chrom_tissue <- dmrs_regions
#   
#   type_df <- chrom_tissue[chrom_tissue$region_chromhmm_new == type,c("chr_dmr", "start_dmr", "end_dmr", "direction_dmr", "meanDiff_dmr", "n_dmr")] %>% distinct()
#   other_type <- chrom_tissue[chrom_tissue$region_chromhmm_new != type,c("chr_dmr", "start_dmr", "end_dmr", "direction_dmr", "meanDiff_dmr", "n_dmr")] %>% distinct()
#   type_diff <- nrow(type_df[(type_df$meanDiff_dmr)>=(0.01) & type_df$n_dmr>=3,])
#   type_notdiff <- nrow(type_df) - type_diff
#   other_type_diff <- nrow(other_type[(other_type$meanDiff_dmr)>=(0.01) & other_type$n_dmr>=3,])
#   other_type_notdiff <- nrow(other_type) - other_type_diff
#   
#   ### test significance
#   m <- matrix(c(type_diff, type_notdiff, other_type_diff, other_type_notdiff), 2,2, byrow = T)
#   print(m)
#   
#   m[is.na(m)] <- 0
#   #m <- m[c(type,paste0('No ',type)),]
#   rownames(m) <- c(type, "Other")
#   colnames(m) <- c("Hyper","Not Hyper")
#   print(m)
#   f <- fisher.test(m)
#   print(f)
#   return(list("f" = f, "m" = type_diff))
# }
# # Two-tailed Fisher test
# #families <- as.vector(unique(shared_cpgs$region_chromhmm))
# families <- c('Enh','EnhBiv','Het','Quies','ReprPC','TSS','TssBiv','Tx','ZNF/Rpts')
# fisher_results <- lapply(families, function(region) my_fisher(region))
# names(fisher_results) <- families
# 
# saveRDS(fisher_results, '~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_hypo_WGBS.001.allsamples.v2.rds')
# 
# ### plot####
# read_data <- function(variables, data, trait){ #Function to prepare data to plot and compute adjusted p value
#   
#   odds_ratio <- lapply(variables, function(type) data[[type]][['f']]$estimate)
#   adj.P.Val <- p.adjust(sapply(variables, function(type) data[[type]][['f']]$p.value), method = "BH")
#   CI_down <- lapply(variables, function(type) data[[type]][['f']]$conf.int[1])
#   CI_up <- lapply(variables, function(type) data[[type]][['f']]$conf.int[2])
#   sample_size <- lapply(variables, function(type) data[[type]][['m']])
#   
#   
#   names(odds_ratio) <- variables
#   names(adj.P.Val) <- variables
#   names(CI_down) <- variables
#   names(CI_up) <- variables
#   names(sample_size) <- variables
#   
#   odds_ratio_df <- as.data.frame(unlist(odds_ratio))
#   odds_ratio_df$label <- variables
#   odds_ratio_df$type <- deparse(substitute(data)) #Either hypo or hyper
#   colnames(odds_ratio_df) <- c('oddsRatio', 'region','type')
#   
#   adj.P.Val_df <- as.data.frame(unlist(adj.P.Val))
#   adj.P.Val_df$label <- variables
#   adj.P.Val_df$type <- deparse(substitute(data))
#   colnames(adj.P.Val_df) <- c('adjPvalue','region','type')
#   
#   CI_down_df <- as.data.frame(unlist(CI_down))
#   CI_down_df$label <- variables
#   CI_down_df$type <- deparse(substitute(data))
#   colnames(CI_down_df) <- c('CI_down','region','type')
#   
#   CI_up_df <- as.data.frame(unlist(CI_up))
#   CI_up_df$label <- variables
#   CI_up_df$type <- deparse(substitute(data))
#   colnames(CI_up_df) <- c('CI_up','region','type')
#   
#   sample_size_df <- as.data.frame(unlist(sample_size))
#   sample_size_df$label <- variables
#   sample_size_df$type <- deparse(substitute(data))
#   colnames(sample_size_df) <- c('sample_size','region','type')
#   
#   all <- Reduce(function(x, y) merge(x, y, all=TRUE), list(odds_ratio_df, adj.P.Val_df, CI_down_df, CI_up_df, sample_size_df))
#   head(all)
#   all$sig <- 'not Sig'
#   all$sig[all$adjPvalue<0.05] <- 'Sig'
#   all <- all[,c("region","oddsRatio","adjPvalue","CI_down","CI_up","sig","type", "sample_size")]
#   return(all)
# }
# 
# hypo <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_hypo_WGBS.001.allsamples.v2.rds')
# hyper <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_hyper_WGBS.001.allsamples.v2.rds')
# 
# colors_traits <- list('SEX'=c('#3B734E','#89AA94'))
# trait <- 'SEX'
# 
# hypo_d <- read_data(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"), hypo, trait)
# hyper_d <- read_data(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"), hyper, trait)
# hyper_hypo <- rbind(hypo_d, hyper_d)
# hyper_hypo$sig[hyper_hypo$sig =="Sig"] <- "FDR < 0.05"
# hyper_hypo$sig[hyper_hypo$sig =="not Sig"] <- "FDR >= 0.05"
# hyper_hypo$sig <- factor(hyper_hypo$sig, levels=c("FDR >= 0.05", "FDR < 0.05"))
# hyper_hypo$region <- factor(hyper_hypo$region, levels=rev(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts")))
# hyper_hypo$type[hyper_hypo$type =="hypo"] <- "Hypomethylation"
# hyper_hypo$type[hyper_hypo$type =="hyper"] <- "Hypermethylation"
# g <- ggplot(hyper_hypo, aes(x=log2(oddsRatio), y=region, colour=type, alpha=sig)) +
#   geom_errorbar(aes(xmin=log2(CI_down), xmax=log2(CI_up)), width=.3) +
#   geom_vline(xintercept = 0) +
#   #xlim(0,20) + #Only for Lung to show the 0
#   geom_point(size=3) + ylab('') + theme_bw() +
#   scale_colour_manual(values=colors_traits[[trait]]) +
#   xlab("log2(Odds ratio)") +
#   scale_alpha_discrete(range = c(0.4, 1), drop = FALSE) +
#   theme(legend.title = element_blank(),
#         axis.text.x = element_text(colour="black", size=13),
#         axis.text.y = element_text(colour="black", size=14),
#         legend.text = element_text(colour="black", size=13),
#         axis.title.x = element_text(size=16),
#         legend.spacing.y = unit(-0.05, "cm"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(colour = "black", linewidth=1)) +
#   scale_y_discrete(breaks=c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"),
#                    labels=c("Enhancer","Enhancer Bivalent","Heterochromatin","Quiescent","Repressed Polycomb","TSS","TSS Bivalent","Transcription","ZNF & Repeats"))# + xlim(0, 3)
# # pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/genomic_location_shared",'_',trait,".pdf"), w = 6, h = 3.5)
# # print(g)
# # dev.off()
# 
# #Plot sample sizes:
# 
# g2 <- ggplot(hyper_hypo) + geom_col(aes(sample_size, region, fill=type), width = 0.6) +
#   theme_classic() + xlab("Number of DMPs") + ylab("") +
#   scale_fill_manual(values=colors_traits[[trait]]) +
#   theme(legend.position = "none",
#         axis.text.x = element_text(colour="black", size=13),
#         axis.text.y=element_blank(),
#         axis.title.x = element_text(size=16)) +
#   scale_y_discrete(breaks=c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"),
#                    labels=c("Enhancer","Enhancer Bivalent","Heterochromatin","Quiescent","Repressed Polycomb","TSS","TSS Bivalent","Transcription","ZNF & Repeats")) #+
# #scale_x_continuous(breaks=c(0, 20000, 40000)) #Only for lung
# # pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/chromhmm/genomic_location_shared",'_',trait,"_sample_size.pdf"), w = 4, h = 3.5)
# # print(g2)
# # dev.off()
# library(ggpubr)
# p <- ggarrange(g, g2, labels = c("A", "B"),
#                common.legend = TRUE, legend = "right", widths = c(0.8,0.3))
# pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/chromhmm/enrichment","WGBS001.vsalldmr.allsamples.pdf"), w = 8, h = 4)
# print(p)
# dev.off()
# 
# 
# ### average meth levels ######
# BS_fit_subset_small <- loadHDF5SummarizedExperiment(paste0(savedir,'/Data/Lung_WGBS_v2_CG', ".small_smooth"))
# cov <- getCoverage(BS_fit_subset_small)
# keepLoci.ex <- which(rowSums(cov>0)==ncol(cov))
# BS_fit_subset_small <- BS_fit_subset_small[keepLoci.ex,]
# # summary(getMeth(BS_fit_subset)) #I need to remove all NAs
# meth <- getMeth(BS_fit_subset_small)
# keepLoci.ex <- which(rowSums(is.na(meth))==0)
# length(keepLoci.ex)
# BS_fit_subset_small <- BS_fit_subset_small[keepLoci.ex,]
# ## final meth data 
# meth <- getMeth(BS_fit_subset_small)
# head(meth)
# summary(as.matrix(meth))
# median <- as.data.frame(colMedians(meth))
# median$sample <- colnames(meth)
# median$third <- matrixStats::colQuantiles(as.matrix(meth), probs = c(0.75))
# median$first <- matrixStats::colQuantiles(as.matrix(meth), probs = c(0.25))
# median$sex <- gsub('_.*','', median$sample)
# 
# #### plot #####
# 
# # Most basic error bar
# pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/all_samples_","WGBSBetas.pdf"), w = 5, h = 4)
# ggplot(median) +
#   geom_bar( aes(x=sample, y=`colMedians(meth)`, fill=sex), stat="identity", alpha=0.7) +
#   geom_errorbar( aes(x=sample, ymin=first, ymax=third), width=0.4, colour="orange", alpha=0.9, size=1.3) +
#   ylab('Median Beta') + theme_bw()
# dev.off()
# 
# ### get methylation levels polycomb ###
# meth <- getMeth(BS_fit_subset_small)