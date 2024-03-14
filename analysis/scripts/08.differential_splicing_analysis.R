#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to parse splicing results to share in Zenodo
# @software version: R=4.2.2

Sys.time()
#Set path 
# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")

tissues <- list.dirs("tissues/", full.names = F)[-1]


print("#######################")
print("Reading data")

# Differential splicing analyses: results tables ----
dsa_res <- lapply(tissues, function(tissue) readRDS(paste0("tissues/", tissue, "/fractional_regression_results.rds")))
names(dsa_res) <- tissues

# Hier.part ----
# All AS events in tissue
hier.part.splic <- lapply(tissues, function(tissue) readRDS(paste0("tissues/", tissue, "/Alternatively_spliced_events.hier_part.rds")))
names(hier.part.splic) <- tissues

print("#######################")
print("Adding R2")
#Adding R2
for(tissue in tissues){
  print(tissue)
  #Renaming variables
  names(dsa_res[[tissue]])[names(dsa_res[[tissue]])=="AncestryAMR"] <- "AncestryAMR-AFR"
  names(dsa_res[[tissue]])[names(dsa_res[[tissue]])=="AncestryEUR"] <- "AncestryEUR-AFR"
  names(dsa_res[[tissue]])[names(dsa_res[[tissue]])=="AncestryAMR-AncestryEUR"] <- "AncestryAMR-EUR"
  names(dsa_res[[tissue]])[names(dsa_res[[tissue]])=="Sex2"] <- "Sex" #The reference level was male
  names(dsa_res[[tissue]])[names(dsa_res[[tissue]])=="Smoking1"] <- "SmokingEX-NEVER"
  names(dsa_res[[tissue]])[names(dsa_res[[tissue]])=="Smoking2"] <- "SmokingSMOKER-NEVER"
  names(dsa_res[[tissue]])[names(dsa_res[[tissue]])=="Smoking1-Smoking2"] <- "SmokingEX-SMOKER"
  names(dsa_res[[tissue]])[names(dsa_res[[tissue]])=="Smoking1-2"] <- "SmokingEX-SMOKER"
  
  for(trait in names(dsa_res[[tissue]])){
    dsa_res[[tissue]][[trait]]$R2 <- NA
    if(grepl("Ancestry", trait)){
      trait_name <- "Ancestry"
    } else if(grepl("Smoking", trait)){
      trait_name <- "Smoking"
    } else{
      trait_name <- trait
    }
    for(event in rownames(dsa_res[[tissue]][[trait]])){
      if(dsa_res[[tissue]][[trait]][event, "adj.P.Val"] < 0.05){
        dsa_res[[tissue]][[trait]][event, "R2"] <- 100*hier.part.splic[[tissue]][event, paste0(trait_name, "_abs")]
      }
    }
  }
}
saveRDS(dsa_res, "output/differential_splicing_results_uncomplete.rds")
dsa_res <- readRDS("output/differential_splicing_results_uncomplete.rds")



print("#######################")
print("Adding gene names")
#Adding gene name
gene_annotation <- read.table("data/public/gencode.v26.GRCh38.genes.bed")[,c(6,7)]
names(gene_annotation) <- c("ensembl_id", "gene.name")

for(tissue in tissues){
  print(tissue)
  # for(trait in traits){
  for(trait in names(dsa_res[[tissue]])){
    dsa_res[[tissue]][[trait]] <- dsa_res[[tissue]][[trait]][, c("Estimate", "P.Value", "adj.P.Val", "Ensembl_id", "R2")]
    colnames(dsa_res[[tissue]][[trait]])[1] <- "beta"
    colnames(dsa_res[[tissue]][[trait]])[4] <- "ensembl_id"
    dsa_res[[tissue]][[trait]]$gene_name <- sapply(dsa_res[[tissue]][[trait]]$ensembl_id, function(i) gene_annotation[gene_annotation$ensembl_id == i, "gene.name"])
    dsa_res[[tissue]][[trait]] <- dsa_res[[tissue]][[trait]][,c(1,2,3,5,4,6)]
  }
}


save.image(file = "output/differential_splicing_results_uncomplete_2.RData")
# load(paste0(first_dir, "/Jose/04_Smoking/Tables/Saving_table_2.RData"))

# add which is the most abundant spliced-in and spliced-out isofomrs
# isoform pair per splicing event ----

#Event info with event and transcript
print("#######################")
print("Most expressed isoforms")
# event_annotation <- readRDS(paste0(first_dir, "Raquel/Draft/SUPPA/gencode.v26.PC_lincRNA.biotype_matched_v38.splicing_events_coordinates.rds"))

for(tissue in tissues){
  file <- paste0("tissues/", tissue, "/events_transcripts_DSE.csv")
  if(file.exists(file)){
    print(tissue)
    event_transcripts <- read.csv(file)
    # event_transcripts <- read.csv(paste0(first_dir, "Jose/03_Models/Transcripts/", tissue, "_events_transcripts_all_ASE.csv"))
    #Transcript TPMs, to later keep the most expressed ones:
    # transcript.tpm <- read.delim(paste0(first_dir, "Raquel/Draft/SUPPA/TranscriptExpressionFiles/", tissue, ".transcript_TPM.txt"))
    # colnames(transcript.tpm) <- gsub("\\.", "-", colnames(transcript.tpm))
    # transcript.tpm <- transcript.tpm[, metadata$Sample]
    for(trait in names(dsa_res[[tissue]])){
    # for(trait in traits){
      dsa_res[[tissue]][[trait]]$spliced_in <- as.character(sapply(rownames(dsa_res[[tissue]][[trait]]), function(event) event_transcripts$isoform.in[event_transcripts$event_id==event]))
      dsa_res[[tissue]][[trait]]$spliced_in[dsa_res[[tissue]][[trait]]$spliced_in=="character(0)"] <- "NA"
      dsa_res[[tissue]][[trait]]$spliced_out <- as.character(sapply(rownames(dsa_res[[tissue]][[trait]]), function(event) event_transcripts$isoform.out[event_transcripts$event_id==event]))
      dsa_res[[tissue]][[trait]]$spliced_out[dsa_res[[tissue]][[trait]]$spliced_out=="character(0)"] <- "NA"
    }
  }else{ #Not DSE, only ASE
    for(trait in names(dsa_res[[tissue]])){
    # for(trait in traits){
      dsa_res[[tissue]][[trait]]$spliced_in <- "NA"
      dsa_res[[tissue]][[trait]]$spliced_out <- "NA"
      }
  }
}

save.image(file = "output/differential_splicing_results_uncomplete_3.RData")
# load(paste0(first_dir, "/Jose/04_Smoking/Tables/Saving_table_3.RData"))

print("#######################")
print("Adding biotype")
# biotype_data <- readRDS(paste0(first_dir, "Raquel/Draft/Analysis/splicing/data/event_isoform_data.rds"))
transcript_annotation <- read.delim("data/public/gencode.v26.GRCh38.transcripts.bed", header = F)
colnames(transcript_annotation) <- c("chr","start", "end","strand", "feature","ensembl.id", "transcript.id","gene.name", "gene_biotype", "transcript.name", "transcript_biotype")
# transcript_annotation <- rbind(transcript_annotation, c("NA", 11 ))
for(tissue in tissues){
  print(tissue)
  for(trait in names(dsa_res[[tissue]])){
  # for(trait in traits){
    subset <- dsa_res[[tissue]][[trait]][dsa_res[[tissue]][[trait]]$adj.P.Val<0.05,]
    in_biotype_s <- sapply(subset$spliced_in, function(e)
      transcript_annotation[transcript_annotation$transcript.id == e, "transcript_biotype"])
    in_biotype_s[in_biotype_s=="protein_coding"] <- "PC"
    in_biotype_s[in_biotype_s!="PC" & in_biotype_s!="NA"] <- "NC"
    in_biotype <- c(in_biotype_s, rep("NA", nrow(dsa_res[[tissue]][[trait]])-length(in_biotype_s)))
    out_biotype_s <- sapply(subset$spliced_out, function(e)
      transcript_annotation[transcript_annotation$transcript.id == e, "transcript_biotype"])
    out_biotype_s[out_biotype_s=="protein_coding"] <- "PC"
    out_biotype_s[out_biotype_s!="PC" & out_biotype_s!="NA"] <- "NC"
    out_biotype <- c(out_biotype_s, rep("NA", nrow(dsa_res[[tissue]][[trait]])-length(out_biotype_s)))
    dsa_res[[tissue]][[trait]]$biotype <- paste0(in_biotype, "-", out_biotype)
  }
}
save.image(file = "output/differential_splicing_results_uncomplete_4.RData")
# load(paste0(first_dir, "/Jose/04_Smoking/Tables/Saving_table_4.RData"))


# add if events disrupts a PFAM domain
pfam_data <- read.delim(paste0("2022_GTExTranscriptome_CellGenomics_fromSplicingEventsToProteinDomains/results/4.2.PFAM-Genomic-Coords-GTF-Merged/PFAM-Alignment-Genomic-Coordinates-GTF-Merged.txt"))
pfam_data <- pfam_data[,-c(5,6,7,11,12)]
pfam_data <- pfam_data[!duplicated(pfam_data),]

print("#######################")
print("Pfam for PC-PC")
for(tissue in tissues){
  print(tissue)
  # for(trait in traits){
  for(trait in names(dsa_res[[tissue]])){
    dsa_res[[tissue]][[trait]]$spliced_in_domains <- as.character(sapply(dsa_res[[tissue]][[trait]]$spliced_in, function(e) pfam_data[pfam_data$transcript_id==e,"pfam_domain"]))
    dsa_res[[tissue]][[trait]]$spliced_in_domains[dsa_res[[tissue]][[trait]]$spliced_in_domains=="character(0)"] <- "NA"
    dsa_res[[tissue]][[trait]]$spliced_out_domains <- as.character(sapply(dsa_res[[tissue]][[trait]]$spliced_out, function(e) pfam_data[pfam_data$transcript_id==e,"pfam_domain"]))
    dsa_res[[tissue]][[trait]]$spliced_out_domains[dsa_res[[tissue]][[trait]]$spliced_out_domains=="character(0)"] <- "NA"
    
  }
}

print("Saving data")
saveRDS(dsa_res, file = "output/00.differential_splicing_analysis.smoking.rds")
