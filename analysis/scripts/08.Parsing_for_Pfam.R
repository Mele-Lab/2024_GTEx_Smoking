#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to analyze events consequences at the functional level
# @software version: R=4.2.2

#Set path 
# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")


event_annotation <- readRDS("SUPPA/gencode.v26.splicing_events_coordinates.rds")

for(tissue in list.dirs("tissues/", full.names = F)[-1]){
    print(tissue)
  
  #Get my DSE
  dsa_res <- readRDS(paste0("tissues/", tissue, "/fractional_regression_results.rds"))
  DSE <- rownames(dsa_res$Smoking2[dsa_res$Smoking2$adj.P.Val<0.05,])
  # my_ASE <- rownames(dsa_res$Age)
  my_ASE <- DSE
  if(length(my_ASE)==0){
    print("no DSG")
    next
  }
  
  # Metadata ----
  metadata <- readRDS(paste0("tissues/",tissue,"/metadata.rds"))

  #Transcript TPMs, to later keep the most expressed ones:
  transcript.tpm <- read.delim(paste0("SUPPA/TranscripExpressionFiles/", tissue, ".transcript_TPM.txt"))
  colnames(transcript.tpm) <- gsub("\\.", "-", colnames(transcript.tpm))

  # Subset to tissue samples --
  transcript.tpm <- transcript.tpm[, metadata$Sample]
  print(paste("finished reading data:", tissue))
  

  #Itereting per DSG
  table <- c(event_id="", isoform.in="", isoform.out="")
  for(event_id in my_ASE){
    spliced.in.isoforms <- unlist(strsplit(event_annotation[event_id, "isoforms.spliced_in"], split = ","))
    spliced.out.isoforms <- unlist(strsplit(event_annotation[event_id, "isoforms.spliced_out"], split = ","))
    spliced.isoforms <- c(spliced.in.isoforms,
                          spliced.out.isoforms)
    
    tpm.data <- cbind.data.frame("Sample" = rep(metadata$Sample, length(spliced.isoforms)),
                                 "TPM" = unlist(lapply(spliced.isoforms, function(isoform) 
                                   as.numeric(transcript.tpm[isoform,])))
    )
    tpm.data$Isoform <- unlist(lapply(spliced.isoforms, function(isoform) rep(isoform, length(metadata$Sample))))
    tpm.data$Class <- ifelse(tpm.data$Isoform %in% spliced.in.isoforms,
                             "Spliced-in",
                             ifelse(tpm.data$Isoform %in% spliced.out.isoforms,
                                    "Spliced-out",       
                                    "Other"))
    
    # Order isoform by median TPM --
    spliced.in.isoforms <- names(sort(sapply(spliced.in.isoforms, function(isoform) median(tpm.data[tpm.data$Isoform==isoform, "TPM"])), decreasing = T))
    spliced.out.isoforms <- names(sort(sapply(spliced.out.isoforms, function(isoform) median(tpm.data[tpm.data$Isoform==isoform, "TPM"])), decreasing = T))
    
    table <- rbind(table, c(event_id, spliced.in.isoforms[1], spliced.out.isoforms[1]))
  }
  
  table <- table[-1,]
  
  #save table in transcripts folder
  write.csv(table, paste0("tissues/", tissue, "/events_transcripts_DSE.csv"), row.names = F)
}


#Part 2
#Only keep events with at least one protein coding transcript
transcript_annotation <- read.delim("data/public/gencode.v26.GRCh38.transcripts.bed", header = F)
colnames(transcript_annotation) <- c("chr","start", "end","strand", "feature","ensembl.id", "transcript.id","gene.name", "gene_biotype", "transcript.name", "transcript_biotype")
transcript_annotation <- transcript_annotation[transcript_annotation$transcript_biotype=="protein_coding",]
transcript_annotation <- transcript_annotation$transcript.id #Transcript ids that are PC

gene_annotation <- read.delim("data/public/gencode.v26.GRCh38.genes.bed", header=F)[,c(6,7)]
colnames(gene_annotation) <- c("gene", "symbol")


# Metadata ----
final_table <- c(event_id="", gene_name="", isoform.in="", isoform.out="")
all_events <- c(event_id="", gene_name="", isoform.in="", isoform.out="")
for(tissue in list.dirs("tissues/", full.names = F)[-1]){
  file <- paste0("tissues/", tissue, "/events_transcripts_DSE.csv")
  if(!file.exists(file)){next}
  print(tissue)
  tab <- read.csv(file)

  if(ncol(tab)==1){ #If only one DSE
    tmp <- c(tab[1,1], tab[2,1], tab[3,1])
    id <- strsplit(tab[1,1], ";")[[1]][1]
    gene_name <- gene_annotation$symbol[gene_annotation$gene==id]
    all_events <- rbind(all_events, c(tmp[1], gene_name, tmp[-1]))
    
    boolean <- (tmp[2] %in% transcript_annotation && tmp[3] %in% transcript_annotation)
    # boolean <- (tmp[2] %in% transcript_annotation | tmp[3] %in% transcript_annotation)
    if(boolean){
      final_table <- rbind(final_table, c(tmp[1], gene_name, tmp[-1]))
    }
    next
  }
  #Subset tab to keep events with both PC transcript
  ids <- sapply(tab[,1], function(event) strsplit(event, ";")[[1]][1])
  names <- sapply(ids, function(id) gene_annotation$symbol[gene_annotation$gene==id])
  all_events <- rbind(all_events, cbind(tab[,1], names, tab[,-1]))
  boolean <- sapply(1:nrow(tab), function(i) (tab[i,2] %in% transcript_annotation && tab[i,3] %in% transcript_annotation))
  # boolean <- sapply(1:nrow(tab), function(i) (tab[i,2] %in% transcript_annotation | tab[i,3] %in% transcript_annotation))
  filtered_tab <- tab[boolean,]
  if(nrow(filtered_tab)!=0){
    final_table <- rbind(final_table, cbind(filtered_tab[,1], names[boolean], filtered_tab[,-1]))
  }
}

final_table <- final_table[-1,]
all_events <- all_events[-1,]
rownames(final_table) <- NULL
rownames(all_events) <- NULL

# nrow(final_table)/nrow(all_events) # Events associated with changes in PC isoforms
# 
# pc_nc <- all_events[all_events$isoform.in %in% transcript_annotation & !all_events$isoform.out %in% transcript_annotation,]
# nc_pc <- all_events[!all_events$isoform.in %in% transcript_annotation & all_events$isoform.out %in% transcript_annotation,]
# (nrow(pc_nc) + nrow(nc_pc))/nrow(all_events) # Events associated with changes between PC-NC isoforms


write.table(final_table, paste0("SUPPA/events_PC-PC.txt"), row.names = F, col.names = F, quote=F, sep=",")

# transcripts <- unique(c(final_table$isoform.in, final_table$isoform.out))
# write.table(transcripts, paste0("SUPPA/transcripts_PC.txt"), row.names = F, col.names = F, quote = F)


#Now create a file with all the events and their coordinates: (using event_annotation)
# event_coord <- as.data.frame(t(sapply(final_table$`filtered_tab[, 1]`, function(event) event_annotation[rownames(event_annotation)==event, c("spliced.in.start", "spliced.in.end", "spliced.out.start", "spliced.out.end")])))
# event_coord$event <- rownames(event_coord)
# event_coord <- event_coord[,c("event", "spliced.in.start", "spliced.in.end", "spliced.out.start", "spliced.out.end")]
# 
# for(trait in colnames(event_coord)){
#   event_coord[[trait]] <- as.character(event_coord[[trait]])
# }
# write.table(event_coord, paste0("SUPPA/events_PC-PC_coordinates.txt"), row.names = F, col.names = F, quote=F, sep=",")



#Can we work with the previously generated results?
pfam_data <- read.delim("../../Pipeline/Transcript-Protein-Domain-Mapper/results_old/4.2.PFAM-Genomic-Coords-GTF-Merged/PFAM-Alignment-Genomic-Coordinates-GTF-Merged.txt")

table(final_table$isoform.in %in% pfam_data$transcript_id)

