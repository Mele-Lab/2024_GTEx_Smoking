#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez and Raquel Garcia-Perez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to get the events coordinates
# @software version: R=4.2.2

start_time <- Sys.time()

# This script retrieves the genomic coordinates of 
# alternative splicing events in PC and lincRNA genes
# Check: https://github.com/comprna/SUPPA

# Libraries ----
library(GenomicRanges)
library(rtracklayer)

#Set path 
setwd(system("pwd", intern = T)) #If in linux
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
# setwd("..")

# Gene annotation ----
gene_annotation <- read.csv("data/public/gene_annotation.csv") #PC and lincRNA genes (no PAR)
autosomal.genes <- gene_annotation[!gene_annotation$chr %in% c("chrX", "chrY", "chrM"), "ensembl.id"]
length(unique(gene_annotation$ensembl.id)) # 26,676 genes

# Transcript annotation ----
transcript_annotation <- read.delim("data/public/gencode.v26.GRCh38.transcripts.bed", header = F)
colnames(transcript_annotation) <- c("chr","start", "end","strand", "feature","ensembl.id", "transcript.id","gene.name", "biotype")
transcript_annotation <- transcript_annotation[transcript_annotation$ensembl.id %in% gene_annotation$ensembl.id,]

# Exon annotation ----
exon_annotation <- read.delim("data/public/gencode.v26.GRCh38.exons.bed", header = F)
colnames(exon_annotation) <- c("chr","start", "end","strand", "feature","ensembl.id", "transcript.id","gene.name", "transcript.name","exon.id","biotype")
exon_annotation <- exon_annotation[exon_annotation$ensembl.id %in% gene_annotation$ensembl.id,]

# Create granges object  ----
exon.gr <- makeGRangesFromDataFrame(exon_annotation, keep.extra.columns = T)

# Splicing events ----
events <- rownames(readRDS("tissues/BreastMammaryTissue/psi_splicing_events.rds")) #Any tissue works
length(events) # 161,085 events

# Create splicing event dataframe ----
event.df <- cbind.data.frame("ensembl.id" = sapply(events, function(e) unlist(strsplit(e,split = ";"))[[1]]),
                             "type" = sapply(events, function(e) unlist(strsplit(unlist(strsplit(e,split = ";"))[[2]], split = ":"))[[1]] ),
                             "chr" = sapply(events, function(e) unlist(strsplit(unlist(strsplit(e,split = ";"))[[2]], split = ":"))[[2]] ),
                             "strand" = sapply(events, function(e) unlist(strsplit(unlist(strsplit(e,split = ";"))[[2]], split = ":"))[[length(unlist(strsplit(unlist(strsplit(e,split = ";"))[[2]], split = ":")))]] ),
                             "coords.full" = sapply(events, function(e) paste(unlist(strsplit(unlist(strsplit(e,split = ";"))[[2]], split = ":"))[3:(length(unlist(strsplit(unlist(strsplit(e,split = ";"))[[2]], split = ":")))-1)],
                                                    collapse=":"))
                             )
event.df$ensembl.gene.id <- sapply(event.df$ensembl.id, function(g) unlist(strsplit(g, split = "\\."))[[1]])
event.df$gene.name <- sapply(event.df$ensembl.id, function(g) gene_annotation[gene_annotation$ensembl.id==g, "gene.name"])
length(unique(event.df$ensembl.id)) #16,183

print("Ready to run")
Sys.time()
# For each event, get the coordinates of the exons/introns involved --
get.event.coordinates <- function(event.id){
    # Skipping exon and Retained Intro only have one coordinate --
    # And it does not matter if it is in the + or - strand --
    type <- event.df[event.id,"type"]
    coords.full <- event.df[event.id, "coords.full"]
    strand <- event.df[event.id, "strand"]
    if(type == "SE"){
        # s2:e2 (e1-s2:e2-s3)
        s2 <- unlist(strsplit(unlist(strsplit(coords.full, split = ":"))[[1]], split = "-"))[[2]] # start
        e2 <- unlist(strsplit(unlist(strsplit(coords.full, split = ":"))[[2]], split = "-"))[[1]] # end
        event.coords.1 <- paste0(s2,"-", e2)    # coordinated of the exon that is spliced-in or spliced-out
        event.coords.2 <- NA
    }else if(type == "RI"){
        # e1:s2 (s1:e1-s2:e2)
        e1 <- unlist(strsplit(unlist(strsplit(coords.full, split = "-"))[[1]], split = ":"))[[2]] # start
        s2 <- unlist(strsplit(unlist(strsplit(coords.full, split = "-"))[[2]], split = ":"))[[1]] # end
        event.coords.1 <- paste0(e1,"-", s2)    # coordinated of the exon that is spliced-in or spliced-out
        event.coords.2 <- NA
    }else  if(type == "MX"){
        # It does not matter if it is in the + or - strand --
        # s2:e2 (e1-s2:e2-s4:e1-s3:e3-s4)
        s2 <- unlist(strsplit(unlist(strsplit(coords.full, split = ":"))[[1]], split = "-"))[[2]] 
        e2 <- unlist(strsplit(unlist(strsplit(coords.full, split = ":"))[[2]], split = "-"))[[1]] 
        event.coords.1 <- paste0(s2,"-", e2)  
        # s3:e3 (e1-s2:e2-s4:e1-s3:e3-s4)
        s3 <- unlist(strsplit(unlist(strsplit(coords.full, split = ":"))[[3]], split = "-"))[[2]] 
        e3 <- unlist(strsplit(unlist(strsplit(coords.full, split = ":"))[[4]], split = "-"))[[1]] 
        event.coords.2 <- paste0(s3,"-", e3)  
    }else if(type == "AF" & strand == "+" | type == "AL" & strand == "-"){
        # It matters if it is in the + or - strand --
        # s1:e1 (s1:e1-s3:s2:e2-s3)
        s1 <- unlist(strsplit(coords.full, split = ":"))[[1]]
        e1 <- unlist(strsplit(unlist(strsplit(coords.full, split = ":"))[[2]], split = "-"))[[1]]
        event.coords.1 <- paste0(s1,"-", e1)  
        # s2:e2 (s1:e1-s3:s2:e2-s3)
        s2 <- unlist(strsplit(coords.full, split = ":"))[[3]]
        e2 <- unlist(strsplit(unlist(strsplit(coords.full, split = ":"))[[4]], split = "-"))[[1]]
        event.coords.2 <- paste0(s2,"-", e2)  
    }else if(type == "AF" & strand == "-" | type == "AL" & strand == "+"){
        # It matters if it is in the + or - strand -- # Anotation esta mal en la pagina web! Correcto en el Draft de SUPPA (no SUPPA2)
        # s3:e3 (e1-s3:e3:e1-s2:e2)
        s3 <- unlist(strsplit(unlist(strsplit(coords.full, split = ":"))[[1]], split = "-"))[[2]]
        e3 <- unlist(strsplit(coords.full, split = ":"))[[2]]
        event.coords.1 <- paste0(s3,"-", e3)  
        # s2:e2 (e1-s3:e3:e1-s2:e2)
        s2 <- unlist(strsplit(unlist(strsplit(coords.full, split = ":"))[[3]], split = "-"))[[2]]
        e2 <- unlist(strsplit(coords.full, split = ":"))[[4]]
        event.coords.2 <- paste0(s2,"-", e2)  
        
        # # s3:e3 (e1-s2:e2:e1-s3:e3)
        # s3 <- unlist(strsplit(unlist(strsplit(coords.full, split = ":"))[[3]], split = "-"))[[2]]
        # e3 <- unlist(strsplit(coords.full, split = ":"))[[4]]
        # event.coords.1 <- paste0(s3,"-", e3)  
        # # s2:e2 (e1-s2:e2:e1-s3:e3)
        # s2 <- unlist(strsplit(unlist(strsplit(coords.full, split = ":"))[[1]], split = "-"))[[2]]
        # e2 <- unlist(strsplit(coords.full, split = ":"))[[2]]
        # event.coords.2 <- paste0(s2,"-", e2)  
    }else if(type == "A5" & strand == "+" | type == "A3" & strand == "-"){
        # It matters if it is in the + or - strand --
        # e1:e2 (e2-s3:e1-s3)
        e1 <- unlist(strsplit(unlist(strsplit(coords.full, split = ":"))[[2]], split = "-"))[[1]]
        e2 <- unlist(strsplit(unlist(strsplit(coords.full, split = ":"))[[1]], split = "-"))[[1]]
        event.coords.1 <- paste0(e1,"-", e2)    # coordinated of the exon that is spliced-in or spliced-out
        event.coords.2 <- NA
        
        # # Find overlapping exon
        # gr <- GRanges(
        #     seqnames = event.df[event.id, "chr"],
        #     ranges = IRanges(start = as.numeric(e1)-1, end = as.numeric(e2), names = event.id),
        #     strand = strand)
        # overlaps.index <- as.data.frame(findOverlaps(gr, exon.gr))
        # overlapping.exons <- as.data.frame(exon.gr[overlaps.index$subjectHits,])
        # candidate.exons <- overlapping.exons[overlapping.exons$end==as.numeric(e2),] # the end coordinate must coindice with the end coordinate of the 'black event'
        # if(nrow(candidate.exons)>1){ # if more than one, keep the longest exon
        #     candidate.exons$width <- candidate.exons$end - candidate.exons$start
        #     candidate.exons <- candidate.exons[order(candidate.exons$width, decreasing = T),]
        #     candidate.exons <- candidate.exons[1,]
        # }
        # event.coords.1 <- paste0(candidate.exons$start+1,"-", candidate.exons$end)
    }else if(type == "A5" & strand == "-" | type == "A3" & strand == "+"){
        # It matters if it is in the + or - strand --
        # s2:s3 (e1-s2:e1-s3)
        s2 <- unlist(strsplit(unlist(strsplit(coords.full, split = ":"))[[1]], split = "-"))[[2]]
        s3 <- unlist(strsplit(unlist(strsplit(coords.full, split = ":"))[[2]], split = "-"))[[2]]
        event.coords.1 <- paste0(s2,"-", s3)    # coordinated of the exon that is spliced-in or spliced-out
        event.coords.2 <- NA 
        
        # # Find overlapping exon
        # gr <- GRanges(
        #     seqnames = event.df[event.id, "chr"],
        #     ranges = IRanges(start = as.numeric(s2)-1, end = as.numeric(s3), names = event.id),
        #     strand = strand)
        # overlaps.index <- as.data.frame(findOverlaps(gr, exon.gr))
        # overlapping.exons <- as.data.frame(exon.gr[overlaps.index$subjectHits,])
        # candidate.exons <- overlapping.exons[overlapping.exons$start==as.numeric(s2)-1,] # the start coordinate must coindice with the end coordinate of the 'black event'
        # if(nrow(candidate.exons)>1){ # if more than one, keep the longest exon
        #     candidate.exons$width <- candidate.exons$end - candidate.exons$start
        #     candidate.exons <- candidate.exons[order(candidate.exons$width, decreasing = T),]
        #     candidate.exons <- candidate.exons[1,]
        # }
        # event.coords.1 <- paste0(candidate.exons$start+1,"-", candidate.exons$end)
    }else{
        stop("Incorrect splicing type")
    }
    return(paste0(event.coords.1, ":", event.coords.2 ))
}

# Get the genomic coordinates of the exon(s)/intron involved in the splicing events --
event.coordinates <- lapply(rownames(event.df), function(event.id) get.event.coordinates(event.id))
names(event.coordinates) <- rownames(event.df)

print("Main function finished")
Sys.time()

# Add information to event.df --
event.df$spliced_in.coords <- sapply(rownames(event.df), function(event.id) unlist(strsplit(event.coordinates[[event.id]], split = ":"))[[1]])
event.df$spliced_out.coords <- sapply(rownames(event.df), function(event.id) unlist(strsplit(event.coordinates[[event.id]], split = ":"))[[2]])

# Check --
# 1-based coordinates of event.df --
exon_coords <- paste0(exon_annotation$chr, ":", (exon_annotation$start+1),"-", exon_annotation$end)
# SE
sum(paste0(event.df[event.df$type=="SE",]$chr, ":", event.df[event.df$type=="SE",]$splice_in.coords) %in% exon_coords)/nrow(event.df[event.df$type=="SE",]) # all exon coordinates match
# MX
sum(paste0(event.df[event.df$type=="MX",]$chr, ":", event.df[event.df$type=="MX",]$splice_in.coords) %in% exon_coords)/nrow(event.df[event.df$type=="MX",]) # all exon coordinates match
sum(paste0(event.df[event.df$type=="MX",]$chr, ":", event.df[event.df$type=="MX",]$splice_out.coords) %in% exon_coords)/nrow(event.df[event.df$type=="MX",]) # all exon coordinates match
#  AF
sum(paste0(event.df[event.df$type=="AF",]$chr, ":", event.df[event.df$type=="AF",]$splice_in.coords) %in% exon_coords)/nrow(event.df[event.df$type=="AF",]) # all exon coordinates match
sum(paste0(event.df[event.df$type=="AF",]$chr, ":", event.df[event.df$type=="AF",]$splice_out.coords) %in% exon_coords)/nrow(event.df[event.df$type=="AF",]) # all exon coordinates match
#  AL
sum(paste0(event.df[event.df$type=="AL",]$chr, ":", event.df[event.df$type=="AL",]$splice_in.coords) %in% exon_coords)/nrow(event.df[event.df$type=="AL",]) # all exon coordinates match
sum(paste0(event.df[event.df$type=="AL",]$chr, ":", event.df[event.df$type=="AL",]$splice_out.coords) %in% exon_coords)/nrow(event.df[event.df$type=="AL",]) # all exon coordinates match
#  A5
sum(paste0(event.df[event.df$type=="A5",]$chr, ":", event.df[event.df$type=="A5",]$splice_in.coords) %in% exon_coords)/nrow(event.df[event.df$type=="A5",]) # 
sum(paste0(event.df[event.df$type=="A5",]$chr, ":", event.df[event.df$type=="A5",]$splice_out.coords) %in% exon_coords)/nrow(event.df[event.df$type=="A5",]) # all exon coordinates match
#  A3
sum(paste0(event.df[event.df$type=="A3",]$chr, ":", event.df[event.df$type=="A3",]$splice_in.coords) %in% exon_coords)/nrow(event.df[event.df$type=="A3",]) # 
sum(paste0(event.df[event.df$type=="A3",]$chr, ":", event.df[event.df$type=="A3",]$splice_out.coords) %in% exon_coords)/nrow(event.df[event.df$type=="A3",]) # all exon coordinates match

# For each event, add the transcript.id of the isoforms that contribute to the splicing event ----
ioe <- lapply(unique(event.df$type), function(e) 
    read.delim(paste0("SUPPA/gencode.v26.annotation_events.ioe_", e, "_strict.ioe")))
names(ioe) <-  unique(event.df$type)
isoforms <- do.call(rbind.data.frame,
                    ioe)
event.df$isoforms.spliced_in <- sapply(rownames(event.df), function(e) isoforms[isoforms$event_id==e,"alternative_transcripts"])
event.df$isoforms.spliced_out <- sapply(rownames(event.df), function(e) 
    paste(unlist(strsplit(isoforms[isoforms$event_id==e,"total_transcripts"], split = ","))[
        ! unlist(strsplit(isoforms[isoforms$event_id==e,"total_transcripts"], split = ",")) %in%
            unlist(strsplit(isoforms[isoforms$event_id==e,"alternative_transcripts"], split = ","))
    ], collapse = ",")
)

print("Transcript id added")   
Sys.time()

# Reorganize data ----
# exclude coords.full
event.df <- event.df[, c("ensembl.id",
                         "type",
                         "chr",
                         "strand",
                         "spliced_in.coords",
                         "spliced_out.coords",
                         "isoforms.spliced_in",
                         "isoforms.spliced_out",
                         "ensembl.gene.id",
                         "gene.name")]

# Create 1-based file
event.df.1base <- event.df
event.df.1base$spliced.in.start <- sapply(event.df.1base$spliced_in.coords, function(i) as.numeric(unlist(strsplit(i,split = "-"))[[1]]))
event.df.1base$spliced.in.end <- sapply(event.df.1base$spliced_in.coords, function(i) as.numeric(unlist(strsplit(i,split = "-"))[[2]]))
event.df.1base$spliced.out.start <- sapply(event.df.1base$spliced_out.coords, function(i) ifelse(i=="NA", 
                                                                                                 NA,
                                                                                                 as.numeric(unlist(strsplit(i,split = "-"))[[1]])))
event.df.1base$spliced.out.end <- sapply(event.df.1base$spliced_out.coords, function(i) ifelse(i=="NA", 
                                                                                                 NA,
                                                                                                 as.numeric(unlist(strsplit(i,split = "-"))[[2]])))
event.df.1base <- event.df.1base[, c("ensembl.id",
                                     "type",
                                     "chr",
                                     "strand",
                                     "spliced.in.start",
                                     "spliced.in.end",
                                     "isoforms.spliced_in",
                                     "spliced.out.start",
                                     "spliced.out.end",
                                     "isoforms.spliced_out",
                                     "ensembl.gene.id",
                                     "gene.name")]

# Create and add event labels ----
labels.tmp <- make.names(paste0(gsub("\\.", "-", event.df.1base$ensembl.gene.id), ":", event.df.1base$type, ":1"), unique = T)
labels <- sapply(labels.tmp, function(i) ifelse(length(unlist(strsplit(i, split = "\\.")))==3,
                                                paste0(unlist(strsplit(i, split = "\\."))[c(1,2,3)], collapse = ":"),
                                                    paste0(c(unlist(strsplit(i, split = "\\."))[c(1,2)], 
                                                           as.numeric(unlist(strsplit(i, split = "\\."))[[4]])+1),
                                                        collapse = ":")))
names(labels) <- NULL
event.labels <- sapply(1:length(labels), function(i) paste0(event.df.1base$gene.name[i],
                                                  ":", 
                                            unlist(strsplit(labels[i],split = ":"))[[2]],
                                            ":",
                                            unlist(strsplit(labels[i],split = ":"))[[3]],
                                            ":",
                                            unlist(strsplit(labels[i],split = ":"))[[1]]))
#event.df.1base$event.label <- 
event.df.1base$event.label <- event.labels
write.table(event.df.1base,
            "SUPPA/gencode.v26.splicing_events_coordinates.txt",
            col.names = T, row.names = T,
            quote = F,
            sep = "\t")
saveRDS(event.df.1base, "SUPPA/gencode.v26.splicing_events_coordinates.rds")

print("Creating bed file")
events.info <- event.df.1base
# Create bed file ----
# For events that are not SE or RI, there are 2 entries per event
# In theory the spliced-out version should be displayed in light grey
make.bed <- function(event.id){
  chr <- events.info[event.id,"chr"]
  start <- events.info[event.id,"spliced.in.start"] - 1
  end <- events.info[event.id,"spliced.in.end"]
  name <- events.info[event.id,"event.label"]
  score <- 1000
  strand <- events.info[event.id,"strand"]
  if(events.info[event.id,"type"] %in% c("SE","RI", "A5", "A3")){
    df <- t(as.data.frame(c(chr, start, end, name, score, strand))) 
    rownames(df) <- NULL
    return(df)
  }else{
    start2 <- events.info[event.id,"spliced.out.start"] - 1
    end2 <- events.info[event.id,"spliced.out.end"]
    score2 <- 100
    df <- rbind.data.frame( t(as.data.frame(c(chr, start, end, name, score, strand))),
                            t(as.data.frame(c(chr, start2, end2, name, score2, strand)))
    )
    rownames(df) <- NULL
    return(df)
  }
}

# BED file --
events.bed <- do.call(rbind.data.frame,
                      lapply(rownames(events.info), function(event.id) make.bed(event.id)))
lines <- nrow(events.info[!events.info$type %in% c("SE","RI", "A5", "A3"),]) * 2 + nrow(events.info[events.info$type %in% c("SE","RI","A5", "A3"),])
nrow(events.bed)
write.table(events.bed,
            "SUPPA/gencode.v26.splicing_events_coordinates.bed",
            col.names = F, row.names = F, quote = F, sep = "\t")
head(events.bed)

# ---------------------- #
end_time <- Sys.time()
# ---------------------- #

# ---------------------- #
print("")
print("# ---- Elapsed time ---- #")
end_time - start_time
print("# ---------------------- #\n")
# ---------------------- #



