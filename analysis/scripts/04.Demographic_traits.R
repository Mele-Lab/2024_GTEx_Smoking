#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to study the common effects of smoking and other demographic traits
# @software version: R=4.2.2

Sys.time()
#Set path 
# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")

#Load libraries 
library(reshape2)


# Reading data
tissues <- list.dirs("tissues/", full.names = F)[-1]
tissue_info <- read.csv("data/public/tissue_abreviation.txt")

#Saving tissue info sorted by sample size
tissue_info <- tissue_info[tissue_info$tissue %in% tissues,]
metadata <- readRDS("output/metadata.rds")
n_samples <- sapply(metadata, function(tissue) nrow(tissue))
tissue_info <- tissue_info[match(names(sort(n_samples)), tissue_info$tissue),]
write.csv(tissue_info, "data/public/tissues_sorted.csv", row.names = F)
tissues <- tissue_info$tissue

sex_tissues <- c("Uterus","Vagina","Ovary","Testis","Prostate")


dea_res <- list()
de.genes <- list()
expression_tissues <- c() #tissues with Smoking-DEGs
for(i in 1:length(tissues)){
  tissue <- tissues[i]
  #Reading differential expression results
  dea <- readRDS(paste0("tissues/", tissue, "/voom_limma_results.rds"))
  
  if(sum(dea$Smoking2$adj.P.Val<0.05)>0){ #Only tissues with smoking-DEGs. No, now we include all
    expression_tissues <- c(expression_tissues, tissue)
    }
  dea_res[[tissue]] <- dea
  de.genes[[tissue]] <- sapply(dea, function(trait) rownames(trait[trait$adj.P.Val<0.05,]))

}

#Save file with the tissues with expression
# write.table(expression_tissues, paste0(first_dir, "04_Smoking/Plots/Hier_Part/expression_signal.txt"), quote = F, row.names = F, col.names = F)


#Now I have the input data as follows: 
#dea_res is a list (one per tissue-disease) of many lists (one per trait)
#de.genes is a list  (one per tissue-disease) of many lists (the genes DEG per trait)


# Overlap between paiwise combinations of traits ----
exprs.genes <- lapply(expression_tissues, function(tissue) rownames(dea_res[[tissue]][["Age"]])) #It does not matter which trait to use here
names(exprs.genes) <- expression_tissues

# Fisher's exact test to see if there are more DEGs with 2 traits than expected --
overlap.enrichment.fun <- function(tissue, trait.1, trait.2){ #test tissue <- tissues[i], trait.1 <- "Age, trait.2 <- acronyms[i]
  if(tissue %in% sex_tissues & trait.1 == "Sex"){
    return(list("overlap" = NA,
                "odds.ratio" = NA,
                "p.value" = NA))
  }else{
    if(trait.1=="Sex"){
      trait.1<-"Sex2" #Changed
    } else if(trait.1=="Ancestry"){
      if("Ancestry2" %in% names(dea_res[[tissue]])){
        trait.1<-"Ancestry2"
      }else{
        trait.1<-"AncestryEUR"  #EUR vs AFR
      }
    }
    trait.2 <- paste0(trait.2, 2) # Smoking2 is never vs smokers
    x11 <- length(intersect(de.genes[[tissue]][[trait.1]], de.genes[[tissue]][[trait.2]]))
    x12 <- length(de.genes[[tissue]][[trait.1]][! de.genes[[tissue]][[trait.1]] %in% de.genes[[tissue]][[trait.2]]])
    x21 <- length(de.genes[[tissue]][[trait.2]][! de.genes[[tissue]][[trait.2]] %in% de.genes[[tissue]][[trait.1]]])
    x22 <- length(exprs.genes[[tissue]][! exprs.genes[[tissue]] %in% unique(c(de.genes[[tissue]][[trait.1]], de.genes[[tissue]][[trait.2]]))])
    m <- matrix(c(x11,x12,x21,x22),2,2,byrow = T)
    rownames(m) <- c(paste0(trait.1,".deg"), paste0(trait.1, ".not_deg"))
    colnames(m) <- c(paste0(trait.2,".deg"), paste0(trait.2, ".not_deg"))
    f <- fisher.test(m, alternative = "greater")
    f$estimate
    return(list("overlap" = x11,
                "odds.ratio" = f$estimate,
                "p.value" = f$p.value))#,
    # "counts.matrix" = m,
    # "lower_CI" = f$conf.int[1],
    # "upper_CI" = f$conf.int[2]))
  }
}

# Pairwise combination of traits --
pw.traits <- c("Age & Smoking", "Ancestry & Smoking", "Sex & Smoking", "BMI & Smoking") 
variables <- c("Age", "Ancestry", "Sex", "BMI")

# Enrichment analysis  --
fisher.results <- lapply(variables, function(variable)
  lapply(1:length(expression_tissues), function(i) 
    overlap.enrichment.fun(expression_tissues[i], variable, "Smoking")
  ))
names(fisher.results) <- variables
for(pw in variables){names(fisher.results[[pw]]) <- expression_tissues}

# Enrichment statistics --
p.value <- sapply(variables, function(pw)
  sapply(expression_tissues, function(tissue) 
    fisher.results[[pw]][[tissue]][["p.value"]]
  )
)

#MUltiple testing correction by demographic trait
fdr <- apply(p.value, 2, function(x) p.adjust(x, method="BH")) #Nas in Sex are not considered in multiple testing
fdr[fdr>0.05] <- NA

# Number of DEGs in overlap --
overlap <- sapply(variables, function(pw)
  sapply(expression_tissues, function(tissue) 
    fisher.results[[pw]][[tissue]][["overlap"]]
  )
)
rownames(overlap) <- expression_tissues

# Odds ratio --
odds.ratio <- sapply(variables, function(pw)
  sapply(expression_tissues, function(tissue) 
    fisher.results[[pw]][[tissue]][["odds.ratio"]]
  )
)
rownames(odds.ratio) <- expression_tissues
or <- odds.ratio
or[is.na(fdr)] <- NA

fdr[is.na(fdr)] <- 1
# to_plot <- log2(or)
to_plot <- or
to_plot[fdr>0.05] <- 0
to_plot[to_plot=="-Inf"] <- 0
to_plot_overlap <- to_plot


### Interaction plot

interactions <- matrix(nrow = length(expression_tissues), ncol=4, dimnames = list(expression_tissues, c("Age", "Angestry", "Sex", "BMI")))

for(tissue in expression_tissues){ 
  print(tissue)
  file <- paste0("tissues/", tissue, "/interactions.rds")
  if(!file.exists(file)){
    print("No interactions are tested")
    next}
  model <- readRDS(file)
  if(!is.null(model$`Age:Smoking2`)){
    interactions[tissue,1] <- sum(model$`Age:Smoking2`$adj.P.Val<0.05)}
  if(!is.null(model$`AncestryEUR:Smoking2`)){
    interactions[tissue,2] <- sum(model$`AncestryEUR:Smoking2`$adj.P.Val<0.05)}
  if(!is.null(model$`Sex2:Smoking2`)){
    interactions[tissue,3] <- sum(model$`Sex2:Smoking2`$adj.P.Val<0.05)}
  if(!is.null(model$`BMI:Smoking2`)){
    interactions[tissue,4] <- sum(model$`BMI:Smoking2`$adj.P.Val<0.05)}
}



# png(paste0(first_dir, "/04_Smoking/Plots/Interactions.png"),
#     units="in", width = 3.5, height = 9, res=200)
# create_heatmap(interactions, tissue_info, tissues_plot = rownames(interactions), diseases_plot = colnames(interactions),
#                label_data=interactions, color_palette = colorRamp2( c(0,1,20), brewer.pal(8, "BuPu")[c(1,4,7)]))
# dev.off()


### Bias Plot
bias.fun <- function(tissue, trait1, trait2){
  print(paste0(tissue, ": ", trait1, "-", trait2))
  if(tissue %in% sex_tissues & trait1 == "Sex"){
    return(list("p-value" = NA, "which" = NA,
                "OvsE" = "NA:NA:NA:NA", "observed" = NA,
                "enough" = F
    ))
  }else{
    if(trait1=="Sex"){
      trait1<-"Sex2"
    } else if(trait1=="Ancestry"){
      trait1<-"AncestryEUR"
    }
    trait2 <- paste0(trait2, 2)
    # Do we observe a higher than expected overlap of DEGs in a particular direction of change?
    genes.overlap <- intersect(de.genes[[tissue]][[trait1]], de.genes[[tissue]][[trait2]])
    trait1.up <- rownames(dea_res[[tissue]][[trait1]][dea_res[[tissue]][[trait1]]$adj.P.Val < 0.05 &
                                                        dea_res[[tissue]][[trait1]]$logFC > 0,])
    trait1.down <- rownames(dea_res[[tissue]][[trait1]][dea_res[[tissue]][[trait1]]$adj.P.Val < 0.05 &
                                                          dea_res[[tissue]][[trait1]]$logFC < 0,])
    
    trait2.up <- rownames(dea_res[[tissue]][[trait2]][dea_res[[tissue]][[trait2]]$adj.P.Val < 0.05 &
                                                        dea_res[[tissue]][[trait2]]$logFC > 0,])
    trait2.down <- rownames(dea_res[[tissue]][[trait2]][dea_res[[tissue]][[trait2]]$adj.P.Val < 0.05 &
                                                          dea_res[[tissue]][[trait2]]$logFC < 0,])
    
    # Observed counts
    counts <- c(sum(trait1.up %in% trait2.up), # upup
                sum(trait1.down %in% trait2.up), # downup
                sum(trait1.up %in% trait2.down), # updown
                sum(trait1.down %in% trait2.down) # downdown
    )
    # Expected proportions
    trait1.up.p <- length(trait1.up)/length(c(trait1.up, trait1.down))
    trait1.down.p <- length(trait1.down)/length(c(trait1.up, trait1.down))
    trait2.up.p <- length(trait2.up)/length(c(trait2.up, trait2.down))
    trait2.down.p <- length(trait2.down)/length(c(trait2.up, trait2.down))
    expected_prob <- c(trait1.up.p * trait2.up.p, trait1.down.p * trait2.up.p, trait1.up.p * trait2.down.p, trait1.down.p * trait2.down.p)
    expected_counts <- expected_prob*length(genes.overlap)
    print(expected_counts)
    if(length(genes.overlap)==0){
      return(list("p-value" = NA, "which" = NA,
                  "OvsE" = "NA:NA:NA:NA",
                  "observed" = paste(as.character(counts), collapse = "/"),
                  "enough" = F
      ))
      # } else if(length(genes.overlap)<20){
    } else if(length(genes.overlap)<10){
      return(list("p-value" = NA, "which" = NA,
                  "OvsE" = "NA:NA:NA:NA",
                  "observed" = paste(as.character(counts), collapse = "/"),
                  "enough" = F
      ))
    }
    
    # Chi-square goodness of fit test
    if(min(expected_counts)<5){
      chi.test <- chisq.test(counts, p = expected_prob, simulate.p.value = T)
      enough <- F
    }else{
      chi.test <- chisq.test(counts, p = expected_prob)
      enough <- T
    }
    return(list("p-value" = chi.test$p.value,
                "which" = paste(c("up_up", "down_up", "up_down", "down_down")[chi.test$observed/round(chi.test$expected) > 1],collapse = ":"),
                "OvsE" = paste0(as.character(round(chi.test$observed/round(chi.test$expected),2)),collapse = ":"),
                "observed" = paste(as.character(counts), collapse = "/"),
                "enough" = enough
    ))
  }
}


# Enrichment analysis  --

chi_square.results <- lapply(variables, function(variable)
  lapply(1:length(tissues), function(i) 
    bias.fun(tissues[i], variable, "Smoking")
  ))
names(chi_square.results) <- variables
for(pw in variables){names(chi_square.results[[pw]]) <- tissues}


# Pairwise combination of traits --
pw.traits <- c("Age & Smoking", "Ancestry & Smoking", "Sex & Smoking", "BMI & Smoking") 
variables <- c("Age", "Ancestry", "Sex", "BMI")


observed <- sapply(variables, function(pw)
  sapply(tissues, function(tissue) 
    chi_square.results[[pw]][[tissue]][["observed"]]))

# p-values
p_values <- sapply(variables, function(pw)
  sapply(tissues, function(tissue) 
    chi_square.results[[pw]][[tissue]][["p-value"]]))


#FDR is computed when plotting


### Example gene plots
gene_annotation <- read.delim("data/public/gencode.v26.GRCh38.genes.bed", header=F)[,c(1,6,7)]
colnames(gene_annotation) <- c("chr", "gene", "symbol")

traits_cols <- c("Age" = "#56B4E9", "Ancestry" = "#E69F00", 
                 "Sex" = "#009E73", "BMI" = "#CC79A7", "Smoking" = "#000000") #696969
my_traits <- names(traits_cols)

get.gene.data <- function(acronym, tissue, gene, combo, g = 2){
  #  Expression data ----
  tpm <- readRDS(paste0("tissues/", tissue, "/tpm.rds")) #The plot with residuals is the same
  metadata <- readRDS(paste0("tissues/", tissue, "/metadata.rds"))
  
  #Excluding ex-smokers
  metadata <- metadata[metadata$Smoking!=1,]
  tpm <- tpm[,colnames(tpm) %in% metadata$Sample]  #Subset for the metadata we want
  
  #Subset with the gene data
  gene.name <- gene_annotation[gene_annotation$gene==gene, "symbol"]
  df <- as.data.frame(t(tpm[gene,]))
  colnames(df) <- "TPM"
  
  df$Age_int <- sapply(rownames(df), function(i) metadata[metadata$Sample==i,"Age"])
  df$BMI_int <- sapply(rownames(df), function(i) metadata[metadata$Sample==i,"BMI"])
  metadata$Age_Class <- cut(metadata$Age, breaks = c(20, 30, 40, 50, 60, 70), include.lowest=T, right=F) #Include_lowest allows 20 to be included
  levels(metadata$Age_Class) <- c("20-29", "30-39", "40-49", "50-59", "60-70")
  df$Age <- sapply(rownames(df), function(i) metadata[metadata$Sample==i,"Age_Class"])
  df$Age <- factor(df$Age, levels = c("20-29",
                                      "30-39",
                                      "40-49",
                                      "50-59",
                                      "60-70"),
                   order = T)
  
  df$Age2 <- ifelse(df$Age_int < 45, "Young", "Old")
  df$Age2 <- factor(df$Age2,
                    levels = c("Young", "Old"),
                    order = T)
  df$Ancestry <- sapply(rownames(df), function(i) metadata[metadata$Sample==i,"Ancestry"])
  df$Ancestry <- gsub("EUR", "EA", df$Ancestry)
  df$Ancestry <- gsub("AFR", "AA", df$Ancestry)
  df$Ancestry <- gsub("AMR", "AMR", df$Ancestry)
  df$Ancestry <- factor(df$Ancestry, levels = c("EA", "AA", "AMR"), order = T)
  df$Sex <- sapply(rownames(df), function(i) metadata[metadata$Sample==i,"Sex"])
  df$Sex <- gsub("1", "Male", df$Sex)
  df$Sex <- gsub("2", "Female", df$Sex)
  df$Sex <- factor(df$Sex, levels = c("Male", "Female"), order = T)
  metadata$BMI_Class <- cut(metadata$BMI, breaks = c(15, 25, 30, 36), include.lowest=T, right=F) #Include_lowest allows 20 to be included
  
  df$BMI <- sapply(rownames(df), function(i) metadata[metadata$Sample==i,"BMI_Class"])
  levels(df$BMI) <- c("Normal", "Overweight", "Obese")
  df$BMI2 <- ifelse(df$BMI_in < 25, "<25", ">=25")
  df$BMI2 <- factor(df$BMI2,
                    levels = c("<25", ">=25"),
                    order = T)
  df$Tissue <- rep(tissue, nrow(df))
  df$Gene <- rep(gene.name, nrow(df))
  
  df$Smoking <- sapply(rownames(df), function(i) metadata[metadata$Sample==i, acronym])
  df$Smoking <- as.character(df$Smoking)
  df$Smoking[df$Smoking==0] <- "Never smokers"
  df$Smoking[df$Smoking==2] <- "Smokers"
  df$Smoking <- factor(df$Smoking,
                       levels = c("Never smokers", "Smokers"),
                       order = T)
  
  if(combo == "Smoking:Age"){
    df$x_dummy <- paste0(df$Age2, "\n", df$Smoking)
    df$x_dummy <- factor(df$x_dummy,
                         levels = unlist(lapply(c("Young", "Old"),
                                                function(i) paste0(i, c("\nNever smokers", "\nSmokers")))),
                         order = T)  
  }else if(combo == "Smoking:Sex"){
    df$x_dummy <- paste0(df$Sex, "\n", df$Smoking)
    df$x_dummy <- factor(df$x_dummy,
                         levels = unlist(lapply(c("Male", "Female"),
                                                function(i) paste0(i, c("\nNever smokers", "\nSmokers")))),
                         order = T)  
  }else if(combo == "Smoking:Ancestry"){
    df$x_dummy <- paste0(df$Ancestry, "\n", df$Smoking)
    df$x_dummy <- factor(df$x_dummy,
                         levels = unlist(lapply(c("EA", "AA", "AMR"),
                                                function(i) paste0(i, c("\nNever smokers", "\nSmokers")))),
                         order = T)  
  }else if(combo == "Smoking:BMI"){
    df$x_dummy <- paste0(df$BMI2, "\n", df$Smoking)
    df$x_dummy <- factor(df$x_dummy,
                         levels = unlist(lapply(c("<25", ">=25"),
                                                function(i) paste0(i, c("\nNever smokers", "\nSmokers")))),
                         order = T)  
  }
  return(df)
}


acronym <- "Smoking"
tissue <- "Thyroid"
variables <- "Smoking:Age"
gene <- "ENSG00000113721.13" #PDGFRB

data_PDGFRB <- get.gene.data(acronym, tissue, gene, variables, g = 2)
# p1 <- my.box_plot(data, gene, "Smoking", "Age2")

#I may want to find another lung example with age, maybe a known one
tissue <- "Lung"
variables <- "Smoking:Age"
gene <- "ENSG00000244694.7" #PTCHD4
data_PTCHD4 <- get.gene.data(acronym, tissue, gene, variables, g = 2)
# p2 <- my.box_plot(data, gene, "Smoking", "Age2")

# png(paste0(first_dir, "Jose/04_Smoking/Plots/Overlaps/", tissue, "_", gene_annotation[gene_annotation$gene == gene, "symbol"], ".png"),
#       units = "in", width = 4, height = 3.5, res=150)
# print(p1)
# dev.off()

#Save data here:
figure_2 <- list(to_plot_overlap, tissue_info, overlap, chi_square.results, tissues, p_values, data_PDGFRB, data_PTCHD4, gene_annotation)
names(figure_2) <- c("to_plot_overlap", "tissue_info", "overlap", "chi_square.results", "tissues", "p_values", "data_PDGFRB", "data_PTCHD4", "gene_annotation")
saveRDS(figure_2, "../figures/data/figure_2.rds")





### Hier Part
autosomal_genes <- gene_annotation[!gene_annotation$chr %in% c("chrX","chrY","chrM"),"gene"] 

# Hier.part:expression ----
# All genes expressed in tissue 
hier.part.exprs <- list()
for(i in 1:length(tissues)){ #iterating over tissues
  data <- readRDS(paste0("tissues/", tissues[i], "/hier.part.rds"))
  hier.part.exprs[[tissues[i]]] <- data
}

dea_res_autosomal <- dea_res
# Keep only autosomal genes
for(i in 1:length(tissues)){
  name <- tissues[i]
  for(trait in names(dea_res_autosomal[[name]])){
    dea_res_autosomal[[name]][[trait]] <- dea_res_autosomal[[name]][[trait]][rownames(dea_res_autosomal[[name]][[trait]]) %in% autosomal_genes,] 
  }
}

#Barplots

# 4.3 Create gene:trait adj.P.Val matrix ----
fdr.matrix <- lapply(1:length(tissues), function(i)
  do.call(cbind.data.frame,
          lapply(names(dea_res_autosomal[[tissues[i]]]), function(trait)
            sapply(rownames(hier.part.exprs[[tissues[i]]]), function(gene)
              dea_res_autosomal[[tissues[i]]][[trait]][gene,"adj.P.Val"]
            )
          )
  )
)

names(fdr.matrix) <- tissues
for(i in 1:length(tissues)){
  colnames(fdr.matrix[[tissues[i]]]) <- names(dea_res[[tissues[i]]])
  fdr.matrix[[tissues[i]]][fdr.matrix[[tissues[i]]] >= 0.05] <- NA
  hier.part.exprs[[tissues[i]]] <- hier.part.exprs[[tissues[i]]][,grepl("_abs", names(hier.part.exprs[[tissues[i]]]))]
}


#Ancestry if signif in any comparison and smoking only smoking2
for(i in 1:length(tissues)){ #take into account sex, but not always the variables are in order, osmetimes one ancestry is missing, or age is missing
  print(tissues[i])
  boolean <- is.na(fdr.matrix[[tissues[i]]])
  if(tissues[i] %in% sex_tissues){
    final_boolean <- boolean[, which(colnames(boolean) %in% c("Age", "AncestryEUR", "BMI", "Smoking2"))] #Smoking will be Smoking2 for now
    if(sum(grepl("Ancestry", colnames(boolean)))==3){
      tmp_boolean <- boolean[,which(colnames(boolean) == c("AncestryEUR"))] & boolean[,which(colnames(boolean) == c("AncestryAMR"))] & boolean[,which(colnames(boolean) == c("AncestryAMR-AncestryEUR"))] #Is the gene affected by Ancestry (at least one comparsion)?
    } else{ #In one case, Pituitary, we only have two ancestries, but I do it like this in case other examples with only 2 ancestries happened
      first <- grep("Ancestry", colnames(boolean))[1]
      second <- grep("Ancestry", colnames(boolean))[2]
      tmp_boolean <- boolean[,first] & boolean[,second] 
    }
    final_boolean[,"AncestryEUR"] <- tmp_boolean #Ancestry
    
    # #Try smoking as ancestry (three levels)
    # tmp_boolean <- boolean[,which(colnames(boolean) == c("Smoking"))] & boolean[,which(colnames(boolean) == c("Smoking2"))] & boolean[,which(colnames(boolean) == c("Smoking1-Smoking2"))] #Is the gene affected by Ancestry (at least one comparsion)?
    # final_boolean[,"Smoking2"] <- tmp_boolean
    
  } else{
    final_boolean <- boolean[,which(colnames(boolean) %in% c("Age", "AncestryEUR", "Sex2", "BMI", "Smoking2"))] #Smoking will be Smoking2 only, 
    if(sum(grepl("Ancestry", colnames(boolean)))==3){
      tmp_boolean <- boolean[,which(colnames(boolean) == c("AncestryEUR"))] & boolean[,which(colnames(boolean) == c("AncestryAMR"))] & boolean[,which(colnames(boolean) == c("AncestryAMR-AncestryEUR"))] #Is the gene affected by Ancestry (at least one comparsion)?
    } else{ #In one case, Pituitary, we only have two ancestries, but I do it like this in case other examples with only 2 ancestries happened
      print(sum(grepl("Ancestry", colnames(boolean))))
      first <- grep("Ancestry", colnames(boolean))[1]
      second <- grep("Ancestry", colnames(boolean))[2]
      tmp_boolean <- boolean[,first] & boolean[,second] 
    }   
    final_boolean[,"AncestryEUR"] <- tmp_boolean #Even though the name is Ancestry EUR (EUR vs AFR), this variable now represents all Ancestries
    #Ancestry boolean will be TRUE if the gene is significant in any comparison (EUR-AFR OR AMR-AFR OR AMR-EUR) 
    
    # #Try smoking as ancestry (three levels with ex-smokers), but it didn't change the results and it makes more sense to only consider never vs smokers
    # tmp_boolean <- boolean[,which(colnames(boolean) == c("Smoking"))] & boolean[,which(colnames(boolean) == c("Smoking2"))] & boolean[,which(colnames(boolean) == c("Smoking1-Smoking2"))] #Is the gene affected by Ancestry (at least one comparsion)?
    # final_boolean[,"Smoking2"] <- tmp_boolean
  }
  if(!"Age" %in% colnames(final_boolean)){ #Necessary for small intestine
    final_boolean <- cbind(rep(NA, nrow(final_boolean)), final_boolean)
    names(final_boolean)[1] <- "Age"
  }
  final_names <- gsub('_abs', '', names(hier.part.exprs[[tissues[i]]]))
  names(hier.part.exprs[[tissues[i]]]) <- final_names
  hier.part.exprs[[tissues[i]]][final_boolean] <- NA
}

# 4.5 Explained variance (ev) ----
get_expln_var <- function(tissue, hier.part.data, i){
  if(is.null(hier.part.data[[tissue]])){
    x <- rep(0, length(names(hier.part.data[[tissue]])))
    names(x) <- names(hier.part.data[[tissue]])
  } else{
    x <- colSums(hier.part.data[[tissue]], na.rm = T)
  }
  return(x)
}

# 4.6 Explained variance (ev) ; autosomal genes ----   #Not necessary, I already filtered by autosomal_genes in a previous step in the last update
hier.part.exprs.autosomal_genes <- lapply(1:length(tissues), function(i)
  hier.part.exprs[[tissues[i]]][rownames(hier.part.exprs[[tissues[i]]]) %in% autosomal_genes,])
names(hier.part.exprs.autosomal_genes) <- tissues

ev.exprs.total.aut_genes <- lapply(1:length(tissues), function(i)
  get_expln_var(tissues[i], hier.part.exprs.autosomal_genes, i)
)
names(ev.exprs.total.aut_genes) <- tissues


# 4.8 Bar plot ordered by sample size ----
traits <- c()
my_tissues <- c()
for(i in 1:length(tissues)){
  traits_t <- names(ev.exprs.total.aut_genes[[tissues[i]]])
  traits <- c(traits, traits_t)
  tissue <- strsplit(tissues[i],".", fixed = T)[[1]][1]
  my_tissues <- c(my_tissues, rep(tissue, length(traits_t)))
}

data_plot <- cbind.data.frame(traits,
                              my_tissues,
                              unlist(ev.exprs.total.aut_genes)) 
colnames(data_plot) <- c("Trait","Tissue", "Value")                 
data_plot$Abbrv <- sapply(data_plot$Tissue, function(i) tissue_info[tissue_info$tissue_ID==i,"tissue_abbrv"])
data_plot$Trait <- as.factor(data_plot$Trait)#, levels=traits, order = T)

data_plot$Trait_type <- data_plot$Trait
data_plot$Trait_type <- as.character(data_plot$Trait_type)
data_plot$Trait_type <- as.factor(data_plot$Trait_type)

final_data <- reshape(data_plot, idvar = "Tissue", timevar = "Trait_type", direction = "wide")
final_data <- final_data[,c(1, grep("Value", colnames(final_data)))]
rownames(final_data) <- final_data[,1]
final_data <- final_data[,-1]
colnames(final_data) <- unique(data_plot$Trait)
final_data <- final_data[,my_traits] #Sorting columns in the order we want -> Smoking last
final_data <- as.matrix(final_data)
final_data[is.na(final_data)] <- 0
final_data_frac <- final_data/rowSums(final_data)
final_data_frac <- final_data_frac[nrow(final_data_frac):1, ]

rownames(final_data_frac) <- sapply(rownames(final_data_frac), function(tissue) tissue_info$TISSUENAMEABREV[tissue_info$tissue==tissue])

order <- colnames(final_data_frac)

test <- final_data_frac
test_2 <- test[order(match(rownames(test),tissue_info$tissue_abbrv), decreasing=T),]
expression_barplot <- test_2


#Pirate Plots


# 5.1.1 Prepare expression data ----
melt.data <- function(i, data){
  tissue <- expression_tissues[i]
  print(tissue)
  df <- data[[tissue]]
  traits <- names(df)
  df$ensembl.id <- rownames(df)
  d <- reshape2::melt(df[,c(traits,"ensembl.id")],
                      id.vars = "ensembl.id",
                      variable.name = 'Trait',
                      value.name = 'R2')
  
  # d$Trait <- gsub(pattern = "_abs",replacement = "", d$Trait)
  tissue <- strsplit(tissue, ".", fixed=T)[[1]][1]
  d$Tissue <- rep(tissue, nrow(d))
  d$name <- rep(tissue, nrow(d))
  return(d)
}

#Maybe here we should use expression_tissues instead of tissues and not autosomal genes

exprs.data <- do.call(rbind.data.frame,
                      lapply(1:length(expression_tissues), function(i) 
                        melt.data(i, hier.part.exprs)
                        # melt.data(i, hier.part.exprs.autosomal_genes)
                      ))
exprs.data <- exprs.data[!is.na(exprs.data$R2),]
exprs.data <- exprs.data[!duplicated(exprs.data),]

exprs.data$Tissue <- factor(exprs.data$Tissue,
                            levels = unique(exprs.data$Tissue),
                            order = T)

exprs.data_disease <- exprs.data[!exprs.data$Trait %in% c("Age", "BMI", "Ancestry", "Sex"),]

exprs.data_disease$Trait <- as.factor(exprs.data_disease$Trait)

exprs.data_disease$Tissue_abbvr <- sapply(exprs.data_disease$Tissue, function(x) {
  tissue_info$TISSUENAMEABREV[tissue_info$tissue==x]
})
# exprs.data_disease$Tissue_abbvr <- factor(exprs.data_disease$Tissue_abbvr, levels = rev(rownames(test_2)))
exprs.data_disease$Tissue_abbvr <- factor(exprs.data_disease$Tissue_abbvr, levels = unique(exprs.data_disease$Tissue_abbvr))
exprs.data_disease$R2 <- as.numeric(exprs.data_disease$R2)
expression_pirate <- exprs.data_disease





### Hier Part Splicing
tissues <- list.dirs("tissues/", full.names = F)[-1]
# Differential splicing results ----
names <- c()
dsa_res <- list()
for(i in 1:length(tissues)){
  tissue <- tissues[i]
  print(tissue)
  file <- list.files(paste0("tissues/", tissue, "/"))[grepl("fractional_regression_results.rds",list.files(paste0("tissues/", tissue, "/")))] #We do this in case there is no file for splicing
  if(length(file)==0){next}
  data <- readRDS(paste0("tissues/", tissue, "/", file))
  # if(sum(data$Smoking2$adj.P.Val<0.05)>0){ #Only tissues with smoking-DSEs
  dsa_res[[tissue]] <- data
  names <- c(names, tissue)
  # }
}

hier.part.splicing <- list()
for(i in 1:length(names)){
  data <- readRDS(paste0("tissues/", names[i], "/Alternatively_spliced_events.hier_part.rds"))
  hier.part.splicing[[names[i]]] <- data
}

#Get DSGs
function_for_dsg <- function(i, trait){
  genes <- dsa_res[[names[i]]][[trait]][dsa_res[[names[i]]][[trait]]$adj.P.Val<0.05,]$Ensembl_id
  genes <- genes[genes %in% autosomal_genes]
  return(genes)
}

ds.genes <- lapply(1:length(names), function(i)
  lapply(names(dsa_res[[names[i]]]), function(trait)
    function_for_dsg(i, trait)
  )
)
names(ds.genes) <- names

for(tissue in names){
  names(ds.genes[[tissue]]) <- names(dsa_res[[tissue]])
}

#Barplot
get.dsg <- function(i){
  genes <- unique(unlist(ds.genes[[names[i]]]))
  return(genes)
}

dsg <- lapply(1:length(names), function(i) get.dsg(i))
names(dsg) <- names

# 4.2 Subset hier.part data ---- This has already been done in the previous script
my_toy_function <- function(hier.part.splicing, tissue, genes){
  if(is.null(hier.part.splicing[[tissue]])){
    return(NULL)
  } else{
    return(hier.part.splicing[[tissue]][unlist(strsplit(rownames(hier.part.splicing[[tissue]]), ";"))[ c(TRUE,FALSE) ] %in% genes,])
  }
}
hier.part.splicing.autosomal_genes <- lapply(names, function(tissue) my_toy_function(hier.part.splicing, tissue, autosomal_genes))
names(hier.part.splicing.autosomal_genes) <- names

# 4.3 Create gene:trait adj.P.Val matrix ----
fdr.matrix_s <- lapply(1:length(names), function(i)
  do.call(cbind.data.frame,
          lapply(names(dsa_res[[names[i]]]), function(trait)
            sapply(rownames(hier.part.splicing.autosomal_genes[[names[i]]]), function(gene)
              dsa_res[[names[i]]][[trait]][gene,"adj.P.Val"]
            )
          )
  )
)
names(fdr.matrix_s) <- names
for(i in 1:length(names)){
  colnames(fdr.matrix_s[[names[i]]]) <- names(dsa_res[[names[i]]])
}

for(i  in 1:length(names)){
  fdr.matrix_s[[names[i]]][fdr.matrix_s[[names[i]]] >= 0.05] <- NA
}


# 4.4 Only consider R2 of DE genes -----
for(i in 1:length(names)){ #keeping only the absolute R2 values
  hier.part.splicing.autosomal_genes[[names[i]]] <- hier.part.splicing.autosomal_genes[[names[i]]][,grepl("_abs", names(hier.part.splicing.autosomal_genes[[names[i]]]))]
}

for(i in 1:length(names)){
  print(names[i])
  boolean <- is.na(fdr.matrix_s[[names[i]]])
  
  if(names[i] %in% sex_tissues){
    final_boolean <- boolean[, which(colnames(boolean) %in% c("Age", "AncestryEUR", "BMI", "Smoking2"))] #Smoking will be Smoking2 for now
    if(sum(grepl("Ancestry", colnames(boolean)))==3){
      tmp_boolean <- boolean[,which(colnames(boolean) == c("AncestryEUR"))] & boolean[,which(colnames(boolean) == c("AncestryAMR"))] & boolean[,which(colnames(boolean) == c("AncestryAMR-EUR"))] #Is the gene affected by Ancestry (at least one comparsion)?
    } else{ #In one case, Pituitary, we only have two ancestries, but I do it like this in case other examples with only 2 ancestries happened
      first <- grep("Ancestry", colnames(boolean))[1]
      second <- grep("Ancestry", colnames(boolean))[2]
      tmp_boolean <- boolean[,first] & boolean[,second] 
    }
    final_boolean[,"AncestryEUR"] <- tmp_boolean #Ancestry
    final_boolean <- final_boolean[,c(1,4,2,3)]
  } else{
    final_boolean <- boolean[,which(colnames(boolean) %in% c("Age", "AncestryEUR", "Sex2", "BMI", "Smoking2"))] #Smoking will be Smoking2 only
    if(sum(grepl("Ancestry", colnames(boolean)))==3){
      tmp_boolean <- boolean[,which(colnames(boolean) == c("AncestryEUR"))] & boolean[,which(colnames(boolean) == c("AncestryAMR"))] & boolean[,which(colnames(boolean) == c("AncestryAMR-EUR"))] #Is the gene affected by Ancestry (at least one comparsion)?
    } else{ #In some cases we only have two ancestries
      first <- grep("Ancestry", colnames(boolean))[1]
      second <- grep("Ancestry", colnames(boolean))[2]
      tmp_boolean <- boolean[,first] & boolean[,second] 
    }   
    final_boolean[,"AncestryEUR"] <- tmp_boolean
    final_boolean <- final_boolean[,c(1,5,2,3,4)]
  }
  if(!"Age" %in% colnames(final_boolean)){ #Necessary for small intestine
    final_boolean <- cbind(rep(NA, nrow(final_boolean)), final_boolean)
    names(final_boolean)[1] <- "Age"
  }
  #Ancestry boolean will be TRUE if the gene is significant in any comparison (EUR-AFR OR AMR-AFR OR AMR-EUR) 
  final_names <- gsub('_abs', '', names(hier.part.splicing.autosomal_genes[[names[i]]]))
  names(hier.part.splicing.autosomal_genes[[names[i]]]) <- final_names
  hier.part.splicing.autosomal_genes[[names[i]]][final_boolean] <- NA
}

# 4.6 Explained variance (ev) ; autosomal genes ----
ev.splicing.total.aut_genes <- lapply(1:length(names), function(i)
  get_expln_var(names[i], hier.part.splicing.autosomal_genes, i)
)
names(ev.splicing.total.aut_genes) <- names

# 4.8 Bar plot order by sample size ----
traits <- c()
my_tissues <- c()
# folder_names <- c()
for(i in 1:length(names)){
  names_t <- names(ev.splicing.total.aut_genes[[names[i]]])
  traits <- c(traits, names_t)
  tissue <- strsplit(names[i],".", fixed = T)[[1]][1]
  my_tissues <- c(my_tissues, rep(tissue, length(names_t)))
}


data_plot <- cbind.data.frame(traits,
                              my_tissues,
                              unlist(ev.splicing.total.aut_genes))   #The fraction (frac.) is not correctly computed for the diseases in the lung
colnames(data_plot) <- c("Trait","Tissue", "Value")                 
data_plot$Abbrv <- sapply(data_plot$Tissue, function(i) tissue_info[tissue_info$tissue==i,"TISSUENAMEABREV"])
data_plot$Trait <- as.factor(data_plot$Trait)#, levels=traits, order = T)

data_plot$Trait_type <- data_plot$Trait
data_plot$Trait_type <- as.character(data_plot$Trait_type)
data_plot$Trait_type <- as.factor(data_plot$Trait_type)

final_data <- reshape(data_plot, idvar = "Tissue", timevar = "Trait_type", direction = "wide")
final_data <- final_data[,c(1, grep("Value", colnames(final_data)))]
rownames(final_data) <- final_data[,1]
final_data <- final_data[,-1]
colnames(final_data) <- unique(data_plot$Trait)
final_data <- final_data[,my_traits] #Sorting columns in the order we want -> Smoking last
final_data <- as.matrix(final_data)
final_data[is.na(final_data)] <- 0
final_data_frac <- final_data/rowSums(final_data)
final_data_frac <- final_data_frac[nrow(final_data_frac):1, ]

# order <- colnames(final_data_frac)

final_data_frac <- final_data/rowSums(final_data)
final_data_frac <- final_data_frac[nrow(final_data_frac):1, ]

rownames(final_data_frac) <- sapply(rownames(final_data_frac), function(tissue) tissue_info$TISSUENAMEABREV[tissue_info$tissue==tissue])

test <- final_data_frac
test_2 <- test[order(match(rownames(test),tissue_info$TISSUENAMEABREV), decreasing=T),]

#Add zeros in the other tissues where expression has a signal -> NO
tissues <- tissues[tissues!="KidneyCortex"]
expression_tissues <- sapply(tissues, function(tissue) tissue_info$TISSUENAMEABREV[tissue_info$tissue==tissue])
test_3 <- matrix(nrow = length(expression_tissues), ncol=5, dimnames = list(expression_tissues, colnames(test_2)))

for(i in 1:nrow(test_3)){
  if(rownames(test_3)[i] %in% rownames(test_2)){
    test_3[i,] <- test_2[which(rownames(test_2)==rownames(test_3)[i]),]
  }
}
test_3 <- test_3[order(match(rownames(test_3),tissue_info$tissue_abbrv), decreasing=T),]
splicing_barplot <- test_3





#Pirate plot
melt.data <- function(i, data){
  tissue <- names(expression_tissues)[i]
  print(tissue)
  df <- data[[tissue]]
  traits <- names(df)
  df$ensembl.id <- rownames(df)
  d <- reshape2::melt(df[,c(traits,"ensembl.id")],
                      id.vars = "ensembl.id",
                      variable.name = 'Trait',
                      value.name = 'R2')
  
  # d$Trait <- gsub(pattern = "_abs",replacement = "", d$Trait)
  tissue <- strsplit(tissue, ".", fixed=T)[[1]][1]
  d$Tissue <- rep(tissue, nrow(d))
  d$name <- rep(tissue, nrow(d))
  return(d)
}

#expression_tissues instead of names. And all genes and not only autosomal
splicing.data <- do.call(rbind.data.frame,
                         lapply(1:length(expression_tissues), function(i) 
                           melt.data(i, hier.part.splicing.autosomal_genes)
                           # melt.data(i, hier.part.splicing)
                         ))
splicing.data <- splicing.data[!is.na(splicing.data$R2),]
splicing.data <- splicing.data[!duplicated(splicing.data),]


splicing.data$Tissue <- factor(splicing.data$Tissue,
                               levels = unique(splicing.data$Tissue),
                               order = T)

#Only smoking
splicing.data_disease <- splicing.data[!splicing.data$Trait %in% c("Age", "BMI", "Ancestry", "Sex"),]

splicing.data_disease$Trait <- as.factor(splicing.data_disease$Trait)

splicing.data_disease$Tissue_abbvr <- sapply(splicing.data_disease$Tissue, function(x) {
  tissue_info$TISSUENAMEABREV[tissue_info$tissue==x]
})

levels(splicing.data_disease$Tissue) <- expression_tissues
#Do the plot with all the tissues as well, so we can easily compare with expression
for(tissue in levels(expression_pirate$Tissue_abbvr)){
  if(!tissue %in% splicing.data_disease$Tissue_abbvr){
    splicing.data_disease <- rbind(splicing.data_disease, c("NA", "Smoking", 0, tissue, tissue, tissue))
  }
}
# tissues_to_plot <- tissue_info$TISSUENAMEABREV[tissue_info$tissue %in% splicing.data_disease$Tissue | tissue_info$TISSUENAMEABREV %in% rownames(test_3)]
# splicing.data_disease$Tissue_abbvr <- factor(splicing.data_disease$Tissue_abbvr, levels = tissues_to_plot)
splicing.data_disease$Tissue_abbvr <- factor(splicing.data_disease$Tissue_abbvr, levels = levels(expression_pirate$Tissue_abbvr))
splicing.data_disease$R2 <- as.numeric(splicing.data_disease$R2)
splicing_pirate <- splicing.data_disease

hier_part <- list(interactions, expression_barplot, expression_pirate, splicing_barplot, splicing_pirate, gene_annotation)
names(hier_part) <- c("interactions", "expression_barplot", "expression_pirate", "splicing_barplot", "splicing_pirate", "gene_annotation")
saveRDS(hier_part, "../figures/data/figure_S3.rds")
