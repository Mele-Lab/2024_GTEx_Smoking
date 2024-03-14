#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez and Raquel Garcia-Perez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to run residuals on gene expression
# @software version: R=4.2.2
rm(list=ls())

# Load libraries ####
suppressMessages(library(edgeR))
library(hier.part)
library(ggplot2)
library(ggpubr)

#Set path 
setwd(system("pwd", intern = T)) #If in linux
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
# setwd("..")

# -------------- #
print(Sys.time())
#-------------- #

tissues <- list.dirs("tissues/", full.names = F)[-1]
# tissues <- tissues[tissues!="KidneyCortex"]
# tissues <- tissues[33:length(tissues)]

#Functions
current.model.mod <- function (y, current.comb, xcan, SStot=0,family = c("gaussian","quasibinomial"), 
                               link = c("logit"), gof = c("Rsqu","RMSPE"), ...)  {
  comb.data <- data.frame(xcan[, current.comb])
  colnames(comb.data) <- colnames(xcan)[current.comb]
  data <- data.frame(y, comb.data)
  depv <- names(data)[1]
  n.comb <- dim(comb.data)[2]
  xs <- vector("character", n.comb)
  for (i in 1:(n.comb - 1)) xs[i] <- paste(names(comb.data)[i], 
                                           "+", sep = "")
  xs[n.comb] <- names(comb.data)[n.comb]
  xss <- paste(xs, collapse = " ", sep = "")
  formu <- stats::formula(paste(depv, "~", xss, sep = ""))
  if (gof == "RMSPE") gf <- sqrt(sum((stats::glm(formu, data = data, family = family,...)$fitted.values - y)^2))
  if (gof == "Rsqu") {
    if (family == "quasibinomial") 
      gf <- (SStot-sum((stats::glm(formu, data = data, family = family,...)$fitted.values - y)^2))/SStot
    if (family == "gaussian") 
      gf <- summary(stats::lm(formu, data = data))$r.squared
  }
  gf
}
all.regs.mod <- function (y, xcan, family = c("gaussian", "quasibinomial"), link = c("logit"), gof = c("Rsqu","RMSPE"),...) { 
  if (sum(is.na(xcan)) > 0) {
    missing <- is.na(apply(xcan, 1, FUN = sum))
    xcan <- xcan[!missing, ]
    y <- y[!missing]
    warning(paste(sum(missing), "observations deleted due to missingness in xcan\n"), 
            call. = FALSE)
  }
  if (sum(is.na(y)) > 0) {
    missing <- is.na(y)
    xcan <- xcan[!missing, ]
    y <- y[!missing]
    warning(paste(sum(missing), "observations deleted due to missingness in y\n"), 
            call. = FALSE)
  }
  pcan <- dim(xcan)[2]
  n <- (2^pcan) - 1
  combs <- combos1(pcan)$ragged
  SStot <- sum((y-mean(y))^2)
  
  if (gof == "RMSPE")  gfs <- sqrt(sum((stats::glm(y ~ 1, family = family,...)$fitted.values - y)^2))
  if (gof == "Rsqu")   gfs <- 0
  
  for (i in 1:n) {
    if (i%%500 == 0) 
      cat(i, "regressions calculated:", n - i, "to go...\n")
    current.comb <- as.vector(combs[i, ][combs[i, ] > 0]) 
    combn <- paste(names(data.frame(xcan)[current.comb]), "", collapse = "")
    # if (gof == "RMSPE") new.line <- current.model.mod(y, current.comb, xcan,family = family, gof = "RMSPE",...)
    # if (gof == "Rsqu")  new.line <- current.model.mod(y, current.comb, xcan,family = family, SStot=SStot,gof = "Rsqu",...)
    if (gof == "RMSPE") new.line <- current.model.mod(y, current.comb, xcan,family = family, gof = "RMSPE",...)
    if (gof == "Rsqu")  new.line <- current.model.mod(y, current.comb, xcan,family = family, SStot=SStot,gof = "Rsqu",...)
    gfs <- c(gfs, new.line)
  }
  gfs
}
hier.part.mod <- function(y,xcan,family='gaussian',gof = "Rsqu", link = "",...) {
  pcan <- dim(xcan)[2] #5
  gfs <- all.regs.mod(y, xcan, family = family, gof = gof, link = link, ...)
  hp <- partition(gfs, pcan, var.names = names(data.frame(xcan)))
  
  params <- list(full.model = paste("y ~", paste(names(xcan),collapse = " + ")), 
                 family = family, link = link, gof = gof)
  list(gfs = gfs, IJ = hp$IJ, I.perc = hp$I.perc, params = params)
}



sexual_tissues <- c("Prostate", "Testis", "Vagina", "Uterus", "Ovary")
#Do a for loop 
for(i in 1:length(tissues)){
  tissue <- tissues[i]
  print(paste0("Computing hier part for ", tissue))
  metadata <- readRDS(paste0("tissues/", tissue, "/metadata.rds"))
  if(tissue %in% sexual_tissues){
    metadata <- metadata[,names(metadata)!="Sex"]
  }
  
  if(sum(metadata$Ancestry=="AMR")==0){
    metadata$Ancestry <- droplevels(metadata$Ancestry)
  }
  
  dea_res <- readRDS(paste0("tissues/", tissue, "/voom_limma_results.rds"))

  exprs_genes <- rownames(dea_res$Age)
  # if(is.null(exprs_genes)){
  #   exprs_genes <- rownames(dea_res$BMI)
  # }
  
  # Counts
  counts <- readRDS(paste0("tissues/", tissue, "/counts.rds")) #Counts of expressed genes in samples of interest for the given tissue
  
  # Create DGEList object
  dge <- DGEList(counts)
  
  # Calculate normalization factors (does not do the normalization, only computes the factors)
  dge <- calcNormFactors(dge)
  
  # Voom
  v <- voom(dge, design = NULL, normalize="quantile", save.plot=F, plot = F) # samples are treated as replicates
  
  if(tissue %in% sexual_tissues){
    individual_traits <- c("Age", "Ancestry", "BMI", "Smoking")
  } else{
    individual_traits <- c("Age", "Ancestry", "Sex", "BMI", "Smoking")
  }
  covariates <- names(metadata)[!names(metadata) %in% c("Donor", "Sample", individual_traits)]
  
  # Expression ~ covariates 
  fml_args_mod <- paste(c(covariates), collapse = " + ")
  mod <- model.matrix( as.formula(paste(" ~  ", paste(fml_args_mod,collapse = " "))), data =  metadata)
  
  # Limma fit
  fit <- lmFit(v, mod)
  
  # Calculate expression residuals
  exprs_residuals <- resid(fit, v)
  
  saveRDS(exprs_residuals, paste0("tissues/", tissue, "/expression_residuals.rds"))
  
  #Including demographic traits for other purposes in other scripts:
  covariates <- c(covariates, individual_traits) 
  fml_args_mod <- paste(c(covariates), collapse = " + ")
  mod_2 <- model.matrix( as.formula(paste(" ~  ", paste(fml_args_mod,collapse = " "))), data =  metadata)
  fit_2 <- lmFit(v, mod_2)
  exprs_residuals_2 <- resid(fit_2, v)
  saveRDS(exprs_residuals_2, paste0("tissues/", tissue, "/expression_residuals_demographic_traits.rds"))
  
  # Differentially expressed genes
  traits <- names(dea_res) 
  de_genes <- unique(unlist(lapply(traits, function(trait) rownames(dea_res[[trait]][dea_res[[trait]][,"adj.P.Val"] < 0.05,]) ) ))

  # Getting the residuals of only the differentially expressed genes with at least one trait
  DE_residuals <- exprs_residuals[de_genes,]
  
  # hier.part
  print("Hier Part")
  hier.part <- lapply(rownames(DE_residuals), function(gene)
    hier.part.mod(y=DE_residuals[gene,], x=metadata[,individual_traits], fam = "gaussian", gof = "Rsqu"))
  names(hier.part) <- rownames(DE_residuals)
  
  # Parse information
  rsq <- sapply(rownames(DE_residuals), function(gene) 
    sum(hier.part[[gene]]$IJ[,1])) 
  names(rsq) <- rownames(DE_residuals)
  
  rel_perc <- do.call(rbind.data.frame,
                      lapply(names(hier.part), function(gene) 
                        as.numeric(unlist(hier.part[[gene]]$I.perc))))
  rownames(rel_perc) <- rownames(DE_residuals)
  colnames(rel_perc) <- individual_traits
  
  abs_perc <-  do.call(rbind.data.frame,
                       lapply(names(hier.part), function(gene) 
                         hier.part[[gene]]$IJ[,1]*100)
  )
  rownames(abs_perc) <- rownames(DE_residuals)
  colnames(abs_perc) <- individual_traits
  
  hier_data <- cbind.data.frame(rsq,rel_perc,abs_perc)
  colnames(hier_data) <- c("R2",paste0(individual_traits,"_rel"),paste0(individual_traits,"_abs"))
  
  saveRDS(hier_data, paste0("tissues/", tissue, "/hier.part.rds"))
  
  print("saved")
}


# -------------- #
print(Sys.time())
#-------------- #