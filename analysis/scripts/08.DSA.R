#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez and Raquel Garcia-Perez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to run differential expression analysis with smoking
# @software version: R=4.2.2

Sys.time()
start_time <- Sys.time()

# Loading libraries
suppressMessages(library(caret))
suppressMessages(library(parallel))
suppressMessages(library(hier.part))
suppressMessages(library(sandwich))
suppressMessages(library(lmtest))
suppressMessages(library(multcomp))
suppressMessages(library(optparse))
options(warn=-1)

#Set path 
setwd(system("pwd", intern = T)) #If in linux
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
# setwd("..")

# Parsing
parser <- OptionParser()
parser <- add_option(parser, opt_str=c("-t", "--tissue"), type="character",
                     dest="tissue",
                     help="Tissue")
options=parse_args(parser)
tissue=options$tissue
# tissue <- "Lung"

# Functions ####
# source(paste0(first_dir, "Raquel/R_functions/DEA_and_DSA.R_functions.R"))

# To save:
# -- PSI & TPM of alternatively spliced events (ASE)
# -- hier.part of ASE
# -- PSI residuals of ASE
# -- Differential splicing analysis (DSA) results


# Number of CPU (cores) to parallelize mclapply
n_cores <- 36
# n_cores <- 2


# Data ####
# Gene annotation ----
# PCG and lincRNA genes with mathched biotype annotation in gencode v38
# PAR genes excluded
gene_annotation <- read.csv("data/public/gene_annotation.csv")

sex.biased.genes <- gene_annotation[gene_annotation$chr=="chrY" |
                                    gene_annotation$gene.name=="XIST", "ensembl.id"]

# Transcript annotation ----
transcript_annotation <- read.delim("data/public/gencode.v26.GRCh38.transcripts.bed", header = F)
colnames(transcript_annotation) <- c("chr","start", "end","strand", "feature","ensembl.id", "transcript.id","gene.name", "g.biotype", "transcript.name","t.biotype")
transcript_annotation <- transcript_annotation[transcript_annotation$ensembl.id %in% gene_annotation$ensembl.id,]

# Exon annotation ----
exon_annotation <- read.delim("data/public/gencode.v26.GRCh38.exons.bed", header = F) #This file is too heavy to share via github, but I can be downloaded from gencode
colnames(exon_annotation) <- c("chr","start", "end","strand", "feature","ensembl.id", "transcript.id","gene.name", "transcript.name","exon.id","exon.numer","g.biotype", "t.biotype")
exon_annotation <- exon_annotation[exon_annotation$ensembl.id %in% gene_annotation$ensembl.id,]

# Event annotation ----
events.info <- readRDS("SUPPA/gencode.v26.splicing_events_coordinates.rds")

# Metadata ----
metadata <- readRDS(paste0("tissues/", tissue, "/metadata.rds"))

if(tissue %in% c("Prostate", "Testis", "Vagina", "Uterus", "Ovary")){
  metadata <- metadata[,names(metadata)!="Sex"]
}
if(length(unique(metadata$HardyScale))<5){metadata$HardyScale <- droplevels(metadata$HardyScale)}
if(length(unique(metadata$Ancestry))<3){metadata$Ancestry <- droplevels(metadata$Ancestry)}

# Transcript TPM ---- Do we need this?
transcript.tpm <- read.delim(paste0("SUPPA/TranscripExpressionFiles/", tissue, ".transcript_TPM.txt"))
colnames(transcript.tpm) <- gsub("\\.", "-", colnames(transcript.tpm))

# Reading in PSI and TPM values for alternative splicing events annotated in PC and lincRNA genes ----
psi <- readRDS(paste0("tissues/", tissue, "/psi_splicing_events.rds"))

tpm <- readRDS(paste0("tissues/", tissue, "/tpm_splicing_events.rds"))

# Expressed genes ----
exprs_genes <- readRDS(paste0("tissues/",tissue,"/expressed_genes.rds"))


# Create functions to be used:
is.event.exprs <- function(event.id, tpm.threshold = 0.5){
  if(which(rownames(psi) == event.id) %% 1000 == 0){ system(paste("echo 'Processed: ", which(rownames(psi) == event.id)," out of ", nrow(psi), "events'"))}
  # The two most abundant isforms in numeratos and denominator[!numerator] median TPM >= 1
  # Isoforms that include the event
  isoforms.in <- unlist(strsplit(events.info[event.id, "isoforms.spliced_in"],split = ","))
  # Isoforms that excluded the event
  isoforms.out <- unlist(strsplit(events.info[event.id, "isoforms.spliced_out"],split = ","))

  # Isoforms exprs value for each tissue sample
  isoforms.tpm <- lapply(metadata$Sample, function(sample)
    sapply(c(isoforms.in, isoforms.out), function(isoform)
      transcript.tpm[isoform,sample]
    ))
  names(isoforms.tpm) <- metadata$Sample

  # At least 20% of the samples express the most abundant isoform above threshold TPM ----
  # Most abundant isoform expression (median per isoform across samples)
  iso.in.most_abundant <- names(which.max(sapply(isoforms.in, function(isoform)
    median(sapply(metadata$Sample, function(sample)
      isoforms.tpm[[sample]][isoform]
    ))
  )))
  iso.out.most_abundant <- names(which.max(sapply(isoforms.out, function(isoform)
    median(sapply(metadata$Sample, function(sample)
      isoforms.tpm[[sample]][isoform]
    ))
  )))

  # 20% of the samples express the most abundant above "threshold" TPM
  condition1 <- sum(sapply(metadata$Sample, function(sample)
    isoforms.tpm[[sample]][iso.in.most_abundant] >= tpm.threshold
  )) >= round(0.2*nrow(metadata))
  condition2 <- sum(sapply(metadata$Sample, function(sample)
    isoforms.tpm[[sample]][iso.out.most_abundant] >= tpm.threshold
  )) >= round(0.2*nrow(metadata))

  return(condition1 & condition2)

}

mountBeta <- function(glmBatch,batchReg,mydata) {
  # mount beta coefficients
  
  nbeta <- 1
  ncoef <- 2
  beta  <- numeric(1)
  
  for (i in batchReg) {
    if (!is.factor(mydata[,i])) {
      beta[nbeta] <- glmBatch$coefficients[ncoef]
      nbeta <- nbeta + 1
      ncoef <- ncoef + 1
    }
    if (is.factor(mydata[,i])) {
      nlev        <- nlevels(mydata[,i])
      b1          <- -sum(glmBatch$coefficients[ncoef:(ncoef+nlev-2)])/nlev
      beta[nbeta] <- b1
      beta[(nbeta+1):(nbeta+nlev-1)] <- b1+glmBatch$coefficients[ncoef:(ncoef+nlev-2)] 
      nbeta <- nbeta + nlev
      ncoef <- ncoef + (nlev-1)
    }
  }
  return(beta)
}

computeResiduals <- function(mydata, batchReg, traitReg) {
  
  regressors <- c(batchReg,traitReg)
  myformeq <- paste("y", paste(colnames(mydata)[regressors], collapse=" + "), sep=" ~ ")
  glmFull  <- glm(myformeq, data = mydata,family = quasibinomial('logit'),control = list(maxit = 100) )
  
  ylogit    <- log(mydata$y/(1-mydata$y))
  myformeqB <- paste("~ ", paste(colnames(mydata)[batchReg], collapse=" + ")," - 1")
  if (length(batchReg) > 1) index   <- sapply(mydata[,batchReg],is.factor)
  if (length(batchReg) == 1) index  <- is.factor(mydata[,batchReg])
  fact      <- batchReg[index]
  if (length(fact) > 1) AUZ  <- model.matrix(as.formula(myformeqB), data=mydata[,regressors],
                                             contrasts.arg = lapply(mydata[,fact], contrasts, contrasts=FALSE))
  if (length(fact) <2 ) AUZ  <- model.matrix(as.formula(myformeqB), data=mydata[,regressors])
  
  #-- from model with only batch effects, taking out these batch effects
  
  myformeq <- as.formula(paste("y", paste(colnames(mydata)[batchReg], collapse=" + "), sep=" ~ "))
  glmBatch <- glm(myformeq , data = mydata,family = quasibinomial('logit'),control = list(maxit = 100) )
  
  # only glm with batch effects
  beta <- mountBeta(glmBatch,batchReg,mydata) 
  predBatch   <- AUZ %*% matrix(beta)
  
  ylogitRes1  <- ylogit - predBatch
  yRes1       <- 1/(1+exp(-ylogitRes1))
  myformeq    <- as.formula(paste("yRes1", paste(colnames(mydata)[traitReg], collapse=" + "), sep=" ~ "))
  glmRes1     <- glm(myformeq, data = mydata,family = quasibinomial('logit'),control = list(maxit = 100) )
  
  # -- from full model, taking out the batch effects
  SStot   <- sum((mydata$y-mean(mydata$y))^2)
  SStotR1 <- sum((yRes1-mean(yRes1))^2)
  R2Full  <- round((SStot - sum((glmFull$fitted.values- mydata$y)^2))/SStot,digits=5)
  R2Batch <- round((SStot - sum((glmBatch$fitted.values- mydata$y)^2))/SStot,digits=5)
  R2Res1  <- round((SStotR1 - sum((glmRes1$fitted.values- yRes1)^2))/SStotR1,digits=5)

  tableCoefs <- data.frame(Full=glmFull$coefficients,Batch=NA,Res1=NA)
  index   <- match(names(glmBatch$coefficients),names(glmFull$coefficients))
  tableCoefs[index,2] <- glmBatch$coefficients
  
  index   <- match(names(glmRes1$coefficients),names(glmFull$coefficients))
  tableCoefs[index,3] <- glmRes1$coefficients
  
  difGlobalMod1 <- sqrt(sum((glmFull$fitted.values-glmRes1$fitted.values)^2))/length(glmFull$fitted.values)

  retObj <- list(c(R2Full=R2Full,R2Batch=R2Batch,R2Res1=R2Res1),
                 c(difGlobalMod1=difGlobalMod1),
                 tableCoefs,
                 cleanedData=data.frame(cleanMod1=yRes1))
  
  return(retObj)  
}

get_residuals <- function(event_id, mdata){
  # Model one event at a time
  y <- pmin(pmax(as.numeric(psi[event_id,]),0),1) 
  glm_data <- cbind.data.frame(mdata, y)
  
  obj <- computeResiduals(glm_data, which(colnames(mdata) %in% covariates), which(colnames(mdata) %in% individual_traits))
  return(obj$cleanedData$cleanMod1)
}

# hier.part functions ####
# ---- methods derived from hier.part
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
    if (gof == "RMSPE") new.line <- current.model.mod(y, current.comb, xcan,family = family, gof = "RMSPE",...)
    if (gof == "Rsqu")  new.line <- current.model.mod(y, current.comb, xcan,family = family, SStot=SStot,gof = "Rsqu",...)
    gfs <- c(gfs, new.line)
  }
  gfs
}

hier.part.mod <- function(y,xcan,family='gaussian',gof = "Rsqu", link = "",...) {
  pcan <- dim(xcan)[2]
  gfs <- all.regs.mod(y, xcan, family = family, gof = gof, link = link, ...)
  hp <- partition(gfs, pcan, var.names = names(data.frame(xcan)))
  
  params <- list(full.model = paste("y ~", paste(names(xcan),collapse = " + ")), 
                 family = family, link = link, gof = gof)
  if(sum(hp$IJ$I<0)>0){
    NA
  }else{
    list(gfs = gfs, IJ = hp$IJ, I.perc = hp$I.perc, params = params)
  }
}

# 1. Filter AS events ----
print("# ---- Selecting alternatively spliced events ---- #")
print(paste0("Number of ASE in PC and lincRNA genes: ", nrow(psi)))

# * Track number of alternative splicing events ----
no.ase <- vector()
no.ase <- nrow(psi)

# 1.1 Subset events in PCG and lincRNA expressed in tissue ----

# Keep events in expressed PC and lincRNA --
psi$ensembl_id <- sapply(rownames(psi), function(gene) unlist(strsplit(gene,split = ";"))[[1]])
tpm$ensembl_id <- sapply(rownames(tpm), function(gene) unlist(strsplit(gene,split = ";"))[[1]])

psi <- psi[psi$ensembl_id %in% exprs_genes,]
tpm <- tpm[tpm$ensembl_id %in% exprs_genes,]
psi <- psi[!psi$ensembl_id %in% sex.biased.genes,]
tpm <- tpm[!tpm$ensembl_id %in% sex.biased.genes,]

psi <- psi[,-ncol(psi)]
tpm <- tpm[,-ncol(tpm)]


#Subset for testing
# psi <- psi[1:50,]

# * Track number of alternative splicing events ----
no.ase <- c(no.ase, nrow(psi))

# 1.2 Retain events with no samples having an NA ----

# Retain events with no samples having an NA
number_na <- rowSums(is.na(psi),na.rm=F) # vector with number of samples with NA per event
noNA_events.psi <- names(number_na)[number_na==0] # events with 0 NAs
psi <- psi[noNA_events.psi, ]

# * Track number of alternative splicing events ----
no.ase <- c(no.ase, nrow(psi))

# 1.3 Calculate event residuals ####
print("# ---- Calculating residuals ---- #")

# 1.3.1 Variables ----  I changed it becaus if I have more than 1 disease, I want them to be individual traits (e.g., Smoking and Pneumonia)
covariates <- names(metadata)[!names(metadata) %in% c("Donor", "Sample", "Smoking", "Age", "Ancestry", "Sex", "BMI")]
individual_traits <- names(metadata)[!names(metadata) %in% c(covariates, "Donor", "Sample")]


# 1.3.2 Compute residuals ----
# -------------- #
print(Sys.time())
# -------------- #
# Residuals    20 min
fr <- mclapply(rownames(psi), function(event_id) get_residuals(event_id, metadata), mc.cores = n_cores )
names(fr) <- rownames(psi)
# -------------- #
print(Sys.time())
# -------------- #

# 1.3.3 Create dataframe ----
psi_residuals <- do.call(rbind.data.frame,
                         fr)
colnames(psi_residuals) <- colnames(psi)
rownames(psi_residuals) <- rownames(psi)
psi_residuals <- round(psi_residuals, 2)

# 1.4 Exclude events with low complexity ----

# exclude events with fewer than max(10, 0.1n) unique values, where n is the sample size
psi.complexity <- apply(psi, 1, function(x) length(unique(x)))
# exclude events with fewer than max(10, 0.1n) unique values, where n is the sample size
psi_residuals.complexity <- apply(psi_residuals, 1, function(x) length(unique(x)))

threshold_complexity <- 15

psi <- psi[intersect(names(psi.complexity[psi.complexity >= threshold_complexity]),
                              names(psi_residuals.complexity[psi_residuals.complexity >= threshold_complexity])
),]
psi_residuals <- psi_residuals[intersect(names(psi.complexity[psi.complexity >= threshold_complexity]),
                                         names(psi_residuals.complexity[psi_residuals.complexity >= threshold_complexity])
                                         ),]


# * Track number of alternative splicing events ----
no.ase <- c(no.ase, nrow(psi))

# 1.5 Exclude events with insufficient variability ----

event_freq <- apply(psi, 1, function(x) sort(table(x),decreasing = T)[1]/sort(table(x),decreasing = T)[2] < 80/20)
event_freq_residuals <- apply(psi_residuals, 1, function(x) sort(table(x),decreasing = T)[1]/sort(table(x),decreasing = T)[2] < 80/20)
psi <- psi[event_freq & event_freq_residuals,]
psi_residuals <- psi_residuals[event_freq & event_freq_residuals,]

# * Track number of alternative splicing events ----
no.ase <- c(no.ase, nrow(psi))

save.image(file = paste0("err/", tissue, "_filter_data_0.RData"))
# load(paste0("err/", tissue,"_filter_data.RData"))

# 1.6 Exclude events not sufficiently expressed ----
print("# ---- Exclude events not sufficiently expressed ---- #")
# -------------- #
print(Sys.time())
# -------------- #
events_exprs <- unlist(mclapply(rownames(psi), function(i) is.event.exprs(i, 0.5),  mc.cores = n_cores  ))
# -------------- #
print(Sys.time())
# -------------- #
save.image(file = paste0("err/", tissue, "_filter_data.RData"))
# load(paste0("err/", tissue,"_filter_data.RData"))

# Subset ASE sufficiently exprs ----
psi <- psi[events_exprs,]
psi_residuals <- psi_residuals[rownames(psi),]
tpm <- tpm[rownames(psi),]

# * Track number of alternative splicing events ----
no.ase <- c(no.ase, nrow(psi))

# 2. hier.part ----
print("# ---- Running hier.part ---- #")

# 2.1 Calculate explained variance ----
# -------------- #
print(Sys.time())
# -------------- #
hier.part.results <- mclapply(rownames(psi_residuals), function(event)
  hier.part.mod(y=as.numeric(psi_residuals[event,]), x=metadata[,individual_traits],
                fam = "quasibinomial", link = "logit", gof = "Rsqu",control = list(maxit = 100)), mc.cores = n_cores)
names(hier.part.results) <- rownames(psi_residuals)
# -------------- #
print(Sys.time())
# -------------- #

# 2.2 Exclude events with negative estimates ----
print(paste0("ASE events with unestimable contributions: ", sum(is.na(hier.part.results))))
if(length(which(is.na(hier.part.results)))>0){
  psi <- psi[-which(is.na(hier.part.results)),]
  psi_residuals <- psi_residuals[-which(is.na(hier.part.results)),]
  hier.part.results <- hier.part.results[-which(is.na(hier.part.results))]
}

# * Track number of alternative splicing events ----
no.ase <- c(no.ase, nrow(psi))

# 2.3 Parse results ----
rsq <- sapply(rownames(psi_residuals), function(event)
  sum(hier.part.results[[event]]$IJ[,1]))
names(rsq) <- rownames(psi_residuals)
rel_perc <- do.call(rbind.data.frame,
                    lapply(names(hier.part.results), function(event)
                      as.numeric(unlist(hier.part.results[[event]]$I.perc))))
rownames(rel_perc) <- rownames(psi_residuals)
colnames(rel_perc) <- individual_traits
abs_perc <-  do.call(rbind.data.frame,
                     lapply(names(hier.part.results), function(event)
                       hier.part.results[[event]]$IJ[,1])
)
rownames(abs_perc) <- rownames(psi_residuals)
colnames(abs_perc) <- individual_traits
hier_data <- cbind.data.frame(rsq,rel_perc,abs_perc)
colnames(hier_data) <- c("R2",paste0(individual_traits,"_rel"),paste0(individual_traits,"_abs"))

save.image(file = paste0("err/", tissue, "_filter_data_2.RData"))
# load(file = paste0("err/", tissue, "_filter_data_2.RData"))
# 3. Differential splicing analysis ----
print("# ---- Running differential splicing analysis ---- #")

get_tables <- function(event_id, glm_model, mult, mult_ancestry, metadata){
  df <-  as.data.frame(coef(summary(glm_model)))
  df <- df[contrast_names,]
  
  df$Trait <- contrast_names
  df$Ensembl_id <- unlist(strsplit(event_id, split = ";"))[[1]]
  df$Event_id <- event_id
  df$Event <- unlist(strsplit(event_id, split = ";"))[[2]]
  df$Type <- unlist(strsplit(df$Event, split = ":"))[[1]]
  df$`t value`<- NULL
  df$`Std. Error` = NULL
  df$glm.p.value = df$`Pr(>|t|)` 
  # ----- robust sandwich error #
  cft <- coeftest(glm_model, vcov.=vcovHC(glm_model, type="HC0"))
  df$Std.Error <- cft[contrast_names, 2]
  df$Z.Value <- cft[contrast_names, 3]
  df$P.Value <- cft[contrast_names, 4]
  df$Iter <- glm_model$iter

  mult_cft <- coeftest(mult) #To add the last comparison needed Smoking1-2. P-values of Smoking 1 and Smoking2 match with the ones we get here. 
  df <- rbind(df, c(mult$test$coefficients[3], mult$test$pvalues[3], "Smoking1-2", 
                    unlist(strsplit(event_id, split = ";"))[[1]], event_id, 
                    unlist(strsplit(event_id, split = ";"))[[2]], unlist(strsplit(df$Event, split = ":"))[[1]],
                    mult$test$pvalues[3], mult_cft[3,2:4], glm_model$iter))
  rownames(df)[nrow(df)] <- "Smoking1-2"
  if("AMR" %in% levels(metadata$Ancestry)){
    mult_cft_ancestry <- coeftest(mult_ancestry) #To add the last comparison needed AncestryAMR-EUR. 
    df <- rbind(df, c(mult_ancestry$test$coefficients[3], mult_ancestry$test$pvalues[3], "AncestryAMR-EUR", 
                      unlist(strsplit(event_id, split = ";"))[[1]], event_id, 
                      unlist(strsplit(event_id, split = ";"))[[2]], unlist(strsplit(df$Event, split = ":"))[[1]],
                      mult_ancestry$test$pvalues[3], mult_cft_ancestry[3,2:4], glm_model$iter))
    rownames(df)[nrow(df)] <- "AncestryAMR-EUR"
  }
  
  # -------------------------------------- #
  # If warning, event might be modelled properly and function might return an NA (NaN)
  # Do this to keep track of warnigns. 
  # Note how events with warnings are associated with extremely large betas
  cft <- tryCatch({
    coeftest(glm_model, vcov.=vcovHC(glm_model, type="HC0"))
  }, warning = function(w) {
    return(NA)
    # "Warning message:
    #   In sqrt(diag(se)) : NaNs produced"
  }  
  )
  # Subset data used in downstream analysis
  if(!is.na(cft[1])){
    df$coeftest_Warning <- rep("0",nrow(df)) # 0 are events with no warnings
  }else{
    df$coeftest_Warning <- rep("1", nrow(df)) # 1 are events with  warnings
  }
  df <- df[,c("Event_id",
              "Trait",
              "Ensembl_id",
              "Type",  
              "Event",
              "Estimate",
              "glm.p.value",
              "Iter",
              "Std.Error",
              "Z.Value", 
              "P.Value",
              "coeftest_Warning")]
  return(df)
}

# Function to fit a glm model per event
model_psi <- function(event_id, mdata){
  # Model one event at a time
  psi_values <- pmin(pmax(as.numeric(psi[event_id,]),0),1) # forzar los límites de PSI a los valores teóricos de mínimo y máximo
  glm_data <- cbind.data.frame(mdata, psi_values)
  
  # Model formula: psi ~ covariates + traits
  mod_formula <- as.formula(paste("psi_values ~ ", paste(c(covariates, individual_traits), collapse = " + "), collapse = " ") )
  myglm <- glm(mod_formula, 
               data = glm_data, 
               family = quasibinomial('logit'),
               control = list(maxit = 100))

  #Compare with ex-smokers:
  mult <- summary(glht(myglm, mcp(Smoking="Tukey"), vcov. = vcovHC(myglm, type="HC0")))
  #Compare different ancestries between them
  if("AMR" %in% levels(metadata$Ancestry)){ #Some tissues do not have enough AMR
    mult_ancestry <- summary(glht(myglm, mcp(Ancestry="Tukey"), vcov. = vcovHC(myglm, type="HC0")))
  }
  # If there is a glm warning cause algorithm did not converge 
  my_warning <- tryCatch({
    glm(mod_formula, 
        data = glm_data, 
        family = quasibinomial('logit'),
        control = list( maxit = 100) )
  }, warning = function(w) {
    return(NA)
    # Warning message:
    #  glm.fit: algorithm did not converge 
  }  
  )
  
  # Create table
  event_results <- get_tables(event_id, myglm, mult, mult_ancestry, mdata)
  
  if(!is.na(my_warning[1])){
    event_results$glm_Warning <- rep("0", nrow(event_results)) # 0 are events with no warnings
  }else{
    event_results$glm_Warning <- rep("1", nrow(event_results)) # 1 are events with  warnings
  }
  
  return(list('Event_id' = event_id,
              'glm' = myglm,
              'res' = event_results))  # table with estimates and P.Values per event
}

contrast_names <- c("Age", "Sex2", "BMI", "Smoking1", "Smoking2", "AncestryAMR", "AncestryEUR")
if(!"Sex" %in% names(metadata)){ contrast_names <- contrast_names[!contrast_names %in% "Sex2"]}
if(!"AMR" %in% levels(metadata$Ancestry)){ contrast_names <- contrast_names[!contrast_names %in% "AncestryAMR"]}

# 3.1 Run PSI models ----
# -------------- #
print(Sys.time())
# -------------- #
# One model per event
fr <- mclapply(rownames(psi), function(event_id) model_psi(event_id, metadata), mc.cores = n_cores )
names(fr) <- rownames(psi)
# -------------- #
print(Sys.time())
# -------------- #

# 3.2 Parse results table ----
contrast_names <- c(contrast_names, "Smoking1-2", "AncestryAMR-EUR")
if(!"AMR" %in% levels(metadata$Ancestry)){ contrast_names <- contrast_names[!contrast_names %in% "AncestryAMR-EUR"]}

results <- do.call(rbind.data.frame,
                   lapply(rownames(psi), function(event_id) fr[[event_id]][['res']]))
dsa_res <- lapply(contrast_names, function(trait)
  results[results$Trait==trait,])
names(dsa_res) <- contrast_names

for(trait in contrast_names){
  # Set event ID as as rownames
  rownames(dsa_res[[trait]]) <- dsa_res[[trait]]$Event_id
  dsa_res[[trait]] <- dsa_res[[trait]][,-1]
}

# 3.3 Exclude events with warnings ----
# -- Beware of warnings -- #
# The coeftest function raises warning at particular events when the variance-covariance matrix cannot be computed
# In extreme instances it returns NA, for example, instances of events completely stratified that we excluded before-hand
# Yet, some events might remain that cannot be modelled
# If coef.test raises a warning we consider those events are not properly modelled and assign a P.Value of NA
# Multiple testing pvalue correction ->
# if in the p-values vector (not corrected) there is an NA, the p.adjust() does not consider it in the n (number of observations)
# Remove events with warnings in coeftest function
# These warnings appear in tissues with low sample size for events with low variance
# We consider these events cannot be modelled
print(paste0("ASE that raised glm warning: ", sum(dsa_res$Age$glm_Warning == "1") ))
print(paste0("ASE that raised coef.test warning: ", sum(dsa_res$Age$coeftest_Warning == "1") ))

# * Track number of alternative splicing events ----
no.ase <- c(no.ase, nrow(psi) - sum(dsa_res$Age$glm_Warning == "1"))
no.ase <- c(no.ase, nrow(psi) - sum(dsa_res$Age$glm_Warning == "1") - sum(dsa_res$Age$coeftest_Warning == "1"))
no.ase <- c(no.ase, nrow(dsa_res[["Age"]]))

if(sum(dsa_res$Age$glm_Warning == "1") > 0 ){
  for(trait in contrast_names){
    dsa_res[[trait]]$P.Value[which(dsa_res[[trait]]$glm_Warning=="1")] <- NA
  }
}
if(sum(dsa_res$Age$coeftest_Warning == "1") > 0 ){
  for(trait in contrast_names){
    dsa_res[[trait]]$P.Value[which(dsa_res[[trait]]$coeftest_Warning=="1")] <- NA
  }
}
for(trait in contrast_names){
  dsa_res[[trait]] <- dsa_res[[trait]][!is.na(dsa_res[[trait]]$P.Value),]
}

# 3.4 FDR correction ---
for(trait in contrast_names){
  dsa_res[[trait]]$adj.P.Val <- p.adjust(dsa_res[[trait]]$P.Value, method = "BH")
}
print(paste0("ASE modelled: ", nrow(dsa_res$Age)))

# 3.5 Save results ----
trait <- "Age"
psi <- psi[rownames(dsa_res[[trait]]),]
psi_residuals <- psi_residuals[rownames(dsa_res[[trait]]),]
tpm <- tpm[rownames(dsa_res[[trait]]),]
hier_data <- hier_data[rownames(dsa_res[[trait]]),]

# PSI values
saveRDS(psi, paste0("tissues/", tissue, "/Alternatively_spliced_events.psi.rds"))
# TPM values
saveRDS(tpm, paste0("tissues/", tissue, "/Alternatively_spliced_events.tpm.rds"))
# Save residuals
saveRDS(psi_residuals, paste0("tissues/", tissue, "/Alternatively_spliced_events.psi_residuals.rds"))
# Save hier.part
saveRDS(hier_data, paste0("tissues/", tissue, "/Alternatively_spliced_events.hier_part.rds"))
# Save filtering
saveRDS(no.ase, paste0("tissues/", tissue, "/Alternatively_spliced_events.Filtering.rds"))
#Save DSA results
# saveRDS(dsa_res, paste0("tissues/", tissue, "/fractional_regression_results.rds")	)

# 4. DeltaPSI and average TPM ----
print("# ---- Calculating deltaPSI ---- #")

# 5.2 Compute average TPM expression and variance for each event
avrg_TPM <- apply(tpm, 1, function(x) mean(log2(x+1)))
var_TPM <- apply(tpm, 1, function(x) var(x))

AFR_samples <- metadata[metadata$Ancestry=="AFR", "Sample"]
EUR_samples <- metadata[metadata$Ancestry=="EUR", "Sample"]
AMR_samples <- metadata[metadata$Ancestry=="AMR", "Sample"]
Female_samples <- metadata[metadata$Sex=="2", "Sample"]
Male_samples <- metadata[metadata$Sex=="1", "Sample"]
obese_samples <- metadata[metadata$BMI>=30, "Sample"]
normal_samples <- metadata[metadata$BMI<25, "Sample"]

# 5.2 Compute deltaPSI ----
get.deltaPSI <- function(trait){
  print(trait)
  # -------------- #
  print(Sys.time())
  # -------------- #)

  if(trait == "Sex" | trait == "Sex2"){
    if(tissue %in% c("Ovary","Uterus","Vagina","Testis","Prostate")){
      delta.psi <- NA
    }else{
      delta.psi <- sapply(rownames(dsa_res[[trait]]), function(event)
        mean(as.numeric(psi[event,Female_samples]),na.rm=T) - mean(as.numeric(psi[event,Male_samples]),na.rm=T) )
    }
  }else if(trait == "Age"){
    younger_samples <- metadata[metadata$Age<45,"Sample"] # [20-45)
    older_samples <- metadata[metadata$Age>=45,"Sample"] # [45-70]
    delta.psi <- sapply(rownames(dsa_res[[trait]]), function(event)
      mean(as.numeric(psi[event,older_samples]),na.rm=T) - mean(as.numeric(psi[event,younger_samples]),na.rm=T))
  }else if(trait == "BMI"){
    delta.psi <-  sapply(rownames(dsa_res[[trait]]), function(event)
      mean(as.numeric(psi[event,obese_samples]),na.rm=T) - mean(as.numeric(psi[event,normal_samples]),na.rm=T))
  } else if(trait == "Smoking1"){
    healthy_samples <- metadata[metadata[["Smoking"]]=="0", "Sample"]
    disease_samples <- metadata[metadata[["Smoking"]]=="1", "Sample"]
    delta.psi <- sapply(rownames(dsa_res[[trait]]), function(event)
      mean(as.numeric(psi[event,disease_samples]),na.rm=T) - mean(as.numeric(psi[event,healthy_samples]),na.rm=T))
  }else if(trait == "Smoking2"){
    healthy_samples <- metadata[metadata[["Smoking"]]=="0", "Sample"]
    disease_samples <- metadata[metadata[["Smoking"]]=="2", "Sample"]
    delta.psi <- sapply(rownames(dsa_res[[trait]]), function(event)
      mean(as.numeric(psi[event,disease_samples]),na.rm=T) - mean(as.numeric(psi[event,healthy_samples]),na.rm=T))
  }else if(trait == "Smoking1-2"){
    healthy_samples <- metadata[metadata[["Smoking"]]=="1", "Sample"]
    disease_samples <- metadata[metadata[["Smoking"]]=="2", "Sample"]
    delta.psi <- sapply(rownames(dsa_res[[trait]]), function(event)
      mean(as.numeric(psi[event,disease_samples]),na.rm=T) - mean(as.numeric(psi[event,healthy_samples]),na.rm=T))
  }else if(trait == "AncestryAMR"){
    delta.psi <- sapply(rownames(dsa_res[[trait]]), function(event)
      mean(as.numeric(psi[event,AMR_samples]),na.rm=T) - mean(as.numeric(psi[event,AFR_samples]),na.rm=T))
  }else if(trait == "AncestryEUR"){
    delta.psi <- sapply(rownames(dsa_res[[trait]]), function(event)
      mean(as.numeric(psi[event,EUR_samples]),na.rm=T) - mean(as.numeric(psi[event,AFR_samples]),na.rm=T))
  }else if(trait == "AncestryAMR-EUR"){
    delta.psi <- sapply(rownames(dsa_res[[trait]]), function(event)
      mean(as.numeric(psi[event,EUR_samples]),na.rm=T) - mean(as.numeric(psi[event,AMR_samples]),na.rm=T))
  }else { #This will never be evaluated, but it is ready for other diseases 
    healthy_samples <- metadata[metadata[[trait]]=="0", "Sample"]
    disease_samples <- metadata[metadata[[trait]]=="1", "Sample"]
      delta.psi <- sapply(rownames(dsa_res[[trait]]), function(event)
        mean(as.numeric(psi[event,disease_samples]),na.rm=T) - mean(as.numeric(psi[event,healthy_samples]),na.rm=T))
  }
  # -------------- #
  print(Sys.time())
  # -------------- #
  return(delta.psi)
}

# Add deltaPSI, AvgExprs, ExprsVar and order data.frame
for(trait in contrast_names){
  if(trait == "Sex" | trait == "Sex2"){
    if(tissue %in% c("Ovary","Uterus","Vagina","Testis","Prostate")){
      dsa_res[[trait]] <- NA
    }else{
      # Add deltaPSI
      dsa_res[[trait]]$deltaPSI <- get.deltaPSI(trait)
      # Add average expression TPM
      dsa_res[[trait]]$AvgExprs <- sapply(rownames(dsa_res[[trait]]), function(event) avrg_TPM[event])
      dsa_res[[trait]]$ExprsVar <- sapply(rownames(dsa_res[[trait]]), function(event) var_TPM[event])
      # Order by FDR
      dsa_res[[trait]] <- dsa_res[[trait]][order(dsa_res[[trait]]$adj.P.Val),]
    }
  }else{
    # Add deltaPSI
    dsa_res[[trait]]$deltaPSI <- get.deltaPSI(trait)
    # Add average expression TPM and variance
    dsa_res[[trait]]$AvgExprs <- sapply(rownames(dsa_res[[trait]]), function(event) avrg_TPM[event])
    dsa_res[[trait]]$ExprsVar <- sapply(rownames(dsa_res[[trait]]), function(event) var_TPM[event])
    # Order by FDR
    dsa_res[[trait]] <- dsa_res[[trait]][order(dsa_res[[trait]]$adj.P.Val),]
  }
}

print("Saving dsa results")
# 5.3 Save data ----
saveRDS(dsa_res, paste0("tissues/", tissue, "/fractional_regression_results.rds"))

# ---------------------- #
end_time <- Sys.time()
print("")
print("# ---- Elapsed time ---- #")
end_time - start_time
print("# ---------------------- #\n")

print("")

