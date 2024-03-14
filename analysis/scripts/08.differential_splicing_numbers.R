#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to analyze events consequences at the functional level part 2
# @software version: R=4.2.2

#Set path 
# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")


data <- readRDS(paste0("output/00.differential_splicing_analysis.smoking.rds"))

table <- rep("1", 12)

for(tissue in names(data)){ #Add in table the disease ones and then filter by PC-PC, get number and check in which the domain is changed
  if(nrow(data[[tissue]]$`SmokingSMOKER-NEVER`[data[[tissue]]$`SmokingSMOKER-NEVER`$adj.P.Val<0.05,])==0){next}
  table <- rbind(table, cbind(data[[tissue]]$`SmokingSMOKER-NEVER`[data[[tissue]]$`SmokingSMOKER-NEVER`$adj.P.Val<0.05,], tissue))
}
table <- table[-1,]

# testing <- rownames(table)
# testing <- gsub("\\+1","+", testing) #There is only one duplicated event
# testing <- gsub("\\-1","-", testing)

test <- table
test <- test[rownames(test)!="ENSG00000159618.15;RI:chr16:57565034:57565150-57566599:57566751:+1",] #The duplicated one
test <- test[test$biotype!="NA-NA",]

test <- test[test$biotype=="PC-PC",] #Only one shared
round(83/173, 2) #48 % of the DSE are Pc-PC
test <- test[!(is.na(test$spliced_in_domains) & is.na(test$spliced_out_domains)),] #If there is a change in known pfam
test[is.na(test)] <- "NA"
test <- test[!test$spliced_in_domains==test$spliced_out_domains,] #34
round(38/83, 2)

test <- table
test <- test[rownames(test)!="ENSG00000159618.15;RI:chr16:57565034:57565150-57566599:57566751:+1",] #The duplicated one
test <- test[test$biotype!="NA-NA",]
test <- test[test$biotype=="NC-PC" | test$biotype=="PC-NC",] #70
round(76/173, 2)

#Are most NC in changes PC-NC associated to smoking? Do a binomial test. Positive beta increase the inclusion
sum(test$beta>0 & test$biotype=="PC-NC") + sum(test$beta<0 & test$biotype=="NC-PC")
sum(test$beta<0 & test$biotype=="PC-NC") + sum(test$beta>0 & test$biotype=="NC-PC") #This is changes to NC due to smoking as PC-NC means more protein coding in smoking
round(49/76, 2)
binom.test(49, 76, 0.5)


#How many events of each types we consider? For methods:
events <- readRDS("SUPPA/gencode.v26.splicing_events_coordinates.rds")

dsa_res <- readRDS("output/00.differential_splicing_analysis.smoking.rds")
test <- lapply(dsa_res, function(sublist) sublist$"SmokingSMOKER-NEVER")
test_df <- data.frame(do.call(rbind, test))
test_df$events <- sub("^[^.]+\\.", "", rownames(test_df))
events <- unique(test_df$events)
length(events) #These are the events have been tested in at least one tissue, we considered 160,000 but only 57,002 passed the filters
event_types <- sub("^[^;]+;([^:]+).*", "\\1", events)
table(event_types)
