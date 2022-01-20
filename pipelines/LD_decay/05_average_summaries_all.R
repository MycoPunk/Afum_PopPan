library(dplyr)
library(stringr)
library(ggplot2)
library(purrr)
library(data.table)

setwd("~/bigdata/pop_genomics/LD_decay_4/LD_out")

##Clade1
#list files to read in
dfr_1 <-list.files(path = ".",
                   pattern = glob2rx("Clade_1_samp*ea*"), 
                   full.names = T) %>% 
  map_df(~fread(.))
colnames(dfr_1) <- c("dist","rsq", "rsq.count")

#get mean distances R2 values for each distance across the replicates
rsq_means_Clade1<- aggregate(rsq~dist, data=dfr_1, FUN=function(x) c(mean=mean(x), count=length(x)))

#clean up memory
rm(dfr_1)

#write to file
write.table(rsq_means_Clade1, "rsq_means_n12_Clade_1.tab", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##Clade2
#list files to read in
dfr_2 <-list.files(path = ".",
                   pattern = glob2rx("Clade_2_samp*ea*"), 
                   full.names = T) %>% 
  map_df(~fread(.))
colnames(dfr_2) <- c("dist","rsq", "rsq.count")

#get mean distances R2 values for each distance across the replicates
rsq_means_Clade2<- aggregate(rsq~dist, data=dfr_2, FUN=function(x) c(mean=mean(x), count=length(x)))

#clean up memory
rm(dfr_2)

#write to file
write.table(rsq_means_Clade2, "rsq_means_n12_Clade_2.tab", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


##Clade3
#list files to read in
dfr_3 <-list.files(path = ".",
                   pattern = glob2rx("Clade_3_samp*ea*"), 
                   full.names = T) %>% 
  map_df(~fread(.))
colnames(dfr_3) <- c("dist","rsq", "rsq.count")

#get mean distances R2 values for each distance across the replicates
rsq_means_Clade3<- aggregate(rsq~dist, data=dfr_3, FUN=function(x) c(mean=mean(x), count=length(x)))

#clean up memory
rm(dfr_3)

#write to file
write.table(rsq_means_Clade3, "rsq_means_n12_Clade_3.tab", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##All
#list files to read in
dfr_all <-list.files(path = ".",
                   pattern = glob2rx("n_12_samp*ea*"), 
                   full.names = T) %>% 
  map_df(~fread(.))
colnames(dfr_all) <- c("dist","rsq", "rsq.count")

#get mean distances R2 values for each distance across the replicates
rsq_means_CladeAll<- aggregate(rsq~dist, data=dfr_all, FUN=function(x) c(mean=mean(x), count=length(x)))

#clean up memory
rm(dfr_all)

#write to file
write.table(rsq_means_CladeAll, "rsq_means_n12_Clade_all.tab", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
