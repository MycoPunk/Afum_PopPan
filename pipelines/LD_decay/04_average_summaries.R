library(dplyr)
library(stringr)
library(ggplot2)
library(purrr)
library(data.table)

setwd("~/bigdata/pop_genomics/LD_decay/LD_out")

##Clade1
#list files to read in
dfr_1 <-list.files(path = ".",
                   pattern = "Clade_1_samp*", 
                   full.names = T) %>% 
  map_df(~fread(.))
colnames(dfr_1) <- c("dist","rsq")

#get mean distances R2 values for each distance across the replicates
rsq_means_Clade1<- aggregate(rsq~dist, data=dfr_1, FUN=function(x) c(mean=mean(x), count=length(x)))

#clean up memory
rm(dfr_1)

#write to file
write.table(rsq_means_Clade1, "rsq_means_Clade1.tab", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#clean up memory
rm(rsq_means_Clade1)


##Clade2
dfr_2 <-list.files(path = ".",
                   pattern = "Clade_2_samp*", 
                   full.names = T) %>% 
  map_df(~fread(.))
colnames(dfr_2) <- c("dist","rsq")

#get mean distances R2 values for each distance across the replicates
rsq_means_Clade2<- aggregate(rsq~dist, data=dfr_2, FUN=function(x) c(mean=mean(x), count=length(x)))

#clean up memory
rm(dfr_2)

#write to file
write.table(rsq_means_Clade2, "rsq_means_Clade2.tab", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#clean up memory
rm(rsq_means_Clade2)


##Clade3
dfr_3 <-list.files(path = ".",
                   pattern = "Clade_3_samp*", 
                   full.names = T) %>% 
  map_df(~fread(.))
colnames(dfr_3) <- c("dist","rsq")

#get mean distances R2 values for each distance across the replicates
rsq_means_Clade3<- aggregate(rsq~dist, data=dfr_3, FUN=function(x) c(mean=mean(x), count=length(x)))

#clean up memory
rm(dfr_3)

#write to file
write.table(rsq_means_Clade3, "rsq_means_Clade3.tab", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
