#A fum pan genome analysis: 
#validate parameters for PIRATE pangenome 
#here, AF293 should have approximately the same number of genes in gene families as there are gene models for Af293
#there are currently 9,840 protein coding gene models for Af293 (as of September 2021)


#set packages 
library(data.table)
library(tidyverse)
library(hrbrthemes)
library(forcats)
library(ggforce)
library(ggimage)
library(ggtree)
library(ape)
library(phytools)
library(ggplot2)
library(stringr)

#set wd
setwd("~/Desktop/Project_Afum_pangenome_2/")

#read in OF gene families
OF.gene_counts<-as.data.frame(fread("Orthogroups.GeneCount.tsv")) 

#clean up names
#fix names
names(OF.gene_counts)<- sapply(names(OF.gene_counts), gsub, pattern = "Aspergillus_fumigatus_", replacement = "" )
names(OF.gene_counts)<- sapply(names(OF.gene_counts), gsub, pattern = ".proteins", replacement = "" )

##note the gene counts file does not include singletons that were not clustered into any orthgroup- 
#-note some of what we call "singletons" in the paper do cluster into orthogroups, but they are singlesons because we only find them in one strain.
#so we need to add that back in the singletons tht didn't cluster at all in case any of these are in the AF293 resequenced strain

#read in the singletons
OF.unassigned<- as.data.frame(fread("Orthogroups_UnassignedGenes.tsv")) #this is where your singletons live
#fix names
names(OF.unassigned)<- sapply(names(OF.unassigned), gsub, pattern = "Aspergillus_fumigatus_", replacement = "" )
names(OF.unassigned)<- sapply(names(OF.unassigned), gsub, pattern = ".proteins", replacement = "" )

##change protein ID to counts (all are 1 because these are singletons)
#first change all empty cells to NA
OF.unassigned<- data.frame(apply(OF.unassigned, 2, function(x) gsub("^$|^ $", NA, x)))

#exclude OG col
cols_to_exclude<- "Orthogroup"
strains_only<- OF.unassigned[,!names(OF.unassigned) %in% cols_to_exclude]
dim(OF.unassigned)
dim(strains_only)

#replace all NAs with 0s
strains_only[is.na(strains_only)] <- 0
#fill in ones
strains_only_ones<- as.data.frame(replace(strains_only, strains_only!="0", 1))
#change to numeric
strains_only_ones_num <- mutate_all(strains_only_ones, function(x) as.numeric(as.character(x)))

#check 
test<- rowSums(strains_only_ones_num)
sum(test != 1) #works as expected

#add back Orthogroups col
singles<- cbind(Orthogroup = OF.unassigned$Orthogroup, strains_only_ones_num, Total = rowSums(strains_only_ones_num))

#fix X's introduced in names starting with numbers
names(singles)<- sapply(names(singles), gsub, pattern = "^X", replacement = "" )

#combine datasets
all_OFs_counts<- rbind(OF.gene_counts,singles)
dim(all_OFs_counts)

#isolate the resequenced strain
AF293_only<- data.frame(cbind(Orthogroup = all_OFs_counts$Orthogroup, AF293_ct = all_OFs_counts$AF293))

#how many genes in total distributed across all gene families?
sum(as.numeric(AF293_only$AF293_ct)) #10167

#how many gene families is this? 
AF293_only_not_zero<- AF293_only[as.numeric(AF293_only$AF293_ct != 0),]
dim(AF293_only_not_zero)
#difference of 
sum(as.numeric(AF293_only$AF293_ct == 0))
