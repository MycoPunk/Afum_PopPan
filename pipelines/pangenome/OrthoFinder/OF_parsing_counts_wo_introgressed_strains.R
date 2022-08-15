#based on code used to look at A fum pan genome analysis, using proteins clustered by Ortho Finder results
#this scrip generates count numbers for clade-specific gene families, but with the ingrogressed strains removed. 
#last updated: 15.Aug.2022

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
library(ggpubr)
library(perm)
library(phonTools)
library(vegan)
packageVersion("vegan") 

setwd("")
OF.gene_counts<-as.data.frame(fread("Orthogroups.GeneCount.tsv")) 
OF.gene_families<-as.data.frame(fread("Orthogroups.tsv")) 
OF.unassigned<- as.data.frame(fread("Orthogroups_UnassignedGenes.tsv")) #this is where your singletons live


#combine OF families with singletons
OF.gene_families_all<- rbind(OF.gene_families, OF.unassigned)


#fix names
names(OF.gene_families_all)<- sapply(names(OF.gene_families_all), gsub, pattern = "Aspergillus_fumigatus_", replacement = "" )
names(OF.gene_families_all)<- sapply(names(OF.gene_families_all), gsub, pattern = ".proteins", replacement = "" )

#list strains to remove
to_remove<- c("F7763", 
              "CF098",
              "AF293",
              "F18304",
              "CM7632",
              "10_01_02_27",
              "12_7505220",
              "RSF2S8")

#remove strains
OF.gene_families_all_no_intro<- OF.gene_families_all[,!names(OF.gene_families_all) %in% to_remove]


cols_to_exclude2<- "Orthogroup"
strains_only2<- OF.gene_families_all_no_intro[,!names(OF.gene_families_all_no_intro) %in% cols_to_exclude2]
strain_names2<- names(strains_only2)

#length(strain_names)
ngenomes2<- length(unique(strain_names2)) 
ngenomes2
#252

#add number_genomes column to get totals
#first replace blank calls with NAs
strains_only2<- apply(strains_only2, 2, function(x) gsub("^$|^ $", NA, x))

OF.gene_families_all_no_intro$number_genomes_no_intro<- rowSums(!is.na(strains_only2)) #cols that are not blank


######

###redo stats w.o introgressed strains
#recalculate n_gene_fams
sum(OF.gene_families_all_no_intro$number_genomes_no_intro == 0) #92 gene fams are exclusive to the introgressed strains
dim(OF.gene_families_all_no_intro)
OF.gene_families_all_no_intro<- OF.gene_families_all_no_intro[!OF.gene_families_all_no_intro$number_genomes_no_intro ==0,]
dim(OF.gene_families_all_no_intro)
n_gene_fams_no_intro<- nrow(OF.gene_families_all_no_intro)

#how many of the gene families are in every genome (of 252)
n_gene_fams_core_all_no_intro<- sum(OF.gene_families_all_no_intro$number_genomes_no_intro == 252)
n_gene_fams_core_all_no_intro
#3946, vs. 3,902 if considering introgressed strains. 

#that's what percent out of the total?
(n_gene_fams_core_all_no_intro*100)/n_gene_fams_no_intro
#25.93152, vs 25.48827 if considering introgressed strains 

  
#present in 95% of genomes
cutoff<- round(.95*252)
n_gene_fams_core_w95per<- sum(OF.gene_families_all_no_intro$number_genomes_no_intro >= cutoff)
n_gene_fams_core_w95per
#8,890, vs 8,866 if considering introgressed strains

#that's what percent out of the total?
(n_gene_fams_core_w95per*100)/n_gene_fams_no_intro
#58.4215%, vs 57.91365% if considering introgressed strains

#how many of the gene families are singletons (accessory)
n_gene_fams_singletons<- sum(as.numeric(OF.gene_families_all_no_intro$number_genomes_no_intro) == 1)
n_gene_fams_singletons
#2,081, vs. 2,109 if considering introgressed strains

#that's what percent out of the total?
(n_gene_fams_singletons*100)/n_gene_fams_no_intro
#13.67549% vs. 13.77621% if considering introgressed strains

#how many accessory?)
n_accessory<- n_gene_fams_no_intro - (n_gene_fams_singletons + n_gene_fams_core_w95per)
n_accessory
#4,246 vs. 4,334 if considering introgressed strains

#that's what percent out of the total?
(n_accessory*100)/n_gene_fams_no_intro
#27.903% vs. 28.31014% if considering introgressed strains

#get average per genome 
n_accessory / ngenomes2
#16.84921 accessory genes per genome
