#A fum pan genome analysis, look as interpro annotation results for Clade specific gene families
#started: 12.May.2021
#last updated: 12.Nov.2021

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
library(topGO)
library(scales)
library(cowplot)
library(biostat)

setwd("~/Desktop/Project_Afum_pangenome_2")

anno_all<-as.data.frame(fread("Af.pan.genome.iprout.tsv", header = FALSE, sep = "\t", strip.white=TRUE, quote="", fill = TRUE))

column_names <-c("OF_fam",
                 "MD5",	
                 "Sequence_length",	
                 "Analysis",	
                 "Signature_accession",	
                 "Signature_description",	
                 "start_location",	
                 "stop_location",	
                 "Score",	
                 "Status",	
                 "Date",	
                 "InterPro_annotations_accession",	
                 "InterPro_annotations_description",	
                 "GO_annotations", 
                 "Pathways_annotations")

setnames(anno_all, column_names)
length(unique(anno_all$InterPro_annotations_accession))
length(unique(anno_all$OF_fam)) #13,869 (out of 15,309 gene families) had some kind of annotation.


#read in OF gene fams
#setwd("~/Desktop/Project_Afum_pangenome_2")
#setwd("~/Desktop/Project_Afum_pangenome_legacy")
#PIRATE.gene_families<-as.data.frame(fread("PIRATE.gene_families.ordered_18May2021.tsv")) 

OF.gene_families<-as.data.frame(fread("Orthogroups.tsv")) 
OF.unassigned<- as.data.frame(fread("Orthogroups_UnassignedGenes.tsv")) #this is where your singletons live

#combine OF families with singletons
OF.gene_families<- rbind(OF.gene_families, OF.unassigned)

length(setdiff(OF.gene_families$Orthogroup,anno_all$OF_fam)) #yes, 1440 don't have annotations (compared to 1402 in PIRATE)

#what genes are missing annotations? Write for later use
missing<- data.frame(gene_fam=setdiff(OF.gene_families$Orthogroup,anno_all$OF_fam))
#write.table(missing, "missing_from_IPR_anno_OF_260_genomes.txt", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#add extra col of totals
cols_to_exclude<- "Orthogroup"
strains_only<- OF.gene_families[,!names(OF.gene_families) %in% cols_to_exclude]
#add NAs where blank
strains_only<- apply(strains_only, 2, function(x) gsub("^$|^ $", NA, x))
OF.gene_families$number_genomes<- rowSums(!is.na(strains_only)) #cols that are not blank

#subset to core, accessory and singletons to get gene lists
core<- OF.gene_families[OF.gene_families$number_genomes >= 247,]
accessory1<- OF.gene_families[OF.gene_families$number_genomes > 1,]
accessory<- accessory1[accessory1$number_genomes < 247,]
singleton<- OF.gene_families[OF.gene_families$number_genomes == 1,]

#make sure these numbers are correct 
length(unique(core$Orthogroup)) #8,866
length(unique(accessory$Orthogroup)) #4,334
length(unique(singleton$Orthogroup)) #2,109. -yes these all match up.

#is product length inflating singletons? - NOTE GO BACK AND DO THIS BY PILLING LENGTH FOR ALL SEQS IN BASH
#mean(core$`average_length(bp)`)
#min(core$`min_length(bp)`)
#mean(accessory$`average_length(bp)`)
#min(accessory$`min_length(bp)`)
#mean(singleton$`average_length(bp)`) #no. 
#min(singleton$`min_length(bp)`)

#is length responsible for missing annotations?
#missing_anno<- PIRATE.gene_families[PIRATE.gene_families$gene_family %in% missing$gene_fam,]
#mean(missing_anno$`average_length(bp)`) #yeah, maybe - they are about half the length. Ave. of 485 bp. 
#nrow(missing_anno)

#print
#write.table(core$gene_family, "core.txt", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
#write.table(accessory$gene_family, "accessory.txt", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
#write.table(singleton$gene_family, "singleton.txt", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

##subset to clade 1, 2 and 3 to get gene lists
#read in gene list files exclusive to ea. clade
exclusive_to_clade1<-as.data.frame(fread("exclusive_to_clade1_all.txt"))
exclusive_to_clade2<-as.data.frame(fread("exclusive_to_clade2_all.txt"))
exclusive_to_clade3<-as.data.frame(fread("exclusive_to_clade3_all.txt"))


#subset the pirate files to just these genes
clade_1<- OF.gene_families[OF.gene_families$Orthogroup %in% exclusive_to_clade1$OF_gene_fam,]
clade_2<- OF.gene_families[OF.gene_families$Orthogroup %in% exclusive_to_clade2$OF_gene_fam,]
clade_3<- OF.gene_families[OF.gene_families$Orthogroup %in% exclusive_to_clade3$OF_gene_fam,]
nrow(clade_1)
nrow(clade_2)
nrow(clade_3)

#check that there's no difference in these sets
length(setdiff(exclusive_to_clade1$OF_gene_fam,OF.gene_families$Orthogroup))
length(setdiff(exclusive_to_clade2$OF_gene_fam, OF.gene_families$Orthogroup))
length(setdiff(exclusive_to_clade3$OF_gene_fam,OF.gene_families$Orthogroup))
#looks good


###Subset the above six gene lists (core, accessory, singelton) and (Clade 1,2,3) gene lists from the interpro annotation file


####
#subset core/accessory/singleton
core_subset<- anno_all[anno_all$OF_fam %in% core$Orthogroup,]
n_core<- length(unique(core$Orthogroup))
length(unique(core_subset$OF_fam)) #got most of these: 8,711 out of 8,866
#what's missing?
missing_core<- length(setdiff(core$Orthogroup,core_subset$OF_fam)) #155 core genes missing any functional annotation
missing_core
setdiff(core$Orthogroup,core_subset$OF_fam) #These look to be numerically random 
#percent of total
core_miss_percent<- missing_core*100 / n_core
core_miss_percent #1.748252% of gene fams missing an annotation in core

accessory_subset<- anno_all[anno_all$OF_fam %in% accessory$Orthogroup,]
n_accessory<- length(unique(accessory$Orthogroup))
length(unique(accessory_subset$OF_fam)) #we got 3,444 gene fams annotations out of 4,334
mising_accessory<- length(setdiff(accessory$Orthogroup,accessory_subset$OF_fam)) #890 gene fams missing any functional annotation
mising_accessory
#percent of total
accessory_miss_percent<- mising_accessory*100 / n_accessory
accessory_miss_percent #20.5353% of gene fams missing an annotation in accessory

singleton_subset<- anno_all[anno_all$OF_fam %in% singleton$Orthogroup,]
n_singleton<- length(unique(singleton$Orthogroup))
length(unique(singleton_subset$OF_fam)) #we got 1,714 gene fams with annotations out of 2,109
missing_singleton<-length(setdiff(singleton$Orthogroup,singleton_subset$OF_fam)) #395 gene fams missing any functional annotation
missing_singleton
#percent of total
singleton_miss_percent<- missing_singleton*100 / n_singleton
singleton_miss_percent


#subset clade1/2/3
clade_1_subset<- anno_all[anno_all$OF_fam %in% clade_1$Orthogroup,]
n_clade1<- length(unique(clade_1$Orthogroup))
length(unique(clade_1_subset$OF_fam)) #we got 1024 gene fams with annotations out of 1024 gene fams
missing_clade1<- length(setdiff(clade_1$Orthogroup,clade_1_subset$OF_fam)) # 232 gene fams missing any functional annotations
missing_clade1
#percent of total
clade1_miss_percent<- missing_clade1*100 / n_clade1
clade1_miss_percent


clade_2_subset<- anno_all[anno_all$OF_fam %in% clade_2$Orthogroup,]
n_clade2<- length(unique(clade_2$Orthogroup))
length(unique(clade_2_subset$OF_fam)) #we got 77 gene fams with annotations out of 94 gene fams
missing_clade2<- length(setdiff(clade_2$Orthogroup,clade_2_subset$OF_fam)) # 18 gene fams missing any functional annotations
missing_clade2
#percent of total
clade2_miss_percent<- missing_clade2*100 / n_clade2
clade2_miss_percent


clade_3_subset<- anno_all[anno_all$OF_fam %in% clade_3$Orthogroup,]
n_clade3<- length(unique(clade_3$Orthogroup))
length(unique(clade_3_subset$OF_fam)) #we got 70 gene fams with annotations out of 115 gene fams
missing_clade3<- length(setdiff(clade_3$Orthogroup,clade_3_subset$OF_fam)) # 45 gene fams missing any functional annotations
missing_clade3
#percent of total
clade3_miss_percent<- missing_clade3*100 / n_clade3
clade3_miss_percent



#subset core/accessory/singleton to capture only gene families with ipr annotations
with_annotations_core<- core_subset[core_subset$InterPro_annotations_accession!="-",]
dim(core_subset) #out of 124089 domains 
dim(with_annotations_core) #46982 with annotations
with_annotations_accessory<- accessory_subset[accessory_subset$InterPro_annotations_accession!="-",]
dim(with_annotations_accessory) #9843 with annotations
with_annotations_singleton<- singleton_subset[singleton_subset$InterPro_annotations_accession!="-",]
dim(with_annotations_singleton) #2416 with annotations

#isolate IPR terms for each group
IPROfInterestCore<- with_annotations_core$InterPro_annotations_accession
IPROfInterestAcessory<- with_annotations_accessory$InterPro_annotations_accession
IPROfInterestSingelton<- with_annotations_singleton$InterPro_annotations_accession

#subset clade_1/clade_2/clade_3 to capture only gene families with ipr annotations
with_annotations_clade_1<- clade_1_subset[clade_1_subset$InterPro_annotations_accession!="-",]
dim(clade_1_subset) #out of 7497 domains 
dim(with_annotations_clade_1) #2771 with annotations
with_annotations_clade_2<- clade_2_subset[clade_2_subset$InterPro_annotations_accession!="-",]
dim(with_annotations_clade_2) #193 with annotations
with_annotations_clade_3<- clade_3_subset[clade_3_subset$InterPro_annotations_accession!="-",]
dim(with_annotations_clade_3) #90 with annotations

#isolate IPR terms for each clade
IPROfInterestClade_1<- with_annotations_clade_1$InterPro_annotations_accession
IPROfInterestClade_2<- with_annotations_clade_2$InterPro_annotations_accession
IPROfInterestClade_3<- with_annotations_clade_3$InterPro_annotations_accession

#isolate IPR terms to use as the background (gene universe)
just_interpro<- anno_all[anno_all$InterPro_annotations_accession!= "-", ]
nrow(just_interpro) #59241 domains with ipr annotations
just_GO<- just_interpro[just_interpro$GO_annotations!= "-", ]
nrow(just_GO) #32857 of the 59241 have GO annotations
length(just_GO$InterPro_annotations_accession)
geneUniverse <- just_GO$InterPro_annotations_accession


#tell topGO where to look for the genes of interest - here 1 is of interest, and zero is not

#core/accessory/singleton
geneList_core <- factor(as.integer(geneUniverse %in% IPROfInterestCore))
names(geneList_core) <- geneUniverse 
geneList_accessory <- factor(as.integer(geneUniverse %in% IPROfInterestAcessory))
names(geneList_accessory) <- geneUniverse 
geneList_singleton <- factor(as.integer(geneUniverse %in% IPROfInterestSingelton))
names(geneList_singleton) <- geneUniverse 


#clade1/2/3
geneList_clade_1 <- factor(as.integer(geneUniverse %in% IPROfInterestClade_1))
names(geneList_clade_1) <- geneUniverse 
geneList_clade_2 <- factor(as.integer(geneUniverse %in% IPROfInterestClade_2))
names(geneList_clade_2) <- geneUniverse 
geneList_clade_3 <- factor(as.integer(geneUniverse %in% IPROfInterestClade_3))
names(geneList_clade_3) <- geneUniverse 


#set database
geneID2GO_new<- as.list(just_GO$GO_annotations)
length(geneID2GO_new)
names(geneID2GO_new)<- (just_GO$InterPro_annotations_accession)


#to put the data into an object of type 'topGOdata'. 
#This will contain the list of genes of interest, the GO annotations, and the GO hierarchy.
#NOTE The 'ontology' argument can be set to 'BP' (biological process), 'MF' (molecular function), or 'CC' (cellular component).
GOdata_core_MF <- new("topGOdata", ontology = "MF", allGenes = geneList_core,
                      annot = annFUN.gene2GO, gene2GO = geneID2GO_new, 
                      nodeSize = 6)
GOdata_accessory_MF <- new("topGOdata", ontology = "MF", allGenes = geneList_accessory,
                           annot = annFUN.gene2GO, gene2GO = geneID2GO_new,
                           nodeSize = 6)
GOdata_singleton_MF <- new("topGOdata", ontology = "MF", allGenes = geneList_singleton,
                           annot = annFUN.gene2GO, gene2GO = geneID2GO_new,
                           nodeSize = 6)

#run fisher exact test
resultFisher_core_MF <- runTest(GOdata_core_MF, algorithm="weight01", statistic="fisher")
resultFisher_core_MF
resultFisher_accessory_MF <- runTest(GOdata_accessory_MF, algorithm="weight01", statistic="fisher")
resultFisher_accessory_MF
resultFisher_singleton_MF <- runTest(GOdata_singleton_MF, algorithm="weight01", statistic="fisher")
resultFisher_singleton_MF

#get top 7 results 
allRes_core_MF <- GenTable(GOdata_core_MF, classicFisher = resultFisher_core_MF, orderBy = "resultFisher", ranksOf = "classicFisher", numChar=1000, topNodes = 10)
allRes_core_MF
allRes_df_core_MF<- as.data.frame(allRes_core_MF)

allRes_accessory_MF <- GenTable(GOdata_accessory_MF, classicFisher = resultFisher_accessory_MF, orderBy = "resultFisher", ranksOf = "classicFisher", numChar=1000, topNodes = 10)
allRes_accessory_MF
allRes_df_accessory_MF<- as.data.frame(allRes_accessory_MF)

allRes_singleton_MF <- GenTable(GOdata_singleton_MF, classicFisher = resultFisher_singleton_MF, orderBy = "resultFisher", ranksOf = "classicFisher", numChar=1000, topNodes = 10)
allRes_singleton_MF
allRes_df_singleton_MF<- as.data.frame(allRes_singleton_MF)

#fix p-vlaues to make numeric (core)
allRes_core_MF$classicFisher
allRes_core_MF[]<- lapply(allRes_core_MF, gsub, pattern='< ', replacement="")
allRes_core_MF$classicFisher<- as.numeric(allRes_core_MF$classicFisher)
#fix p-vlaues to make numeric (accessory)
allRes_accessory_MF$classicFisher
allRes_accessory_MF[]<- lapply(allRes_accessory_MF, gsub, pattern='< ', replacement="")
allRes_accessory_MF$classicFisher<- as.numeric(allRes_accessory_MF$classicFisher)
#fix p-vlaues to make numeric (singleton)
allRes_singleton_MF$classicFisher
allRes_singleton_MF[]<- lapply(allRes_singleton_MF, gsub, pattern='< ', replacement="")
allRes_singleton_MF$classicFisher<- as.numeric(allRes_singleton_MF$classicFisher)

#subset to only significant results
allRes_core_MF<- allRes_core_MF[allRes_core_MF$classicFisher<0.05,]
allRes_accessory_MF<- allRes_accessory_MF[allRes_accessory_MF$classicFisher<0.05,]
allRes_singleton_MF<- allRes_singleton_MF[allRes_singleton_MF$classicFisher<0.05,]

#set Term as factor for graphing
allRes_core_MF$Term<- as.factor(allRes_core_MF$Term)
allRes_accessory_MF$Term<- as.factor(allRes_accessory_MF$Term)
allRes_singleton_MF$Term<- as.factor(allRes_singleton_MF$Term)

##run the same for biological process
GOdata_core_BP <- new("topGOdata", ontology = "BP", allGenes = geneList_core,
                      annot = annFUN.gene2GO, gene2GO = geneID2GO_new, 
                      nodeSize = 6)

GOdata_accessory_BP <- new("topGOdata", ontology = "BP", allGenes = geneList_accessory,
                           annot = annFUN.gene2GO, gene2GO = geneID2GO_new,
                           nodeSize = 6)
GOdata_singleton_BP <- new("topGOdata", ontology = "BP", allGenes = geneList_singleton,
                           annot = annFUN.gene2GO, gene2GO = geneID2GO_new,
                           nodeSize = 6)

#run fisher exact test
resultFisher_core_BP <- runTest(GOdata_core_BP, algorithm="weight01", statistic="fisher")
resultFisher_core_BP
resultFisher_accessory_BP <- runTest(GOdata_accessory_BP, algorithm="weight01", statistic="fisher")
resultFisher_accessory_BP
resultFisher_singleton_BP <- runTest(GOdata_singleton_BP, algorithm="weight01", statistic="fisher")
resultFisher_singleton_BP

#get top 10 results 
allRes_core_BP <- GenTable(GOdata_core_BP, classicFisher = resultFisher_core_BP, orderBy = "resultFisher", ranksOf = "classicFisher", numChar=1000, topNodes = 10)
allRes_core_BP
allRes_df_core_BP<- as.data.frame(allRes_core_BP)


allRes_accessory_BP <- GenTable(GOdata_accessory_BP, classicFisher = resultFisher_accessory_BP, orderBy = "resultFisher", ranksOf = "classicFisher", numChar=1000, topNodes = 10)
allRes_accessory_BP
allRes_df_accessory_BP<- as.data.frame(allRes_accessory_BP)

allRes_singleton_BP <- GenTable(GOdata_singleton_BP, classicFisher = resultFisher_singleton_BP, orderBy = "resultFisher", ranksOf = "classicFisher", numChar=1000, topNodes = 10)
allRes_singleton_BP
allRes_df_singleton_BP<- as.data.frame(allRes_singleton_BP)

#fix p-vlaues to make numeric (core)
allRes_core_BP$classicFisher
allRes_core_BP[]<- lapply(allRes_core_BP, gsub, pattern='< ', replacement="")
allRes_core_BP$classicFisher<- as.numeric(allRes_core_BP$classicFisher)
#fix p-vlaues to make numeric (accessory)
allRes_accessory_BP$classicFisher
allRes_accessory_BP[]<- lapply(allRes_accessory_BP, gsub, pattern='< ', replacement="")
allRes_accessory_BP$classicFisher<- as.numeric(allRes_accessory_BP$classicFisher)
#fix p-vlaues to make numeric (singleton)
allRes_singleton_BP$classicFisher
allRes_singleton_BP[]<- lapply(allRes_singleton_BP, gsub, pattern='< ', replacement="")
allRes_singleton_BP$classicFisher<- as.numeric(allRes_singleton_BP$classicFisher)

#subset to only significant results
allRes_core_BP<- allRes_core_BP[allRes_core_BP$classicFisher<0.05,]
allRes_accessory_BP<- allRes_accessory_BP[allRes_accessory_BP$classicFisher<0.05,]
allRes_singleton_BP<- allRes_singleton_BP[allRes_singleton_BP$classicFisher<0.05,]

#set Term as factor for graphing
allRes_core_BP$Term<- as.factor(allRes_core_BP$Term)
allRes_accessory_BP$Term<- as.factor(allRes_accessory_BP$Term)
allRes_singleton_BP$Term<- as.factor(allRes_singleton_BP$Term)


###Clade1/2/3
GOdata_clade_1_MF <- new("topGOdata", ontology = "MF", allGenes = geneList_clade_1,
                         annot = annFUN.gene2GO, gene2GO = geneID2GO_new, 
                         nodeSize = 6)
GOdata_clade_2_MF <- new("topGOdata", ontology = "MF", allGenes = geneList_clade_2,
                         annot = annFUN.gene2GO, gene2GO = geneID2GO_new,
                         nodeSize = 6)
GOdata_clade_3_MF <- new("topGOdata", ontology = "MF", allGenes = geneList_clade_3,
                         annot = annFUN.gene2GO, gene2GO = geneID2GO_new,
                         nodeSize = 6)

#run fisher exact test
resultFisher_clade_1_MF <- runTest(GOdata_clade_1_MF, algorithm="weight01", statistic="fisher")
resultFisher_clade_1_MF
resultFisher_clade_2_MF <- runTest(GOdata_clade_2_MF, algorithm="weight01", statistic="fisher")
resultFisher_clade_2_MF
resultFisher_clade_3_MF <- runTest(GOdata_clade_3_MF, algorithm="weight01", statistic="fisher")
resultFisher_clade_3_MF

#get top 7 results 
allRes_clade_1_MF <- GenTable(GOdata_clade_1_MF, classicFisher = resultFisher_clade_1_MF, orderBy = "resultFisher", ranksOf = "classicFisher", numChar=1000, topNodes = 10)
allRes_clade_1_MF
allRes_df_clade_1_MF<- as.data.frame(allRes_clade_1_MF)

allRes_clade_2_MF <- GenTable(GOdata_clade_2_MF, classicFisher = resultFisher_clade_2_MF, orderBy = "resultFisher", ranksOf = "classicFisher", numChar=1000, topNodes = 10)
allRes_clade_2_MF
allRes_df_clade_2_MF<- as.data.frame(allRes_clade_2_MF)

allRes_clade_3_MF <- GenTable(GOdata_clade_3_MF, classicFisher = resultFisher_clade_3_MF, orderBy = "resultFisher", ranksOf = "classicFisher", numChar=1000, topNodes = 10)
allRes_clade_3_MF
allRes_df_clade_3_MF<- as.data.frame(allRes_clade_3_MF)

#fix p-vlaues to make numeric (clade_1)
allRes_clade_1_MF$classicFisher
allRes_clade_1_MF[]<- lapply(allRes_clade_1_MF, gsub, pattern='< ', replacement="")
allRes_clade_1_MF$classicFisher<- as.numeric(allRes_clade_1_MF$classicFisher)
#fix p-vlaues to make numeric (clade_2)
allRes_clade_2_MF$classicFisher
allRes_clade_2_MF[]<- lapply(allRes_clade_2_MF, gsub, pattern='< ', replacement="")
allRes_clade_2_MF$classicFisher<- as.numeric(allRes_clade_2_MF$classicFisher)
#fix p-vlaues to make numeric (clade_3)
allRes_clade_3_MF$classicFisher
allRes_clade_3_MF[]<- lapply(allRes_clade_3_MF, gsub, pattern='< ', replacement="")
allRes_clade_3_MF$classicFisher<- as.numeric(allRes_clade_3_MF$classicFisher)

#subset to only significant results
allRes_clade_1_MF<- allRes_clade_1_MF[allRes_clade_1_MF$classicFisher<0.05,]
allRes_clade_2_MF<- allRes_clade_2_MF[allRes_clade_2_MF$classicFisher<0.05,]
allRes_clade_3_MF<- allRes_clade_3_MF[allRes_clade_3_MF$classicFisher<0.05,]

#set Term as factor for graphing
allRes_clade_1_MF$Term<- as.factor(allRes_clade_1_MF$Term)
allRes_clade_2_MF$Term<- as.factor(allRes_clade_2_MF$Term)
allRes_clade_3_MF$Term<- as.factor(allRes_clade_3_MF$Term)

##run the same for biological process
GOdata_clade_1_BP <- new("topGOdata", ontology = "BP", allGenes = geneList_clade_1,
                         annot = annFUN.gene2GO, gene2GO = geneID2GO_new, 
                         nodeSize = 6)

GOdata_clade_2_BP <- new("topGOdata", ontology = "BP", allGenes = geneList_clade_2,
                         annot = annFUN.gene2GO, gene2GO = geneID2GO_new,
                         nodeSize = 6)
GOdata_clade_3_BP <- new("topGOdata", ontology = "BP", allGenes = geneList_clade_3,
                         annot = annFUN.gene2GO, gene2GO = geneID2GO_new,
                         nodeSize = 6)

#run fisher exact test
resultFisher_clade_1_BP <- runTest(GOdata_clade_1_BP, algorithm="weight01", statistic="fisher")
resultFisher_clade_1_BP
resultFisher_clade_2_BP <- runTest(GOdata_clade_2_BP, algorithm="weight01", statistic="fisher")
resultFisher_clade_2_BP
resultFisher_clade_3_BP <- runTest(GOdata_clade_3_BP, algorithm="weight01", statistic="fisher")
resultFisher_clade_3_BP

#get top 10 results 
allRes_clade_1_BP <- GenTable(GOdata_clade_1_BP, classicFisher = resultFisher_clade_1_BP, orderBy = "resultFisher", ranksOf = "classicFisher", numChar=1000, topNodes = 10)
allRes_clade_1_BP
allRes_df_clade_1_BP<- as.data.frame(allRes_clade_1_BP)

allRes_clade_2_BP <- GenTable(GOdata_clade_2_BP, classicFisher = resultFisher_clade_2_BP, orderBy = "resultFisher", ranksOf = "classicFisher", numChar=1000, topNodes = 10)
allRes_clade_2_BP
allRes_df_clade_2_BP<- as.data.frame(allRes_clade_2_BP)

allRes_clade_3_BP <- GenTable(GOdata_clade_3_BP, classicFisher = resultFisher_clade_3_BP, orderBy = "resultFisher", ranksOf = "classicFisher", numChar=1000, topNodes = 10)
allRes_clade_3_BP
allRes_df_clade_3_BP<- as.data.frame(allRes_clade_3_BP)


#fix p-vlaues to make numeric (clade_1)
allRes_clade_1_BP$classicFisher
allRes_clade_1_BP[]<- lapply(allRes_clade_1_BP, gsub, pattern='< ', replacement="")
allRes_clade_1_BP$classicFisher<- as.numeric(allRes_clade_1_BP$classicFisher)
#fix p-vlaues to make numeric (clade_2)
allRes_clade_2_BP$classicFisher
allRes_clade_2_BP[]<- lapply(allRes_clade_2_BP, gsub, pattern='< ', replacement="")
allRes_clade_2_BP$classicFisher<- as.numeric(allRes_clade_2_BP$classicFisher)
#fix p-vlaues to make numeric (clade_3)
allRes_clade_3_BP$classicFisher
allRes_clade_3_BP[]<- lapply(allRes_clade_3_BP, gsub, pattern='< ', replacement="")
allRes_clade_3_BP$classicFisher<- as.numeric(allRes_clade_3_BP$classicFisher)

#subset to only significant results
allRes_clade_1_BP<- allRes_clade_1_BP[allRes_clade_1_BP$classicFisher<0.05,]
allRes_clade_2_BP<- allRes_clade_2_BP[allRes_clade_2_BP$classicFisher<0.05,]
allRes_clade_3_BP<- allRes_clade_3_BP[allRes_clade_3_BP$classicFisher<0.05,]

#set Term as factor for graphing
allRes_clade_1_BP$Term<- as.factor(allRes_clade_1_BP$Term)
allRes_clade_2_BP$Term<- as.factor(allRes_clade_2_BP$Term)
allRes_clade_3_BP$Term<- as.factor(allRes_clade_3_BP$Term)



##get significant genes (IPR doms) contributing to each GO term 

#Clade1 
Significant_IPR_in_sig_GO_Clade1 <- sapply(allRes_clade_1_BP$GO.ID, function(x)
{
  genes<-genesInTerm(GOdata_clade_1_BP, x) 
  genes[[1]][genes[[1]] %in% IPROfInterestClade_1]
})
Significant_IPR_in_sig_GO_Clade1

#Clade2
Significant_IPR_in_sig_GO_Clade2 <- sapply(allRes_clade_2_BP$GO.ID, function(x)
{
  genes<-genesInTerm(GOdata_clade_2_BP, x) 
  genes[[1]][genes[[1]] %in% IPROfInterestClade_2]
})
Significant_IPR_in_sig_GO_Clade2

#Clade3
Significant_IPR_in_sig_GO_Clade3 <- sapply(allRes_clade_3_BP$GO.ID, function(x)
{
  genes<-genesInTerm(GOdata_clade_3_BP, x) 
  genes[[1]][genes[[1]] %in% IPROfInterestClade_3]
})


                                              
#print these to file for supplemental 
C1_list_a<- paste(names(Significant_IPR_in_sig_GO_Clade1),Significant_IPR_in_sig_GO_Clade1,sep="=")
C1_list_b<- strsplit(C1_list_a, "=")
C1_list_c<- t(as.data.frame(C1_list_b))
C1_list_d<- gsub("c\\(", "", C1_list_c)
C1_list_e<- data.frame(gsub("\\)", "", C1_list_d))
C1_list_f<- cbind(allRes_clade_1_BP, sig_contrib_IPR_doms=C1_list_e$X2)

C2_list_a<- paste(names(Significant_IPR_in_sig_GO_Clade2),Significant_IPR_in_sig_GO_Clade2,sep="=")
C2_list_b<- strsplit(C2_list_a, "=")
C2_list_c<- t(as.data.frame(C2_list_b))
C2_list_d<- gsub("c\\(", "", C2_list_c)
C2_list_e<- data.frame(gsub("\\)", "", C2_list_d))
C2_list_f<- cbind(allRes_clade_2_BP, sig_contrib_IPR_doms=C2_list_e$X2)

C3_list_a<- paste(names(Significant_IPR_in_sig_GO_Clade3),Significant_IPR_in_sig_GO_Clade3,sep="=")
C3_list_b<- strsplit(C3_list_a, "=")
C3_list_c<- t(as.data.frame(C3_list_b))
C3_list_d<- gsub("c\\(", "", C3_list_c)
C3_list_e<- data.frame(gsub("\\)", "", C3_list_d))
C3_list_f<- cbind(allRes_clade_3_BP, sig_contrib_IPR_doms=C3_list_e$X2)

#write.table(C1_list_f, "significat_IPR_domains_Clade1", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(C2_list_f, "significat_IPR_domains_Clade2", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(C3_list_f, "significat_IPR_domains_Clade3", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



###plot dot plots of GO enrichment 
#plot function
dot_plot_fun <- function(input_df,name,color.low, color.high){
input_df %>%
  arrange(classicFisher) %>%
  mutate(Term=factor(Term, levels=Term)) %>% 
  ggplot(aes(x = Term, 
             y = -(classicFisher), 
             size = as.numeric(Annotated), 
             fill = -(classicFisher))) +
  labs(fill = "Enrichment (p-value)", size = "n domains assigned GO term")+
  #expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = color.low, high = color.high) +
  
  xlab('') + ylab("P-value (Fisher Exact)") +
  labs(
  title = name) +
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
    
    axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 12, face = 'bold'),
    axis.title.x = element_text(size = 12, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(.4, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 9, face = "bold"), # Text size
    title = element_text(size = 9, face = "bold")) +
  
  coord_flip()}

#NOTE:#colors
#singelton = "#316A6E",
#accessory = "#BA9141",
#core = "#6E572C")

###plot
#molecular function
#Test
#dot_plot_fun(input_df = allRes_core_MF, name = "Core MF", color.low = "#6E572C00", color.high = "#6E572C")
#ggplot2::ggsave("test_goterms.pdf",
#                device = NULL,
#                height = 8.5,
#                width = 14)

#####
#graph together
#NOTE: colors:
#clade1 = "#56326E"
#clade2 ="#ED7F6F"
#clade3 = "#ABA778"

#by core/accessory/singleton
core_mf<- dot_plot_fun(input_df = allRes_core_MF, name = "Core MF", color.low = "#6E572C00", color.high = "#6E572C")
accessory_mf<- dot_plot_fun(input_df = allRes_accessory_MF, name = "Accessory MF", color.low = "#BA914100", color.high = "#BA9141")
singleton_mf<- dot_plot_fun(input_df = allRes_singleton_MF, name = "Singleton MF", color.low = "#316A6E00", color.high = "#316A6E")

core_bp<- dot_plot_fun(input_df = allRes_core_BP, name = "Core BP", color.low = "#6E572C00", color.high = "#6E572C")
accessory_bp<- dot_plot_fun(input_df = allRes_accessory_BP, name = "Accessory BP", color.low = "#BA914100", color.high = "#BA9141")
singleton_bp<- dot_plot_fun(input_df = allRes_singleton_BP, name = "Singleton BP", color.low = "#316A6E00", color.high = "#316A6E")

#by clade
clade1_mf<- dot_plot_fun(input_df = allRes_clade_1_MF, name = "Clade 1 MF", color.low = "#56326E00", color.high = "#56326E")
clade2_mf<- dot_plot_fun(input_df = allRes_clade_2_MF, name = "Clade 2 MF", color.low = "#ED7F6F00", color.high = "#ED7F6F")
clade3_mf<- dot_plot_fun(input_df = allRes_clade_3_MF, name = "Clade 3 MF", color.low = "#ABA77800", color.high = "#ABA778")

clade1_bp<- dot_plot_fun(input_df = allRes_clade_1_BP, name = "Clade 1 BP", color.low = "#56326E00", color.high = "#56326E")
clade2_bp<- dot_plot_fun(input_df = allRes_clade_2_BP, name = "Clade 2 BP", color.low = "#ED7F6F00", color.high = "#ED7F6F")
clade3_bp<- dot_plot_fun(input_df = allRes_clade_3_BP, name = "Clade 3 BP", color.low = "#ABA77800", color.high = "#ABA778")

#together
p_bp <- plot_grid(core_bp,
                  accessory_bp,
                  singleton_bp,
                  clade1_bp,
                  clade2_bp,
                  clade3_bp,
                  labels = c("A", "D", "B", "E", "C", "F"), align = "hv", axis = "rlbt", nrow = 3, ncol =2, byrow = FALSE)
ggsave("GO_graph_BP_top10.pdf", p_bp, 
       height = 12,
       width = 20)

p_mf <- plot_grid(core_mf,
                  accessory_mf,
                  singleton_mf,
                  clade1_mf,
                  clade2_mf,
                  clade3_mf,
                  labels = c("A", "D", "B", "E", "C", "F"), align = "hv", axis = "rlbt", nrow = 3, ncol =2, byrow = FALSE)
#ggsave("GO_graph_MF.pdf", p_mf, 
#       height = 12,
#       width = 28)







##########
#Secreted Protein Analysis
#########
#core
#get Phobius SP annotations
Phobius_core<-core_subset[core_subset$Analysis == "Phobius",]
#remove gene fams that have both phobius SP and TM predictions. 
Phobius_core_SP<- Phobius_core[Phobius_core$Signature_accession == "SIGNAL_PEPTIDE" | Phobius_core$Signature_accession == "TRANSMEMBRANE", ]
Phobius_core_SP_table<- as.data.frame(table(Phobius_core_SP$OF_fam, Phobius_core_SP$Signature_accession))
Phobius_core_SP_table_wide<- 
  Phobius_core_SP_table %>%
  pivot_wider(names_from = Var2, values_from = Freq)
#get gene fams with at least 1 SP, and 0 TM doms.   
Phobius_core_SP_table_wide_no_TM1<- Phobius_core_SP_table_wide[(Phobius_core_SP_table_wide$SIGNAL_PEPTIDE >= 1 ),]
Phobius_core_SP_table_wide_no_TM<- Phobius_core_SP_table_wide_no_TM1[(Phobius_core_SP_table_wide_no_TM1$TRANSMEMBRANE == 0 ),]
unique_core_Phobius_gene_fams<- data.frame(gene_fam =Phobius_core_SP_table_wide_no_TM$Var1)
nrow(unique_core_Phobius_gene_fams) #885
#get SignalP annotations 
Signal_P_core<-core_subset[core_subset$Analysis == "SignalP_EUK",]
unique_core_SignalP_gene_fams<- data.frame(gene_fam =unique(Signal_P_core$OF_fam))
nrow(unique_core_SignalP_gene_fams) #768
#intersect- how many gene fame are annotated in both 
intersect_core<- length(intersect(unique_core_Phobius_gene_fams$gene_fam, unique_core_SignalP_gene_fams$gene_fam))
intersect_core
#592 intersect. 

###

#accessory
#get Phobius SP annotations
Phobius_accessory<-accessory_subset[accessory_subset$Analysis == "Phobius",]
#remove gene fams that have both phobius SP and TM predictions. 
Phobius_accessory_SP<- Phobius_accessory[Phobius_accessory$Signature_accession == "SIGNAL_PEPTIDE" | Phobius_accessory$Signature_accession == "TRANSMEMBRANE", ]
Phobius_accessory_SP_table<- as.data.frame(table(Phobius_accessory_SP$OF_fam, Phobius_accessory_SP$Signature_accession))
Phobius_accessory_SP_table_wide<- 
  Phobius_accessory_SP_table %>%
  pivot_wider(names_from = Var2, values_from = Freq)
#get gene fams with at least 1 SP, and 0 TM doms.   
Phobius_accessory_SP_table_wide_no_TM1<- Phobius_accessory_SP_table_wide[(Phobius_accessory_SP_table_wide$SIGNAL_PEPTIDE >= 1 ),]
Phobius_accessory_SP_table_wide_no_TM<- Phobius_accessory_SP_table_wide_no_TM1[(Phobius_accessory_SP_table_wide_no_TM1$TRANSMEMBRANE == 0 ),]
unique_accessory_Phobius_gene_fams<- data.frame(gene_fam =Phobius_accessory_SP_table_wide_no_TM$Var1)
nrow(unique_accessory_Phobius_gene_fams) #296
#get SignalP annotations 
Signal_P_accessory<-accessory_subset[accessory_subset$Analysis == "SignalP_EUK",]
unique_accessory_SignalP_gene_fams<- data.frame(gene_fam =unique(Signal_P_accessory$OF_fam))
nrow(unique_accessory_SignalP_gene_fams) #224
#intersect- how many gene fame are annotated in both 
intersect_accessory<- length(intersect(unique_accessory_Phobius_gene_fams$gene_fam, unique_accessory_SignalP_gene_fams$gene_fam))
intersect_accessory
#177 intersect. 

###

#singleton
#get Phobius SP annotations
Phobius_singleton<-singleton_subset[singleton_subset$Analysis == "Phobius",]
#remove gene fams that have both phobius SP and TM predictions. 
Phobius_singleton_SP<- Phobius_singleton[Phobius_singleton$Signature_accession == "SIGNAL_PEPTIDE" | Phobius_singleton$Signature_accession == "TRANSMEMBRANE", ]
Phobius_singleton_SP_table<- as.data.frame(table(Phobius_singleton_SP$OF_fam, Phobius_singleton_SP$Signature_accession))
Phobius_singleton_SP_table_wide<- 
  Phobius_singleton_SP_table %>%
  pivot_wider(names_from = Var2, values_from = Freq)
#get gene fams with at least 1 SP, and 0 TM doms.   
Phobius_singleton_SP_table_wide_no_TM1<- Phobius_singleton_SP_table_wide[(Phobius_singleton_SP_table_wide$SIGNAL_PEPTIDE >= 1 ),]
Phobius_singleton_SP_table_wide_no_TM<- Phobius_singleton_SP_table_wide_no_TM1[(Phobius_singleton_SP_table_wide_no_TM1$TRANSMEMBRANE == 0 ),]
unique_singleton_Phobius_gene_fams<- data.frame(gene_fam =Phobius_singleton_SP_table_wide_no_TM$Var1)
nrow(unique_singleton_Phobius_gene_fams) #118
#get SignalP annotations 
Signal_P_singleton<-singleton_subset[singleton_subset$Analysis == "SignalP_EUK",]
unique_singleton_SignalP_gene_fams<- data.frame(gene_fam =unique(Signal_P_singleton$OF_fam))
nrow(unique_singleton_SignalP_gene_fams) #81
#intersect- how many gene fame are annotated in both 
intersect_singleton<- length(intersect(unique_singleton_Phobius_gene_fams$gene_fam, unique_singleton_SignalP_gene_fams$gene_fam))
intersect_singleton
#68

###

#clade_1
#get Phobius SP annotations
Phobius_clade_1<-clade_1_subset[clade_1_subset$Analysis == "Phobius",]
#remove gene fams that have both phobius SP and TM predictions. 
Phobius_clade_1_SP<- Phobius_clade_1[Phobius_clade_1$Signature_accession == "SIGNAL_PEPTIDE" | Phobius_clade_1$Signature_accession == "TRANSMEMBRANE", ]
Phobius_clade_1_SP_table<- as.data.frame(table(Phobius_clade_1_SP$OF_fam, Phobius_clade_1_SP$Signature_accession))
Phobius_clade_1_SP_table_wide<- 
  Phobius_clade_1_SP_table %>%
  pivot_wider(names_from = Var2, values_from = Freq)
#get gene fams with at least 1 SP, and 0 TM doms.   
Phobius_clade_1_SP_table_wide_no_TM1<- Phobius_clade_1_SP_table_wide[(Phobius_clade_1_SP_table_wide$SIGNAL_PEPTIDE >= 1 ),]
Phobius_clade_1_SP_table_wide_no_TM<- Phobius_clade_1_SP_table_wide_no_TM1[(Phobius_clade_1_SP_table_wide_no_TM1$TRANSMEMBRANE == 0 ),]
unique_clade_1_Phobius_gene_fams<- data.frame(gene_fam =Phobius_clade_1_SP_table_wide_no_TM$Var1)
nrow(unique_clade_1_Phobius_gene_fams) #65
#get SignalP annotations 
Signal_P_clade_1<-clade_1_subset[clade_1_subset$Analysis == "SignalP_EUK",]
unique_clade_1_SignalP_gene_fams<- data.frame(gene_fam =unique(Signal_P_clade_1$OF_fam))
nrow(unique_clade_1_SignalP_gene_fams) #45
#intersect- how many gene fame are annotated in both 
intersect_clade_1<- length(intersect(unique_clade_1_Phobius_gene_fams$gene_fam, unique_clade_1_SignalP_gene_fams$gene_fam))
intersect_clade_1
#36 intersect. 

#clade_2

#get Phobius SP annotations
Phobius_clade_2<-clade_2_subset[clade_2_subset$Analysis == "Phobius",]
#remove gene fams that have both phobius SP and TM predictions. 
Phobius_clade_2_SP<- Phobius_clade_2[Phobius_clade_2$Signature_accession == "SIGNAL_PEPTIDE" | Phobius_clade_2$Signature_accession == "TRANSMEMBRANE", ]
Phobius_clade_2_SP_table<- as.data.frame(table(Phobius_clade_2_SP$OF_fam, Phobius_clade_2_SP$Signature_accession))
Phobius_clade_2_SP_table_wide<- 
  Phobius_clade_2_SP_table %>%
  pivot_wider(names_from = Var2, values_from = Freq)
#get gene fams with at least 1 SP, and 0 TM doms.   
Phobius_clade_2_SP_table_wide_no_TM1<- Phobius_clade_2_SP_table_wide[(Phobius_clade_2_SP_table_wide$SIGNAL_PEPTIDE >= 1 ),]
Phobius_clade_2_SP_table_wide_no_TM<- Phobius_clade_2_SP_table_wide_no_TM1[(Phobius_clade_2_SP_table_wide_no_TM1$TRANSMEMBRANE == 0 ),]
unique_clade_2_Phobius_gene_fams<- data.frame(gene_fam =Phobius_clade_2_SP_table_wide_no_TM$Var1)
nrow(unique_clade_2_Phobius_gene_fams) #8
#get SignalP annotations 
Signal_P_clade_2<-clade_2_subset[clade_2_subset$Analysis == "SignalP_EUK",]
unique_clade_2_SignalP_gene_fams<- data.frame(gene_fam =unique(Signal_P_clade_2$OF_fam))
nrow(unique_clade_2_SignalP_gene_fams) #
#intersect- how many gene fame are annotated in both 
intersect_clade_2<- length(intersect(unique_clade_2_Phobius_gene_fams$gene_fam, unique_clade_2_SignalP_gene_fams$gene_fam))
intersect_clade_2
#4 intersect. 

###

#clade_3
#get Phobius SP annotations
Phobius_clade_3<-clade_3_subset[clade_3_subset$Analysis == "Phobius",]
#remove gene fams that have both phobius SP and TM predictions. 
Phobius_clade_3_SP<- Phobius_clade_3[Phobius_clade_3$Signature_accession == "SIGNAL_PEPTIDE" | Phobius_clade_3$Signature_accession == "TRANSMEMBRANE", ]
Phobius_clade_3_SP_table<- as.data.frame(table(Phobius_clade_3_SP$OF_fam, Phobius_clade_3_SP$Signature_accession))
Phobius_clade_3_SP_table_wide<- 
  Phobius_clade_3_SP_table %>%
  pivot_wider(names_from = Var2, values_from = Freq)
#get gene fams with at least 1 SP, and 0 TM doms.   
Phobius_clade_3_SP_table_wide_no_TM1<- Phobius_clade_3_SP_table_wide[(Phobius_clade_3_SP_table_wide$SIGNAL_PEPTIDE >= 1 ),]
Phobius_clade_3_SP_table_wide_no_TM<- Phobius_clade_3_SP_table_wide_no_TM1[(Phobius_clade_3_SP_table_wide_no_TM1$TRANSMEMBRANE == 0 ),]
unique_clade_3_Phobius_gene_fams<- data.frame(gene_fam =Phobius_clade_3_SP_table_wide_no_TM$Var1)
nrow(unique_clade_3_Phobius_gene_fams) #11
#get SignalP annotations 
Signal_P_clade_3<-clade_3_subset[clade_3_subset$Analysis == "SignalP_EUK",]
unique_clade_3_SignalP_gene_fams<- data.frame(gene_fam =unique(Signal_P_clade_3$OF_fam))
nrow(unique_clade_3_SignalP_gene_fams) #5
#intersect- how many gene fame are annotated in both 
intersect_clade_3<- length(intersect(unique_clade_3_Phobius_gene_fams$gene_fam, unique_clade_3_SignalP_gene_fams$gene_fam))
intersect_clade_3
#6


#prep data for plot
Totals<- data.frame(Core = length(unique(core$Orthogroup)),
                       Accessory = length(unique(accessory$Orthogroup)),
                       Singleton = length(unique(singleton$Orthogroup)),
                       Clade_1 = nrow(clade_1), 
                       Clade_2 = nrow(clade_2),
                       Clade_3 = nrow(clade_3))
Phobius<- data.frame(Core = nrow(unique_core_Phobius_gene_fams),
                    Accessory = nrow(unique_accessory_Phobius_gene_fams),
                    Singleton = nrow(unique_singleton_Phobius_gene_fams),
                    Clade_1 = nrow(unique_clade_1_Phobius_gene_fams), 
                    Clade_2 = nrow(unique_clade_2_Phobius_gene_fams),
                    Clade_3 = nrow(unique_clade_3_Phobius_gene_fams))
SignalP<- data.frame(Core = nrow(unique_core_SignalP_gene_fams),
                    Accessory = nrow(unique_accessory_SignalP_gene_fams),
                    Singleton = nrow(unique_singleton_SignalP_gene_fams),
                    Clade_1 = nrow(unique_clade_1_SignalP_gene_fams), 
                    Clade_2 = nrow(unique_clade_2_SignalP_gene_fams),
                    Clade_3 = nrow(unique_clade_3_SignalP_gene_fams))
Consensus<- data.frame(Core = intersect_core,
                    Accessory = intersect_accessory,
                    Singleton = intersect_singleton,
                    Clade_1 = intersect_clade_1, 
                    Clade_2 = intersect_clade_2,
                    Clade_3 = intersect_clade_3)

#combine data
for_fig_df<-as.data.frame(t(rbind(Totals = Totals, Phobius = Phobius, SignalP = SignalP, Consensus = Consensus)))                        
for_fig_df$group<- rownames(for_fig_df)

for_fig_df_long<- as.data.frame(for_fig_df %>% 
  pivot_longer(!group, names_to = "source", values_to = "count"))

#order
for_fig_df_long$group <- factor(for_fig_df_long$group,                           
                  levels = c("Clade_3", "Clade_2", "Clade_1", "Singleton", "Accessory", "Core"))



###run stats - are there significantly more SPs in any of the groups?
per_Phobius<- (for_fig_df$Phobius)*100 / (for_fig_df$Totals)
per_SignalP<- (for_fig_df$SignalP)*100 / (for_fig_df$Totals)
per_Consensus<- (for_fig_df$Consensus)*100 / (for_fig_df$Totals)
stats_df<- cbind(for_fig_df, per_Phobius, per_SignalP, per_Consensus)
stats_df_group<- stats_df[1:3,]
stats_df_clade<- stats_df[4:6,] 

#run Pairwise comparisons for proportions test
Phobius_test<- prop.test(x = stats_df_group$Phobius, n = stats_df_group$Totals)
Phobius_test #is significant
Phobius_test_df<- as.matrix(cbind(Phobius = stats_df_group$Phobius, Totals = stats_df_group$Totals))
rownames(Phobius_test_df)<- c("Core", "Accessory", "Singleton")
Phobius_results<- pairwise.prop.test(Phobius_test_df, p.adjust.method = "bonferroni")
make_cld(Phobius_results)

SignalP_test<- prop.test(x = stats_df_group$SignalP, n = stats_df_group$Totals)
SignalP_test #is significant
SignalP_test_df<- as.matrix(cbind(SignalP = stats_df_group$SignalP, Totals = stats_df_group$Totals))
rownames(SignalP_test_df)<- c("Core", "Accessory", "Singleton")
SignalP_results<- pairwise.prop.test(SignalP_test_df, p.adjust.method = "bonferroni")
make_cld(SignalP_results)

Consensus_test<- prop.test(x = stats_df_group$Consensus, n = stats_df_group$Totals)
Consensus_test #is significant
Consensus_test_df<- as.matrix(cbind(Consensus = stats_df_group$Consensus, Totals = stats_df_group$Totals))
rownames(Consensus_test_df)<- c("Core", "Accessory", "Singleton")
Consensus_results<- pairwise.prop.test(Consensus_test_df, p.adjust.method = "bonferroni")
make_cld(Consensus_results)

#by clade
Phobius_test<- prop.test(x = stats_df_clade$Phobius, n = stats_df_clade$Totals)
Phobius_test #is NOT significant

SignalP_test<- prop.test(x = stats_df_clade$SignalP, n = stats_df_clade$Totals)
SignalP_test #is NOT significant

Consensus_test<- prop.test(x = stats_df_clade$Consensus, n = stats_df_clade$Totals)
Consensus_test #is NOT significant



#make and set colors
#clade1 = "#56326E"
#clade2 ="#ED7F6F"
#clade3 = "#ABA778"

#singelton = "#316A6E",
#accessory = "#BA9141",
#core = "#6E572C")
# Stacked + percent
#ggplot(for_fig_df_long, aes(fill=group, y=count, x=group, alpha=source)) + 
#  geom_bar(position="fill", stat="identity", color = "black") +
#  scale_alpha_manual(values=c(1, 0.8, 0.5, 0)) +
#  coord_flip()+
#  theme_ipsum() +
#  scale_fill_manual(values = c(Core = "#6E572C", Accessory = "#BA9141", Singleton = "#316A6E", Clade_1 = "#56326E", Clade_2 = "#ED7F6F", Clade_3 = "#ABA778"))
 
#adjust counts to reflect proportions
for_fig_df$Totals_updated<- for_fig_df$Totals - ((for_fig_df$Phobius +for_fig_df$SignalP) - for_fig_df$Consensus)
for_fig_df$Phobius_updated<-  for_fig_df$Phobius - for_fig_df$Consensus
for_fig_df$SignalP_updated<-  for_fig_df$SignalP - for_fig_df$Consensus

for_fig_df_subset<- for_fig_df[,4:8]

for_fig_df_long_adjusted<- as.data.frame(for_fig_df_subset %>% 
                                  pivot_longer(!group, names_to = "source", values_to = "count"))
#order
for_fig_df_long_adjusted$group <- factor(for_fig_df_long_adjusted$group,                           
                                levels = c("Clade_3", "Clade_2", "Clade_1", "Singleton", "Accessory", "Core"))
#plot with proportions
sp<- ggplot(for_fig_df_long_adjusted, aes(fill=group, y=count, x=group, alpha=source)) + 
  geom_bar(position="fill", stat="identity", color = "black") +
  scale_alpha_manual(values=c(1, 0.8, 0.5, 0)) +
  coord_flip()+
  theme_light() +
  scale_fill_manual(values = c(Core = "#6E572C", Accessory = "#BA9141", Singleton = "#316A6E", Clade_1 = "#56326E", Clade_2 = "#ED7F6F", Clade_3 = "#ABA778"))
sp

#ggsave(file="SP.pdf",sp, device="pdf", width=6.9, height=3, units="in", scale = 1, dpi = 600)
