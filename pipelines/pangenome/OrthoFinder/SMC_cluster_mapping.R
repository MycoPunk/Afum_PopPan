#mapping of Afu (reference Af293) annotated genes to Orthofinder defined OGs was conducted using BLASTP
#GENE lists  for the 26 clusters were procured from Bignal 2016

#last updated: 2.Feb.2022

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
library(drawProteins)
library(ggplot2)
library(dplyr)
library(ggtreeExtra)
library(gtools)
sessionInfo()

#set wd
setwd("~/Desktop/Project_Afum_pangenome_3/")

#load data files - hash as needed - the rest of the script does not need editing, except to set the number of top hits to display in the tree (if there are too many results to easily visualize)
SMC_genes<- as.data.frame(fread("A_fum_SMCs.csv", header = TRUE)) #this is the SMC mapping file

OF.gene_families<-as.data.frame(fread("Orthogroups.tsv")) 
OF.unassigned<- as.data.frame(fread("Orthogroups_UnassignedGenes.tsv")) #this is where your singletons live
OFtoAfu<- as.data.frame(fread("OF_blastP_results_filtered.txt", header = FALSE)) #this is the OG to Afu designations
colnames(OFtoAfu)<- c("Afu_name", "OG_name")
Afum_grp<-read.delim("clade_map_K3_3Jan2022.txt", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE, check.names = TRUE)

#combine OF families with singletons
OF.gene_families_all<- rbind(OF.gene_families, OF.unassigned)

#fix names
names(OF.gene_families_all)<- sapply(names(OF.gene_families_all), gsub, pattern = "Aspergillus_fumigatus_", replacement = "" )
names(OF.gene_families_all)<- sapply(names(OF.gene_families_all), gsub, pattern = ".proteins", replacement = "" )

#how many strains?
#names of cols to exclude
cols_to_exclude<- "Orthogroup"

strains_only<- OF.gene_families_all[,!names(OF.gene_families_all) %in% cols_to_exclude]
strain_names<- names(strains_only)

#add number_genomes column to get totals
#first replace blank calls with NAs
strains_only<- apply(strains_only, 2, function(x) gsub("^$|^ $", NA, x))
OF.gene_families_all$number_genomes<- rowSums(!is.na(strains_only)) #cols that are not blank
#replace blanks with NAs for all
OF.gene_families_all<- data.frame(apply(OF.gene_families_all, 2, function(x) gsub("^$|^ $", NA, x)))

#create 0/1 df 
gene_fam_by_strain_w_anno<-as.data.frame(OF.gene_families_all[,c(2:(ncol(OF.gene_families_all)-1))]) 
#fix X introduced in names
names(gene_fam_by_strain_w_anno)<- sapply(names(gene_fam_by_strain_w_anno), gsub, pattern = "X", replacement = "" )

##make binary (if gene = 1, if not = 0)
#replace all NAs with 0s
gene_fam_by_strain_w_anno[is.na(gene_fam_by_strain_w_anno)] <- 0
#fill in ones
gene_fam_by_strain_w_anno_ones<- as.data.frame(replace(gene_fam_by_strain_w_anno[,1:ncol(gene_fam_by_strain_w_anno)], gene_fam_by_strain_w_anno!=0, 1))
#change to numeric
gene_fam_by_strain_w_anno_ones_num <- mutate_all(gene_fam_by_strain_w_anno_ones, function(x) as.numeric(as.character(x)))
#bind annotations
gene_fam_by_strain_w_anno_ones_num2<- cbind(OF.gene_families_all = OF.gene_families_all$Orthogroup, gene_fam_by_strain_w_anno_ones_num)


#bind Afu annotations, NA if the OG is not in Afu genes.
#load data from OG to Afu BLAST search
OFtoAfu<- as.data.frame(fread("OF_blastP_results_filtered.txt", header = FALSE)) #this is the OG to Afu designations

#assign col names
fmt6_names<- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
colnames(OFtoAfu)<- fmt6_names

#attach Afu annotations
colnames(gene_fam_by_strain_w_anno_ones_num2)[1] <- "OF_gene_fam"
gene_fam_by_strain_w_anno_ones_num2$Af293_reference_gene_name<- OFtoAfu$qseqid[match(gene_fam_by_strain_w_anno_ones_num2$OF_gene_fam, OFtoAfu$sseqid)]

#subset large df to only Afu genes in the BGCs
SMC_genes_all_counts<- gene_fam_by_strain_w_anno_ones_num2[gene_fam_by_strain_w_anno_ones_num2$Af293_reference_gene_name %in% SMC_genes$Af293_gene_name,]

missing_genes<- setdiff(SMC_genes$Af293_gene_name, SMC_genes_all_counts$Af293_reference_gene_name)
#there are 10 genes that are missing (probably because they're co-clustering in other gene fams)

#make presence / absence graph of abundance 

#load tree data
tree <- read.tree("Afum_260_iq_tree_newick_v2.tre")
#remove the reference from the tree:
tree<- drop.tip(tree, "Af293-REF", trim.internal = TRUE, subtree = FALSE,
                   root.edge = 0, rooted = is.rooted(tree), collapse.singles = TRUE,
                   interactive = FALSE)
#root tree
tree<- root(tree, node = 505, resolve.root = TRUE,
               interactive = FALSE, edgelabel = FALSE)
#function to rename tips
rename.tips <- function(phy, old_names, new_names) {
  mpos <- match(old_names,phy$tip.label)
  phy$tip.label[mpos] <- new_names
  return(phy)
}

tree<- rename.tips(tree, old_names = Afum_grp$name_pop_genome_new, new_names = Afum_grp$name_to_use_in_paper)
tree$tip.label

#set colors by group for clade
my_cols_me <- c(clade1 = "#56326E",
                clade2 ="#ED7F6F",
                clade3 = "#ABA778")
grA_me<- split(Afum_grp$name_pop_genome_new, Afum_grp$clade)
tree_grA_me <- ggtree::groupOTU(tree, grA_me)
str(tree_grA_me)
levels(attributes(tree_grA_me)$group) 
levels(attributes(tree_grA_me)$group)[1] <- "1"

attributes(tree_grA_me)$group <- factor(x = attributes(tree_grA_me)$group, 
                                        levels = c("1", "2", "3"))
names(my_cols_me) <- levels(attributes(tree_grA_me)$group)


#plot tree
tree_plot_me <- 
  ggtree(tr = tree_grA_me, 
         mapping = aes(color = group), 
         layout  = 'rectangular') + 
  geom_treescale(x=5, y=NULL, color = NA) +
  # set line thickness
  # adjust coloring of main groups
  scale_color_manual(name = 'Clade', values = my_cols_me) + 
  theme(legend.title=element_text(size=9), # The title of legend 
        legend.text=element_text(size=7))+
  guides(color = guide_legend(override.aes = list(size = 4)))

#asthetics
tree_plot_me<- tree_plot_me +geom_tiplab(size = 0, align = TRUE, linesize = .25, linetype = 3)
tree_plot_me

#plot genes onto tree
#subset to only counts
all_vars<- SMC_genes_all_counts[,2:261]
#change to presence / absence
all_vars_presence<- data.frame(sapply(all_vars, gsub, pattern = "1", replacement = "present"))
all_vars_presence_absence<- data.frame(sapply(all_vars_presence, gsub, pattern = "0", replacement = "absent"))
#fix col names X's introduced in last step because R likes to make life difficult
colnames(all_vars_presence_absence)<- colnames(all_vars)
#add rowames
rownames(all_vars_presence_absence)<- SMC_genes_all_counts$Af293_reference_gene_name
#transpose for graphing
all_vars_presence_absence_t<- data.frame(t(all_vars_presence_absence))

#fix names to match OG names
row.names(all_vars_presence_absence_t) <- Afum_grp$name_to_use_in_paper[match(row.names(all_vars_presence_absence_t), Afum_grp$OF_name)]

#get counts
rownames(all_vars)<- SMC_genes_all_counts$Af293_reference_gene_name
rowSums(all_vars)

#if sums are in all, remove and don't map 
#all_vars$sums_temp<- rowSums(all_vars)
#all_vars_subset<- all_vars[all_vars$sums_temp < 250,]
#dim(all_vars)
#dim(all_vars_subset)
#subset input file 
#all_vars_presence_absence_t_subset<- all_vars_presence_absence_t[,colnames(all_vars_presence_absence_t) %in% rownames(all_vars_subset)]
#dim(all_vars_presence_absence_t)
#dim(all_vars_presence_absence_t_subset)

#split the cols by gene cluster 
cluster_genes_split<-split(SMC_genes$Af293_gene_name, SMC_genes$Cluster_number)

subdf1<- all_vars_presence_absence_t[,colnames(all_vars_presence_absence_t) %in% cluster_genes_split$`BGC 1`]
subdf2<- all_vars_presence_absence_t[,colnames(all_vars_presence_absence_t) %in% cluster_genes_split$`BGC 2`]
subdf3<- all_vars_presence_absence_t[,colnames(all_vars_presence_absence_t) %in% cluster_genes_split$`BGC 3`]
subdf4<- all_vars_presence_absence_t[,colnames(all_vars_presence_absence_t) %in% cluster_genes_split$`BGC 4`]
subdf5<- all_vars_presence_absence_t[,colnames(all_vars_presence_absence_t) %in% cluster_genes_split$`BGC 5`]
subdf6<- all_vars_presence_absence_t[,colnames(all_vars_presence_absence_t) %in% cluster_genes_split$`BGC 6`]
subdf7<- all_vars_presence_absence_t[,colnames(all_vars_presence_absence_t) %in% cluster_genes_split$`BGC 7`]
subdf8<- all_vars_presence_absence_t[,colnames(all_vars_presence_absence_t) %in% cluster_genes_split$`BGC 8`]
subdf9<- all_vars_presence_absence_t[,colnames(all_vars_presence_absence_t) %in% cluster_genes_split$`BGC 9`]
subdf10<- all_vars_presence_absence_t[,colnames(all_vars_presence_absence_t) %in% cluster_genes_split$`BGC 10`]
subdf11<- all_vars_presence_absence_t[,colnames(all_vars_presence_absence_t) %in% cluster_genes_split$`BGC 11`]
subdf12<- all_vars_presence_absence_t[,colnames(all_vars_presence_absence_t) %in% cluster_genes_split$`BGC 12`]
subdf13<- all_vars_presence_absence_t[,colnames(all_vars_presence_absence_t) %in% cluster_genes_split$`BGC 13`]
subdf14<- all_vars_presence_absence_t[,colnames(all_vars_presence_absence_t) %in% cluster_genes_split$`BGC 14`]
subdf15<- all_vars_presence_absence_t[,colnames(all_vars_presence_absence_t) %in% cluster_genes_split$`BGC 15`]
subdf16<- all_vars_presence_absence_t[,colnames(all_vars_presence_absence_t) %in% cluster_genes_split$`BGC 16`]
subdf17<- all_vars_presence_absence_t[,colnames(all_vars_presence_absence_t) %in% cluster_genes_split$`BGC 17`]
subdf18<- all_vars_presence_absence_t[,colnames(all_vars_presence_absence_t) %in% cluster_genes_split$`BGC 18`]
subdf19<- all_vars_presence_absence_t[,colnames(all_vars_presence_absence_t) %in% cluster_genes_split$`BGC 19`]
subdf20<- all_vars_presence_absence_t[,colnames(all_vars_presence_absence_t) %in% cluster_genes_split$`BGC 20`]
subdf21<- all_vars_presence_absence_t[,colnames(all_vars_presence_absence_t) %in% cluster_genes_split$`BGC 21`]
subdf22<- all_vars_presence_absence_t[,colnames(all_vars_presence_absence_t) %in% cluster_genes_split$`BGC 22`]
subdf23<- all_vars_presence_absence_t[,colnames(all_vars_presence_absence_t) %in% cluster_genes_split$`BGC 23`]
subdf24<- all_vars_presence_absence_t[,colnames(all_vars_presence_absence_t) %in% cluster_genes_split$`BGC 24`]
subdf25<- all_vars_presence_absence_t[,colnames(all_vars_presence_absence_t) %in% cluster_genes_split$`Super BGC 25`]
subdf26<- all_vars_presence_absence_t[,colnames(all_vars_presence_absence_t) %in% cluster_genes_split$`BGC 26`]

#attach cluster number to each df to keep track 
colnames(subdf1)<-paste(colnames(subdf1),"BGC_1_(Uncharacterized polyketide)",sep="_")
colnames(subdf2)<-paste(colnames(subdf2),"BGC_2_(Nidulanin)",sep="_")
colnames(subdf3)<-paste(colnames(subdf3),"BGC_3_(Ferricrocin)",sep="_")
colnames(subdf4)<-paste(colnames(subdf4),"BGC_4_(Fusarielin-like polyketide)",sep="_")
colnames(subdf5)<-paste(colnames(subdf5),"BGC_5_(DHN Melanin)",sep="_")
colnames(subdf6)<-paste(colnames(subdf6),"BGC_6_(Fumigaclavine)",sep="_")
colnames(subdf7)<-paste(colnames(subdf7),"BGC_7_(Uncharacterized polyketide)",sep="_")
colnames(subdf8)<-paste(colnames(subdf8),"BGC_8_(Fusarine C)",sep="_")
colnames(subdf9)<-paste(colnames(subdf9),"BGC_9_(Hexadehydroastechrome)",sep="_")
colnames(subdf10)<-paste(colnames(subdf10),"BGC_10_(Uncharacterized non-ribosomal peptide)",sep="_")
colnames(subdf11)<-paste(colnames(subdf11),"BGC_11_(Uncharacterized polyketide)",sep="_")
colnames(subdf12)<-paste(colnames(subdf12),"BGC_12_(Uncharacterized non-ribosomal peptide)",sep="_")
colnames(subdf13)<-paste(colnames(subdf13),"BGC_13_(Endocrocin)",sep="_")
colnames(subdf14)<-paste(colnames(subdf14),"BGC_14_(Trypacidin)",sep="_")
colnames(subdf15)<-paste(colnames(subdf15),"BGC_15_(Protostadienol, helvolic acid)",sep="_")
colnames(subdf16)<-paste(colnames(subdf16),"BGC_16_(Uncharacterized non-ribosomal peptide-like)",sep="_")
colnames(subdf17)<-paste(colnames(subdf17),"BGC_17_(uncharacterized non-ribosomal peptide)",sep="_")
colnames(subdf18)<-paste(colnames(subdf18),"BGC_18_(Fumisoquin)",sep="_")
colnames(subdf19)<-paste(colnames(subdf19),"BGC_19_(Uncharacterized non-ribosomal peptide)",sep="_")
colnames(subdf20)<-paste(colnames(subdf20),"BGC_20_(Gliotoxin)",sep="_")
colnames(subdf21)<-paste(colnames(subdf21),"BGC_21_(Fumiquinozalines)",sep="_")
colnames(subdf22)<-paste(colnames(subdf22),"BGC_22_(Pyripyropene A)",sep="_")
colnames(subdf23)<-paste(colnames(subdf23),"BGC_23_(Neosartoricin)",sep="_")
colnames(subdf24)<-paste(colnames(subdf24),"BGC_24_(Fumitremorgin)",sep="_")
colnames(subdf25)<-paste(colnames(subdf25),"BGC_25_(Intertwined Fumagillin and pseurotin)",sep="_")
colnames(subdf26)<-paste(colnames(subdf26),"BGC_26_(Uncharacterized pks and terpene hybrid)",sep="_")

#bind all into one DF
all_vars_df<- subdf1
for (i in 2:26) {
  all_vars_df <- cbind(all_vars_df, eval(parse(text=paste("subdf", i, sep = ''))))
}


#check that the names match
setdiff(tree$tip.label,row.names(all_vars_df)) #match


#plotall_tree_plot <-  gheatmap(tree_plot_me, 
#plot
all_tree_plot <-  gheatmap(tree_plot_me, 
                           all_vars_df,
                           offset=.1, width=1, low="white", high="black", 
                           colnames = T, 
                           colnames_angle = 90,
                           colnames_position = "top",
                           font.size = .8,
                           colnames_offset_y = 18,
                           #colnames_offset_x = -.1,
                           color="white") +
  scale_fill_manual(values=c("white", "black")) 
all_tree_plot
all_tree<- all_tree_plot + geom_tiplab(size = .8, align = TRUE, linesize = .25, offset = 1.4, linetype = 0) + ggtree::vexpand(.15, 1)
all_tree
ggsave(file="all_SMCs_2.pdf",device="pdf", all_tree, width=31, height=8, units="in")


##find missing genes (are these clustered with other families? or unannotated?)
#find missing genes in OFtoAfu
missing_genes

missing_genes_all_results<- OFtoAfu[,OFtoAfu %in% missing_genes]

missing_genes_all_results<-OFtoAfu[rowSums(sapply(missing_genes, grepl, OFtoAfu$qseqid)) > 0, , drop = FALSE]


#genes with a match 
genes_that_match<- intersect(missing_genes_all_results$qseqid, missing_genes)
genes_that_match

#View(missing_genes_all_results)
#"Afu3g13700" is in BGC 10 (OG0005716)
#"Afu4g00220" is in BCG 13 (OG0003471)
#"Afu4g14560" is in BGC 14 (OG0006238)
#"Afu6g09730" is in BGC 20 (OG0007452)
#"Afu6g13980" is in BGC 22 (OG0007936)
#"Afu8g00490" is in BGC 25 (OG0000692)
#"Afu8g02380" is in BGC 26 (OG0004431)


#what genes are they convoluted with? 
OG0005716<-OFtoAfu[OFtoAfu$sseqid == "OG0005716",]
OG0005716 # Afu3g13700 is convoluted w.  Afu3g13690 (BGC 10)

OG0003471<-OFtoAfu[OFtoAfu$sseqid == "OG0003471",]
OG0003471 # Afu4g00220 is convoluted w.  Afu4g00210 (BGC 13)

OG0006238<-OFtoAfu[OFtoAfu$sseqid == "OG0006238",]
OG0006238 # Afu4g14560 is convoluted w.  Afu4g14550 (BGC 14)

OG0007452<-OFtoAfu[OFtoAfu$sseqid == "OG0007452",]
OG0007452 # Afu6g09730 is convoluted w.  Afu6g09720 (BGC 20)

OG0007936<-OFtoAfu[OFtoAfu$sseqid == "OG0007936",]
OG0007936 # Afu6g13980 is convoluted w.  Afu6g13970 (BGC 22)

OG0000692<-OFtoAfu[OFtoAfu$sseqid == "OG0000692",]
OG0000692 # Afu8g00490 is convoluted w.  Afu8g00480 (BGC 25)

OG0004431<-OFtoAfu[OFtoAfu$sseqid == "OG0004431",]
OG0004431 # Afu8g02380 is convoluted w.  Afu8g02360 (BGC 26)


#genes that don't match
genes_no_match<- setdiff(missing_genes, missing_genes_all_results$qseqid)
genes_no_match

#View(SMC_genes)
#"Afu1g10275" is in BGC 2
#"Afu7g00140" is in BGC 23
#"Afu8g00450" is in BGC 25

#10 out of how many genes? 
length(unique(SMC_genes$Af293_gene_name)) #230


#what are the most conserved clusters? the most variable clusters?
#write function to get frequency
absent_in_x_strains <- function (df) {
df_presence<- data.frame(sapply(df, gsub, pattern = "present", replacement = 1))
df_presence_absence<- data.frame(sapply(df_presence, gsub, pattern = "absent", replacement = 0))
df_presence_absence[] <- lapply(df_presence_absence, function(x) {
  if(is.character(x)) as.numeric(as.character(x)) else x
})
df_presence_absence[] <- lapply(df_presence_absence, as.numeric)
df_presence_absence$sum <- rowSums(df_presence_absence, na.rm = TRUE)
n_absent<- sum(df_presence_absence$sum < ncol(df_presence_absence)-1)
n_absent
}

C1<- absent_in_x_strains(df = subdf1)
C2<- absent_in_x_strains(df = subdf2)
C3<- absent_in_x_strains(df = subdf3)
C4<- absent_in_x_strains(df = subdf4)
C5<- absent_in_x_strains(df = subdf5)
C6<- absent_in_x_strains(df = subdf6)
C7<- absent_in_x_strains(df = subdf7)
C8<- absent_in_x_strains(df = subdf8)
C9<- absent_in_x_strains(df = subdf9)
C10<- absent_in_x_strains(df = subdf10)
C11<- absent_in_x_strains(df = subdf11)
C12<- absent_in_x_strains(df = subdf12)
C13<- absent_in_x_strains(df = subdf13)
C14<- absent_in_x_strains(df = subdf14)
C15<- absent_in_x_strains(df = subdf15)
C16<- absent_in_x_strains(df = subdf16)
C17<- absent_in_x_strains(df = subdf17)
C18<- absent_in_x_strains(df = subdf18)
C19<- absent_in_x_strains(df = subdf9)
C20<- absent_in_x_strains(df = subdf20)
C21<- absent_in_x_strains(df = subdf21)
C22<- absent_in_x_strains(df = subdf22)
C23<- absent_in_x_strains(df = subdf23)
C24<- absent_in_x_strains(df = subdf24)
C25<- absent_in_x_strains(df = subdf25)
C26<- absent_in_x_strains(df = subdf26)

totals<- c(C1,C2, C3, C4, C5, C6, C7, C8, C9, C10, C11, C12, C13, C14, C15,
           C16, C17, C18, C19, C20, C21, C22, C23, C24, C25, C26)
totals_df<- data.frame(cbind(Cluster= paste("Cluster",1:26,sep="_"), n_strains_w_genes_absent = as.numeric(totals)))

totals_df[mixedorder(totals_df$n_strains_w_genes_absent),]



