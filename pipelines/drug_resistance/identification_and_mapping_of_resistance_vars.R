##This is a test to look at the distribution of known resistance variants in the pop-genome data for 260 strains of A.fumigatus
#last updated 30.Jan.2021

#load modules
require(data.table)
require(splitstackshape)
library(plyr)
library(dplyr)
library(tidyverse)
library(rlist)
library(gdata)
library(ggtree)
library(ape)
library(phytools)
library(ggplot2)
library(stringr)

#set wd
setwd("~/Desktop/Project_Afum_pangenome_3/")
options(stringsAsFactors = FALSE)

#sessionInfo() 

##set seed for reproducibility
set.seed(666)


##read in data
snpEff_all<-read.delim("Pop_for_pan_260.v2.snpEff.matrix.tsv", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE, check.names = FALSE)
Afum_grp<-read.delim("clade_map_K3_3Jan2022.txt", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE, check.names = TRUE)
resistance_db<-read.delim("Afum_azole_mutations.csv", header = TRUE, sep = ",", fill = TRUE, strip.white = TRUE)
tree <- read.tree("Afum_260_iq_tree_newick_v2.tre")

#remove the reference from the tree:
tree<- drop.tip(tree, "Af293-REF", trim.internal = TRUE, subtree = FALSE,
                   root.edge = 0, rooted = is.rooted(tree), collapse.singles = TRUE,
                   interactive = FALSE)
#root tree
tree<- root(tree, node = 505, resolve.root = TRUE,
               interactive = FALSE, edgelabel = FALSE)


###clean the SNP data for easy processing###
#remove the brackets form ALT
snpEff_all$ALT<-gsub("\\[|\\]", "", snpEff_all$ALT)

#remove the "/" from all the the SNP calls 
snpEff_all[] <- lapply(snpEff_all, gsub, pattern='/', replacement="")

#get set size
dim(snpEff_all[,11:(ncol(snpEff_all) -1)])
#260 strains, 441,000 variants 

#subset to only missense variants 
snpEff_no_intergenic<- snpEff_all[snpEff_all$TYPE == "missense_variant",]
dim(snpEff_no_intergenic)
#75,730 missence variants

###sep each gene of interest in the resistance_db, sep snpEff_no_intergenic into their own dfs
#set genes of interest
GOI<- unique(resistance_db$A_fum_gene_name)
length(GOI)

#subset df into just genes of interest
snpEff_no_intergenic_GOI <- snpEff_no_intergenic[snpEff_no_intergenic$GENE %in% GOI, ]
#how many genes have variants in them?
genes_w_vars<- unique(snpEff_no_intergenic_GOI$GENE)
genes_w_vars
length(genes_w_vars)
#variants in all 13 genes 

#split each gene into it's own df (keep in a list)
listDf_GOI_variants <- split(snpEff_no_intergenic_GOI, f = snpEff_no_intergenic_GOI$GENE)

##split resistance_db by gene
#subset to only genes with variants in the dataset
resistance_db_GOI <- resistance_db[resistance_db$A_fum_gene_name %in% genes_w_vars, ]
#check
unique(resistance_db$A_fum_gene_name)
unique(resistance_db_GOI$A_fum_gene_name) #works
#split each database entry by gene into it's own df (keep in a list)
listDf_resistance_db_GOI <- split(resistance_db_GOI, f = resistance_db_GOI$A_fum_gene_name)


#function to subset all variants in genes of interest to only known resistance mutations
get_resistance_mutations <- function(variants_OI, resistance_db_OI){
  output<- list()
  #for each gene, subset the rows if variants are known resistance mutations 
  for (i in 1:length(variants_OI)){
      cyp51_db<- resistance_db_OI[[i]]
      cyp51_var<- variants_OI[[i]]
      output[[i]] <- cyp51_var[grep(paste(cyp51_db$mutation_code2,collapse="|"), cyp51_var$CHANGEPEP), ]
    }
  names(output)<- names(resistance_db_OI)
  no_empties<- output[sapply(output, function(x) nrow(x)) > 0]
  return(no_empties)
}

#alternate function to find variable mutation at the same position as known resistance variants
get_resistance_mutations_variable <- function(variants_OI, resistance_db_OI){
  output<- list()
  #for each gene, subset the rows if variants are known resistance mutations 
  for (i in 1:length(variants_OI)){
    cyp51_db<- resistance_db_OI[[i]]
    cyp51_var<- variants_OI[[i]]
    output[[i]] <- cyp51_var[grep(paste(cyp51_db$mutation_code2_start_only,collapse="|"), cyp51_var$CHANGEPEP), ]
  }
  names(output)<- names(resistance_db_OI)
  no_empties<- output[sapply(output, function(x) nrow(x)) > 0]
  return(no_empties)
}

#run function
list_of_dfs_resistance_vars_in_dataset<- get_resistance_mutations(variants_OI = listDf_GOI_variants, resistance_db_OI = listDf_resistance_db_GOI)
names(list_of_dfs_resistance_vars_in_dataset)
#charicterized resistance variants in Cyp51A only

#isolate Cyp51A
Cyp51a<- list_of_dfs_resistance_vars_in_dataset[[1]]
#View(Cyp51a)

#How about if it's a variable change at the same position?
list_of_dfs_resistance_vars_in_dataset_variable<- get_resistance_mutations_variable(variants_OI = listDf_GOI_variants, resistance_db_OI = listDf_resistance_db_GOI)
names(list_of_dfs_resistance_vars_in_dataset_variable)
# - now there's also Afu2g03700
#how many have this mutation? 
Afu2g03700<- list_of_dfs_resistance_vars_in_dataset_variable$Afu2g03700
cols_to_exclude<- c("CHROM",
                    "POS",
                    "FLANKING",
                    "TYPE",
                    "IMPACT",
                    "GENE",
                    "CHANGEDNA",
                    "CHANGEPEP",
                    "REF",
                    "ALT",
                    "ANN")
                    
                    

Afu2g03700_only<- Afu2g03700[,!names(Afu2g03700) %in% cols_to_exclude]
rowSums(Afu2g03700_only == "C")#three strains have this variant


#isolate variable (uncharacterized) Cyp51A variants at known positions
Cyp51a_variable<- list_of_dfs_resistance_vars_in_dataset_variable[[2]]

#check if they're the same
Cyp51a$CHANGEPEP == Cyp51a_variable$CHANGEPEP
#yes - no new aa changes at positions of known resistance variants.


##make presence/absence dataframe for presence of variants in each gene and map onto tree file. 

#Function to get the number of variants in ea gene for each strain
#over all: for each strain, for each gene, are the calls real variants (they don't match the REF?) If so, TRUE, if no FALSE. 
n_variants_per_gene_per_strain<- function(gene_df){
  #shrink the df 
  REF<- data.frame(gene_df$REF)
  strains<- gene_df[, 11:(ncol(gene_df) -1)]
  #create true/false data frame: true if not == reference, and not == "." (possible misscalls)
  T_F_df<- data.frame(apply(strains, 2, function(x) (x != REF$gene_df.REF) & (x != ".")))
  if (nrow(strains) >1) {
    T_F_df<- data.frame(apply(strains, 2, function(x) (x != REF$gene_df.REF) & (x != ".")))
  } else {
    T_F_df<- t(data.frame(apply(strains, 2, function(x) (x != REF$gene_df.REF) & (x != "."))))
  }
  #get col. sums (the number of TRUE values == the number of variants in that gene for that strain)
  n_vars_in_gene<- data.frame(colSums(T_F_df))
  #rename for later
  df.name<- deparse(substitute(gene_df))
  colnames(n_vars_in_gene)<- df.name 
  return(n_vars_in_gene)
}


#use the above function in it's own function to loop over each mutation to make a df for graphing and get totals across all strains.
make_by_variant_output <- function(input){
  output<- data.frame(matrix(NA, nrow = ncol(input[, 11:(ncol(input) -1)]), ncol = 1))
  #for each row (variant) get presence/ absence for which strains have the variant
  for (i in 1:nrow(input)){
    output[[i]]<- n_variants_per_gene_per_strain(input[i,])
    this_name<- input[i,8]
    names(output)[[i]]<- this_name
  }
  return(output)
}

#run function
by_variant_output_df<- make_by_variant_output(input = Cyp51a)
colSums(by_variant_output_df)

Cyp51a

#Gly448Ser Met220Val Pro216Leu His147Tyr Gly138Cys  Leu98His  Gly54Glu 
#1         1         4         1         7        19         2 


#do the same for all non-syn changes in Cyp51A (not just known resistance aleles)
Cyp51a_all<- snpEff_no_intergenic_GOI[snpEff_no_intergenic_GOI$GENE == "Afu4g06890",]

#remove the known virulence positions identified above 
Cyp51a_all_but_vir<- Cyp51a_all[!Cyp51a_all$CHANGEPEP %in% names(by_variant_output_df),]

#make presence/absence df by strain
by_variant_output_df_all_but_vir<- make_by_variant_output(input = Cyp51a_all_but_vir)
colSums(by_variant_output_df_all_but_vir)

#write.table(by_variant_output_df, "variants_in_resistance_aleles.csv", sep=",", row.names = TRUE, col.names = TRUE, quote = FALSE)
#write.table(by_variant_output_df_all_but_vir, "all_nonsyn_variants_in_cyp51a.csv", sep=",", row.names = TRUE, col.names = TRUE, quote = FALSE)

colSums(by_variant_output_df)
##map onto tree

#set colors by group for clade
grA_me<- split(Afum_grp$name_pop_genome_new, Afum_grp$clade)
tree_grA_me <- ggtree::groupOTU(tree, grA_me)
str(tree_grA_me)
levels(attributes(tree_grA_me)$group) 
levels(attributes(tree_grA_me)$group)[1] <- "1"
# Reorder factor levels if needed (this controls the order in the legend)
#attributes(tree_grA_me)$group <- factor(x = attributes(tree_grA_me)$group, 
#                                        levels = c("clade1", "clade2", "clade3"))

attributes(tree_grA_me)$group <- factor(x = attributes(tree_grA_me)$group, 
                                        levels = c("1", "2", "3"))


#set colors
my_cols_me <- c(clade1 = "#56326E",
                clade2 ="#ED7F6F",
                clade3 = "#ABA778")

names(my_cols_me) <- levels(attributes(tree_grA_me)$group)
scales::show_col(my_cols_me); my_cols_me


#format resistance allele data for graphing
#colapse df of df because you're an idiot and didn't append these as a list in the first place
by_variant_output_df_colapse<- bind_rows(lapply(by_variant_output_df,function(i)do.call(cbind,i)))
row.names(by_variant_output_df_colapse)<- row.names(by_variant_output_df[[1]])

#same for unknown function vars
by_variant_output_df_colapse_unknown<- bind_rows(lapply(by_variant_output_df_all_but_vir,function(i)do.call(cbind,i)))
row.names(by_variant_output_df_colapse_unknown)<- row.names(by_variant_output_df_all_but_vir[[1]])


#order columns by abundance
colSums(by_variant_output_df_colapse)
by_abundance<- by_variant_output_df_colapse[,order(colSums(-by_variant_output_df_colapse))]
colSums(by_abundance) #works
#get total genomes influenced 
list<- colSums(by_abundance)
sum(list)
#same for unknowns
colSums(by_variant_output_df_colapse_unknown)
by_abundance_unknown<- by_variant_output_df_colapse_unknown[,order(colSums(-by_variant_output_df_colapse_unknown))]
colSums(by_abundance_unknown) #works

#set presence/absence
by_variant_output_df_binary_presence<- data.frame(sapply(by_abundance, gsub, pattern = "1", replacement = "present"))
by_variant_output_df_binary_presence_absence<- data.frame(sapply(by_variant_output_df_binary_presence, gsub, pattern = "0", replacement = "absent"))
row.names(by_variant_output_df_binary_presence_absence)<- row.names(by_variant_output_df[[1]])

#set presence/absence for unknowns
by_variant_output_df_binary_presence_unknowns<- data.frame(sapply(by_abundance_unknown, gsub, pattern = "1", replacement = "present"))
by_variant_output_df_binary_presence_absence_unknowns<- data.frame(sapply(by_variant_output_df_binary_presence_unknowns, gsub, pattern = "0", replacement = "absent"))
row.names(by_variant_output_df_binary_presence_absence_unknowns)<- row.names(by_variant_output_df[[1]])


#check that all the names match and you're not missing any data
#tree_node_names<- tree_grA_me$tip.label
#df_names<- rownames(by_variant_output_df_binary_presence_absence)
#setdiff(tree_node_names, df_names)
#none missing


#plot tree
tree_plot_me <- 
  ggtree(tr = tree_grA_me, 
         mapping = aes(color = group), 
         layout  = 'rectangular') + 
         #ladderize = TRUE) +
         #branch.length = 'none', 
         geom_treescale(x=5, y=NULL, color = NA) +
         # set line thickness
         #size = .3) +
  # adjust coloring of main groups
  scale_color_manual(name = 'Clade', values = my_cols_me) + 
  theme(legend.title=element_text(size=9), # The title of legend 
        legend.text=element_text(size=7))+
  guides(color = guide_legend(override.aes = list(size = 4)))

#asthetics
tree_plot_me<- tree_plot_me +geom_tiplab(size = 0, align = TRUE, linesize = .25, linetype = 3)
tree_plot_me

ncol(by_variant_output_df_binary_presence_absence)
#plot resistance variants onto tree
#plot
#bind and highlight - can't get width right when adjacent. 
all_vars<- cbind(by_variant_output_df_binary_presence_absence, by_variant_output_df_binary_presence_absence_unknowns)
cyp51A_tree_plot <-  gheatmap(tree_plot_me, 
                              #by_variant_output_df_binary_presence_absence_unknowns,
                              all_vars,
                              offset=0.02, width=0.7, low="white", high="black", 
                              colnames = T, 
                              colnames_angle = 45,
                              colnames_position = "top",
                              font.size = 1,
                              colnames_offset_y = -6,
                              colnames_offset_x = .05,
                              color="white") +
  scale_fill_manual(values=c("white", "black")) 
#+ggtitle("Cyp51 variants")

p<- cyp51A_tree_plot + geom_tiplab(size = .8, align = TRUE, linesize = .25, offset = 1, linetype = 0)
p
#export
#ggsave(file="Cyp51A_variants.pdf",device="pdf", p, width=8, height=8, units="in")



#make supplemental figure of all missence variants from all strains (including those for which only expression data is available)
genes_w_vars

#isolate the four other genes with non-syn changes in resistance genes:
fks1_all<- snpEff_no_intergenic_GOI[snpEff_no_intergenic_GOI$GENE == "Afu6g12400",]
nrow(fks1_all)# 5 positions
atrF_all<- snpEff_no_intergenic_GOI[snpEff_no_intergenic_GOI$GENE == "Afu6g04360",]
nrow(atrF_all)# 12 positions
mdr1_all<- snpEff_no_intergenic_GOI[snpEff_no_intergenic_GOI$GENE == "Afu5g06070",]
nrow(mdr1_all)# 9 positions
mdr3_all<- snpEff_no_intergenic_GOI[snpEff_no_intergenic_GOI$GENE == "Afu3g03500",]
nrow(mdr3_all)# 15 positions
hapE_all<- snpEff_no_intergenic_GOI[snpEff_no_intergenic_GOI$GENE == "Afu2g14720",]
nrow(hapE_all)# 1 positions
hmg1_all<- snpEff_no_intergenic_GOI[snpEff_no_intergenic_GOI$GENE == "Afu2g03700",]
nrow(hmg1_all)# 10 positions
cdr1B_abcC_all<- snpEff_no_intergenic_GOI[snpEff_no_intergenic_GOI$GENE == "Afu1g14330",]
nrow(cdr1B_abcC_all)# 8 positions
mdr4_all<- snpEff_no_intergenic_GOI[snpEff_no_intergenic_GOI$GENE == "Afu1g12690",]
nrow(mdr4_all)# 15 positions
cox10_all<- snpEff_no_intergenic_GOI[snpEff_no_intergenic_GOI$GENE == "Afu4g08340",]
nrow(cox10_all)# 12 positions
mdr2_all<- snpEff_no_intergenic_GOI[snpEff_no_intergenic_GOI$GENE == "Afu4g10000",]
nrow(mdr2_all)# 6 positions
cyp51B_all<- snpEff_no_intergenic_GOI[snpEff_no_intergenic_GOI$GENE == "Afu7g03740",]
nrow(cyp51B_all)# 2 positions
AFUA_7G01960_all<- snpEff_no_intergenic_GOI[snpEff_no_intergenic_GOI$GENE == "Afu7g01960",]
nrow(AFUA_7G01960_all)# 7 positions


#make presence/absence df by strain
fks1_all_PA<- make_by_variant_output(input = fks1_all)
colSums(fks1_all_PA)
atrF_all_PA<- make_by_variant_output(input = atrF_all)
colSums(atrF_all_PA)
mdr1_all_PA<- make_by_variant_output(input = mdr1_all)
colSums(mdr1_all_PA)
mdr3_all_PA<- make_by_variant_output(input = mdr3_all)
colSums(mdr3_all_PA)
hapE_all_PA<- make_by_variant_output(input = hapE_all)
colSums(hapE_all_PA)
hmg1_all_PA<- make_by_variant_output(input = hmg1_all)
colSums(hmg1_all_PA)
cdr1B_abcC_all_PA<- make_by_variant_output(input = cdr1B_abcC_all)
colSums(cdr1B_abcC_all_PA)
mdr4_all_PA<- make_by_variant_output(input = mdr4_all)
colSums(mdr4_all_PA)
cox10_all_PA<- make_by_variant_output(input = cox10_all)
colSums(cox10_all_PA)
mdr2_all_PA<- make_by_variant_output(input = mdr2_all)
colSums(mdr2_all_PA)
cyp51B_all_PA<- make_by_variant_output(input = cyp51B_all)
colSums(cyp51B_all_PA)
AFUA_7G01960_all_PA<- make_by_variant_output(input = AFUA_7G01960_all)
colSums(AFUA_7G01960_all_PA)

##format resistance allele data for graphing
#colapse df of df because you're an idiot and didn't append these as a list in the first place
fks1_colapse<- bind_rows(lapply(fks1_all_PA,function(i)do.call(cbind,i)))
row.names(fks1_colapse)<- row.names(fks1_all_PA[[1]])
atrF_colapse<- bind_rows(lapply(atrF_all_PA,function(i)do.call(cbind,i)))
row.names(atrF_colapse)<- row.names(atrF_all_PA[[1]])
mdr1_colapse<- bind_rows(lapply(mdr1_all_PA,function(i)do.call(cbind,i)))
row.names(mdr1_colapse)<- row.names(mdr1_all_PA[[1]])
mdr3_colapse<- bind_rows(lapply(mdr3_all_PA,function(i)do.call(cbind,i)))
row.names(mdr3_colapse)<- row.names(mdr3_all_PA[[1]])
hapE_colapse<- bind_rows(lapply(hapE_all_PA,function(i)do.call(cbind,i)))
row.names(hapE_colapse)<- row.names(hapE_all_PA[[1]])
hmg1_colapse<- bind_rows(lapply(hmg1_all_PA,function(i)do.call(cbind,i)))
row.names(hmg1_colapse)<- row.names(hmg1_all_PA[[1]])
cdr1B_abcC_colapse<- bind_rows(lapply(cdr1B_abcC_all_PA,function(i)do.call(cbind,i)))
row.names(cdr1B_abcC_colapse)<- row.names(cdr1B_abcC_all_PA[[1]])
mdr4_colapse<- bind_rows(lapply(mdr4_all_PA,function(i)do.call(cbind,i)))
row.names(mdr4_colapse)<- row.names(mdr4_all_PA[[1]])
cox10_colapse<- bind_rows(lapply(cox10_all_PA,function(i)do.call(cbind,i)))
row.names(cox10_colapse)<- row.names(cox10_all_PA[[1]])
mdr2_colapse<- bind_rows(lapply(mdr2_all_PA,function(i)do.call(cbind,i)))
row.names(mdr2_colapse)<- row.names(mdr2_all_PA[[1]])
cyp51B_colapse<- bind_rows(lapply(cyp51B_all_PA,function(i)do.call(cbind,i)))
row.names(cyp51B_colapse)<- row.names(cyp51B_all_PA[[1]])
AFUA_7G01960_colapse<- bind_rows(lapply(AFUA_7G01960_all_PA,function(i)do.call(cbind,i)))
row.names(AFUA_7G01960_colapse)<- row.names(AFUA_7G01960_all_PA[[1]])


#order columns by abundance
fks1_by_abundance<- fks1_colapse[,order(colSums(-fks1_colapse))]
atrF_by_abundance<- atrF_colapse[,order(colSums(-atrF_colapse))]
mdr1_by_abundance<- mdr1_colapse[,order(colSums(-mdr1_colapse))]
mdr3_by_abundance<- mdr3_colapse[,order(colSums(-mdr3_colapse))]
hapE_by_abundance<- hapE_colapse[,order(colSums(-hapE_colapse))]
hmg1_by_abundance<- hmg1_colapse[,order(colSums(-hmg1_colapse))]
cdr1B_abcC_by_abundance<- cdr1B_abcC_colapse[,order(colSums(-cdr1B_abcC_colapse))]
mdr4_by_abundance<- mdr4_colapse[,order(colSums(-mdr4_colapse))]
cox10_by_abundance<- cox10_colapse[,order(colSums(-cox10_colapse))]
mdr2_by_abundance<- mdr2_colapse[,order(colSums(-mdr2_colapse))]
cyp51B_by_abundance<- cyp51B_colapse[,order(colSums(-cyp51B_colapse))]
AFUA_7G01960_by_abundance<- AFUA_7G01960_colapse[,order(colSums(-AFUA_7G01960_colapse))]


#set presence/absence
fks1_by_abundance_presence<- data.frame(sapply(fks1_by_abundance, gsub, pattern = "1", replacement = "present"))
fks1_by_abundance_presence_absence<- data.frame(sapply(fks1_by_abundance_presence, gsub, pattern = "0", replacement = "absent"))
row.names(fks1_by_abundance_presence_absence)<- row.names(fks1_colapse)

atrF_by_abundance_presence<- data.frame(sapply(atrF_by_abundance, gsub, pattern = "1", replacement = "present"))
atrF_by_abundance_presence_absence<- data.frame(sapply(atrF_by_abundance_presence, gsub, pattern = "0", replacement = "absent"))
row.names(atrF_by_abundance_presence_absence)<- row.names(atrF_colapse)

mdr1_by_abundance_presence<- data.frame(sapply(mdr1_by_abundance, gsub, pattern = "1", replacement = "present"))
mdr1_by_abundance_presence_absence<- data.frame(sapply(mdr1_by_abundance_presence, gsub, pattern = "0", replacement = "absent"))
row.names(mdr1_by_abundance_presence_absence)<- row.names(mdr1_colapse)

mdr3_by_abundance_presence<- data.frame(sapply(mdr3_by_abundance, gsub, pattern = "1", replacement = "present"))
mdr3_by_abundance_presence_absence<- data.frame(sapply(mdr3_by_abundance_presence, gsub, pattern = "0", replacement = "absent"))
row.names(mdr3_by_abundance_presence_absence)<- row.names(mdr3_colapse)

hapE_by_abundance_presence<- data.frame(sapply(hapE_by_abundance, gsub, pattern = "1", replacement = "present"))
hapE_by_abundance_presence_absence<- data.frame(sapply(hapE_by_abundance_presence, gsub, pattern = "0", replacement = "absent"))
row.names(hapE_by_abundance_presence_absence)<- row.names(hapE_colapse)

hmg1_by_abundance_presence<- data.frame(sapply(hmg1_by_abundance, gsub, pattern = "1", replacement = "present"))
hmg1_by_abundance_presence_absence<- data.frame(sapply(hmg1_by_abundance_presence, gsub, pattern = "0", replacement = "absent"))
row.names(hmg1_by_abundance_presence_absence)<- row.names(hmg1_colapse)

cdr1B_abcC_by_abundance_presence<- data.frame(sapply(cdr1B_abcC_by_abundance, gsub, pattern = "1", replacement = "present"))
cdr1B_abcC_by_abundance_presence_absence<- data.frame(sapply(cdr1B_abcC_by_abundance_presence, gsub, pattern = "0", replacement = "absent"))
row.names(cdr1B_abcC_by_abundance_presence_absence)<- row.names(cdr1B_abcC_colapse)

mdr4_by_abundance_presence<- data.frame(sapply(mdr4_by_abundance, gsub, pattern = "1", replacement = "present"))
mdr4_by_abundance_presence_absence<- data.frame(sapply(mdr4_by_abundance_presence, gsub, pattern = "0", replacement = "absent"))
row.names(mdr4_by_abundance_presence_absence)<- row.names(mdr4_colapse)

cox10_by_abundance_presence<- data.frame(sapply(cox10_by_abundance, gsub, pattern = "1", replacement = "present"))
cox10_by_abundance_presence_absence<- data.frame(sapply(cox10_by_abundance_presence, gsub, pattern = "0", replacement = "absent"))
row.names(cox10_by_abundance_presence_absence)<- row.names(cox10_colapse)

mdr2_by_abundance_presence<- data.frame(sapply(mdr2_by_abundance, gsub, pattern = "1", replacement = "present"))
mdr2_by_abundance_presence_absence<- data.frame(sapply(mdr2_by_abundance_presence, gsub, pattern = "0", replacement = "absent"))
row.names(mdr2_by_abundance_presence_absence)<- row.names(mdr2_colapse)

cyp51B_by_abundance_presence<- data.frame(sapply(cyp51B_by_abundance, gsub, pattern = "1", replacement = "present"))
cyp51B_by_abundance_presence_absence<- data.frame(sapply(cyp51B_by_abundance_presence, gsub, pattern = "0", replacement = "absent"))
row.names(cyp51B_by_abundance_presence_absence)<- row.names(cyp51B_colapse)

AFUA_7G01960_by_abundance_presence<- data.frame(sapply(AFUA_7G01960_by_abundance, gsub, pattern = "1", replacement = "present"))
AFUA_7G01960_by_abundance_presence_absence<- data.frame(sapply(AFUA_7G01960_by_abundance_presence, gsub, pattern = "0", replacement = "absent"))
row.names(AFUA_7G01960_by_abundance_presence_absence)<- row.names(AFUA_7G01960_colapse)


#add suffix to collumn names
colnames(fks1_by_abundance_presence_absence)<-paste(colnames(fks1_by_abundance_presence_absence),"fks1",sep="_")
colnames(atrF_by_abundance_presence_absence)<-paste(colnames(atrF_by_abundance_presence_absence),"atrF",sep="_")
colnames(mdr1_by_abundance_presence_absence)<-paste(colnames(mdr1_by_abundance_presence_absence),"mdr1",sep="_")
colnames(mdr3_by_abundance_presence_absence)<-paste(colnames(mdr3_by_abundance_presence_absence),"mdr3",sep="_")
colnames(hapE_by_abundance_presence_absence)<-paste(colnames(hapE_by_abundance_presence_absence),"hapE",sep="_")
colnames(hmg1_by_abundance_presence_absence)<-paste(colnames(hmg1_by_abundance_presence_absence),"hmg1",sep="_")
colnames(cdr1B_abcC_by_abundance_presence_absence)<-paste(colnames(cdr1B_abcC_by_abundance_presence_absence),"cdr1B_abcC",sep="_")
colnames(mdr4_by_abundance_presence_absence)<-paste(colnames(mdr4_by_abundance_presence_absence),"mdr4",sep="_")
colnames(cox10_by_abundance_presence_absence)<-paste(colnames(cox10_by_abundance_presence_absence),"cox10",sep="_")
colnames(mdr2_by_abundance_presence_absence)<-paste(colnames(mdr2_by_abundance_presence_absence),"mdr2",sep="_")
colnames(cyp51B_by_abundance_presence_absence)<-paste(colnames(cyp51B_by_abundance_presence_absence),"cyp51B",sep="_")
colnames(AFUA_7G01960_by_abundance_presence_absence)<-paste(colnames(AFUA_7G01960_by_abundance_presence_absence),"AFUA_7G01960",sep="_")

###
#bind and plot together (alphabetically)
all_all_vars<- cbind(AFUA_7G01960_by_abundance_presence_absence,
                     atrF_by_abundance_presence_absence,
                     cdr1B_abcC_by_abundance_presence_absence,
                     cox10_by_abundance_presence_absence,
                     cyp51B_by_abundance_presence_absence,
                     fks1_by_abundance_presence_absence,
                     hapE_by_abundance_presence_absence,
                     hmg1_by_abundance_presence_absence,
                     mdr1_by_abundance_presence_absence,
                     mdr2_by_abundance_presence_absence,
                     mdr3_by_abundance_presence_absence,
                     mdr4_by_abundance_presence_absence)



all_tree_plot <-  gheatmap(tree_plot_me, 
                           all_all_vars,
                                    offset=.1, width=1, low="white", high="black", 
                                    colnames = T, 
                                    colnames_angle = 90,
                                    colnames_position = "top",
                                    font.size = .8,
                                    colnames_offset_y = -8,
                                    #colnames_offset_x = -.1,
                                    color="white") +
  scale_fill_manual(values=c("white", "black")) 

all_tree<- all_tree_plot + geom_tiplab(size = .8, align = TRUE, linesize = .25, offset = 1.4, linetype = 0)
all_tree
#export
ggsave(file="non_Cyp51A_variants.pdf",device="pdf", all_tree, width=18, height=8, units="in")

