#mapping of Afu (reference Af293) annotated genes to Orthofinder defined OGs was conducted using BLASTP
#GENE lists were procured from FungiDB, using all GO terms available for:
#organonitrogen compound metabolic process (GO: 1901564)
#carbohydrate metabolic process (GO: 0005975)
#organophosphate biosynthetic process (GO: 0090407)

#last updated: 31.Jan.2022

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
library(wesanderson)
library(dplyr)

#set wd
setwd("~/Desktop/Project_Afum_pangenome_3/")

#load data files - hash as needed - the rest of the script does not need editing, except to set the number of top hits to display in the tree (if there are too many results to easily visualize, which is changed at line 366)
#GOI_genes<- as.data.frame(fread("GENES_carbohydrate_metabolic_process.csv", header = TRUE)) 
#GOI_genes<- as.data.frame(fread("GENES_organonitrogen_metabolic_process.csv", header = TRUE))
GOI_genes<- as.data.frame(fread("GENES_organophosphate_metabolic_process.csv", header = TRUE)) 

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

#subset large df to only Afu genes in the target gene list
GOI_genes_all_counts<- gene_fam_by_strain_w_anno_ones_num2[gene_fam_by_strain_w_anno_ones_num2$Af293_reference_gene_name %in% GOI_genes$Af293_gene_name,]

setdiff(GOI_genes$Af293_gene_name, GOI_genes_all_counts$Af293_reference_gene_name)

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

#plot GOI genes onto tree
#subset to only counts
all_vars<- GOI_genes_all_counts[,2:261]
#change to presence / absence
all_vars_presence<- data.frame(sapply(all_vars, gsub, pattern = "1", replacement = "present"))
all_vars_presence_absence<- data.frame(sapply(all_vars_presence, gsub, pattern = "0", replacement = "absent"))
#fix col names X's introduced in last step because R likes to make life difficult
colnames(all_vars_presence_absence)<- colnames(all_vars)
#add rowames
rownames(all_vars_presence_absence)<- GOI_genes_all_counts$Af293_reference_gene_name
#transpose for graphing
all_vars_presence_absence_t<- data.frame(t(all_vars_presence_absence))

#fix names to match OG names
row.names(all_vars_presence_absence_t) <- Afum_grp$name_to_use_in_paper[match(row.names(all_vars_presence_absence_t), Afum_grp$OF_name)]

#get counts
rownames(all_vars)<- GOI_genes_all_counts$Af293_reference_gene_name
rowSums(all_vars)

#if sums are in all, remove and don't map 
all_vars$sums_temp<- rowSums(all_vars)
all_vars_subset<- all_vars[all_vars$sums_temp < 250,]
dim(all_vars)
dim(all_vars_subset)
#subset input file 
all_vars_presence_absence_t_subset<- all_vars_presence_absence_t[,colnames(all_vars_presence_absence_t) %in% rownames(all_vars_subset)]

dim(all_vars_presence_absence_t)
dim(all_vars_presence_absence_t_subset)

#plot
gli_tree_plot <-  gheatmap(tree_plot_me, 
                           all_vars_presence_absence_t_subset,
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

p<- gli_tree_plot + geom_tiplab(size = .8, align = TRUE, linesize = .25, offset = 1, linetype = 0)
p
#ggsave(file="Nitrogen_presence_absence.pdf",device="pdf", p, width=8, height=8, units="in")


###get Organonitrogen_GOterms_significantly differentially abundant between groups. Afu genes significantly different between groups 


#remove cols with zero variance 

#which cols have 0 variance?
#change to 1/0
all_vars_presence_absence_t_presence<- data.frame(sapply(all_vars_presence_absence_t, gsub, pattern = "present", replacement = 1))
all_vars_presence_absence_t_absence<- data.frame(sapply(all_vars_presence_absence_t_presence, gsub, pattern = "absent", replacement = 0))
#add back rownames
rownames(all_vars_presence_absence_t_absence)<- rownames(all_vars_presence_absence_t)



#change to numeric
class(all_vars_presence_absence_t_absence[1,1])
all_vars_presence_absence_t_absence<-as.data.frame(lapply(all_vars_presence_absence_t_absence,as.numeric))
class(all_vars_presence_absence_t_absence[1,1])
rownames(all_vars_presence_absence_t_absence)<- rownames(all_vars_presence_absence_t)


#which(apply(all_vars_presence_absence_t_absence, 2, var) == 0)

#remove cols with 0 variance
all_vars_presence_absence_t_no_zero <- all_vars_presence_absence_t_absence[ - as.numeric(which(apply(all_vars_presence_absence_t_absence, 2, var) == 0))]
ncol(all_vars_presence_absence_t)
ncol(all_vars_presence_absence_t_no_zero)
colVars(all_vars_presence_absence_t_no_zero) #looks good


#attach clade designations
clade_1_names<- Afum_grp[Afum_grp$clade == "1",]
clade_2_names<- Afum_grp[Afum_grp$clade == "2",]
clade_3_names<- Afum_grp[Afum_grp$clade == "3",]


#separate the clades 
clade_1<- all_vars_presence_absence_t_no_zero[rownames(all_vars_presence_absence_t_no_zero) %in% clade_1_names$name_to_use_in_paper,]
clade_2<- all_vars_presence_absence_t_no_zero[rownames(all_vars_presence_absence_t_no_zero) %in% clade_2_names$name_to_use_in_paper,]
clade_3<- all_vars_presence_absence_t_no_zero[rownames(all_vars_presence_absence_t_no_zero) %in% clade_3_names$name_to_use_in_paper,]

#make sure these are the right size
dim(clade_1)
dim(clade_2)
dim(clade_3)

#add clade designators
clade_1$clade<- "clade1"
clade_2$clade<- "clade2"
clade_3$clade<- "clade3"

#merge
input_ON_set<- data.frame(rbind(clade_1, clade_2, clade_3))
dim(input_ON_set)

class(input_ON_set[1,2])

###check variance and normalcy assumptions 
#leveneTest from the car package, look for homogenity of variance.
levene_results<- lapply(subset(input_ON_set, select = -clade), leveneTest, group = as.factor(input_ON_set$clade))
#extract p-vals
extract1<- lapply(levene_results, `[[`, c("Pr(>F)")) 
#pval- less than 0.05? if so, reject= variance is not equal. 
extract2<- data.frame(sapply(extract1, "[[", 1) < 0.05)
table(extract2) #19 are less than 0.05, so not all categories have equal variance. = need to use non parametric test.

#test normalcy of residuals with shapiro.test (if < 0.05, data is normal)
shapiro_results<- apply(input_ON_set[,1:ncol(input_ON_set)-1],2,shapiro.test)
#get p-values
shapiro_results2 <- sapply(shapiro_results, `[`, c("statistic","p.value"))
shapiro_results3<- t(shapiro_results2)
shapiro_results3 #data is not normal - need to run non-parametric test


##run non-parametric (kruskal walliace) test as a loop on each column. 

#need to move clade col to the front
input_ON_set <- input_ON_set[,c(ncol(input_ON_set),1:(ncol(input_ON_set)-1))]

#loop
all_kw_p_vals<- data.frame(p.val =apply(input_ON_set[,-1], 2, function(x) kruskal.test(x,input_ON_set[,1])$p.value))
#bonf adjustment for multiple comp. 
all_kw_p_vals$p.adjust<- p.adjust(all_kw_p_vals$p.val, method = "bonferroni", 
                                  n = nrow(all_kw_p_vals))

p_vals_kw_sig<- all_kw_p_vals[all_kw_p_vals$p.adjust < 0.001,]

names_of_sig_results<- rownames(p_vals_kw_sig)
length(names_of_sig_results)

#top x results
p_vals_kw_sig_ordered<- p_vals_kw_sig[order(p_vals_kw_sig$p.adjust),]
p_vals_kw_sig_top10<- p_vals_kw_sig_ordered[1:7,]
names_of_sig_results_top_10<- rownames(p_vals_kw_sig_top10)

#####
##conduct posthoc tests. 

#subset only significant results
df.subset_and_clade <- input_ON_set[, c("clade",names_of_sig_results)]


#as loop
DT_loop_vals<- apply(df.subset_and_clade[,-1], 2, function(x) dunnTest(x,df.subset_and_clade[,1], method="bonferroni"))
length(DT_loop_vals)

### Compact letter display loop
#instantiate empty data frame
ON_table<- data.frame(matrix(NA, nrow = 3, ncol = 0))

for (i in seq_along(DT_loop_vals)){
  
  name <- (names(DT_loop_vals)[i])
  PT = DT_loop_vals[[i]]$res
  x<- cldList(comparison = PT$Comparison,
              p.value    = PT$P.adj,
              threshold  = 0.05)
  row_to_add<-x$Letter
  ON_table[[name]] <- row_to_add
}
ON_table

######
#make means table 
C1_mean_temp2<- df.subset_and_clade[df.subset_and_clade$clade == "clade1",]
C2_mean_temp2<- df.subset_and_clade[df.subset_and_clade$clade == "clade2",]
C3_mean_temp2<- df.subset_and_clade[df.subset_and_clade$clade == "clade3",]

clade1_means<- round(colMeans(C1_mean_temp2[,-1]),1)
clade2_means<- round(colMeans(C2_mean_temp2[,-1]),1)
clade3_means<- round(colMeans(C3_mean_temp2[,-1]),1)

df_of_ON_means<- as.data.frame(rbind("Clade 1" = clade1_means, "Clade 2" = clade2_means, "Clade 3" = clade3_means))
rownames(ON_table)<- c("Clade 1", "Clade 2", "Clade 3")

#combine post-hoc and mean tables for printing 
mypaste <- function(x,y) paste(x, " (", y, ")", sep="")
table_for_print<- data.frame(mapply(mypaste, df_of_ON_means, ON_table))
rownames(table_for_print)<- c("Clade 1", "Clade 2", "Clade 3")

#write.table(table_for_print, "ON_table.txt", sep="\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

####
####### 

#make figure of the significantly different groups
#get names of significant results. 

names_of_sig_results<- rownames(p_vals_kw_sig)
#subset from the main df
#all
df.subset <- input_ON_set[, names_of_sig_results]
#for top x results

#top x results
#df.subset <- input_ON_set[, names_of_sig_results_top_10]

#Scale if not zero one
#scaled<- apply(df.subset, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))

#map these onto tree. 
ON_tree_plot <-  gheatmap(tree_plot_me, 
                          df.subset,
                          low = "white",
                          high = "black",
                          offset=0.02, width=2, 
                          colnames = T, 
                          colnames_angle = 45,
                          colnames_position = "top",
                          font.size = 2,
                          colnames_offset_y = -6,
                          colnames_offset_x = .05,
                          color="white")
#scale_fill_GOIidis_c(option="F", name="continuous\nvalue") 
#scale_fill_GOIidis_c(option="F", name="normalized\nabundance") 
p<- ON_tree_plot + geom_tiplab(size = .8, align = TRUE, linesize = .25, offset = 2.75, linetype = 0)
p
ggsave(file="PHOS_tree_top10.pdf",device="pdf", p, width=8, height=8, units="in")
