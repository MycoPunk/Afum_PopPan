#mapping of Afu (reference Af293) annotated genes to Orthofinder defined OGs to get gli cluster distribution
#last updated: 4.Dec.2021

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

#set wd
#setwd("~")

#load data files
OF.gene_families<-as.data.frame(fread("Orthogroups.tsv")) 
OF.unassigned<- as.data.frame(fread("Orthogroups_UnassignedGenes.tsv")) #this is where your singletons live
OFtoAfu<- as.data.frame(fread("OF_blastP_results_filtered.txt", header = FALSE)) #this is the OG to Afu designations
colnames(OFtoAfu)<- c("Afu_name", "OG_name")
Gli_genes<- as.data.frame(fread("Gli_genes.txt", header = TRUE)) #this is the gli gene mapping file
Afum_grp<-read.delim("clade_map_K3_20Jan2021.txt", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE, check.names = TRUE)

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

View(OFtoAfu)

#assign col names
fmt6_names<- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
colnames(OFtoAfu)<- fmt6_names

#attach Afu annotations
colnames(gene_fam_by_strain_w_anno_ones_num2)[1] <- "OF_gene_fam"
gene_fam_by_strain_w_anno_ones_num2$Af293_reference_gene_name<- OFtoAfu$qseqid[match(gene_fam_by_strain_w_anno_ones_num2$OF_gene_fam, OFtoAfu$sseqid)]

#subset large df to only Afu genes in the gli cluster
Gli_genes_all_counts<- gene_fam_by_strain_w_anno_ones_num2[gene_fam_by_strain_w_anno_ones_num2$Af293_reference_gene_name %in% Gli_genes$Af293_gene_name,]

setdiff(Gli_genes$Af293_gene_name, Gli_genes_all_counts$Af293_reference_gene_name)
#note gliN nad gliF assigned to same cluster 
#Afu6g09730 and Afu6g09720 are both in OG0007452, so we'll consider both here. 

#make presence / absence graph of abundance 
#Gli_genes_all_counts

#load tree data
tree <- read.tree("Afum_260_iq_tree_newick.tre")
#remove the reference from the tree:
tree<- drop.tip(tree, "Af293-REF", trim.internal = TRUE, subtree = FALSE,
                root.edge = 0, rooted = is.rooted(tree), collapse.singles = TRUE,
                interactive = FALSE)
#match tree (for the couple sp. that are  miss-named)
tree$tip.label[tree$tip.label=="F18149"] <- "F18149-JCVI"
tree$tip.label[tree$tip.label=="Afu_343-P/11"] <- "Afu_343-P-11"
tree$tip.label[tree$tip.label=="AFIS1435CDC_6"] <- "AFIS1435_CDC_6"
#root tree based on small outgroup tree (at node 266)
tree<- root(tree, node = 266, resolve.root = FALSE,
            interactive = FALSE, edgelabel = FALSE)

#set colors by group for clade
grA_me<- split(Afum_grp$name_pop_genome, Afum_grp$clade)
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

#plot gli genes onto tree
#subset to only counts
all_vars<- Gli_genes_all_counts[,2:261]
#change to presence / absence
all_vars_presence<- data.frame(sapply(all_vars, gsub, pattern = "1", replacement = "present"))
all_vars_presence_absence<- data.frame(sapply(all_vars_presence, gsub, pattern = "0", replacement = "absent"))
#fix col names X's introduced in last step because R likes to make life difficult
colnames(all_vars_presence_absence)<- colnames(all_vars)
#add rowames
rownames(all_vars_presence_absence)<- Gli_genes_all_counts$Af293_reference_gene_name
#transpose for graphing
all_vars_presence_absence_t<- data.frame(t(all_vars_presence_absence))

#fix names to match OG names
row.names(all_vars_presence_absence_t) <- Afum_grp$name_pop_genome[match(row.names(all_vars_presence_absence_t), Afum_grp$OF_name)]

#get counts
rownames(all_vars)<- Gli_genes_all_counts$Af293_reference_gene_name
rowSums(all_vars)







#plot
gli_tree_plot <-  gheatmap(tree_plot_me, 
                              all_vars_presence_absence_t,
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
ggsave(file="gli_presence_absence.pdf",device="pdf", p, width=8, height=8, units="in")



###protein domain mapping for gliP 
#there are 13 genes in the gli gene cluster
#In addition to the total deletion of gliZ, gliI, and gliJ in 8 strains, gilP has a partial deletion. 
#visualize this relative to known domains in uniprot. 

#get protein features form uniprot via accession number
rel_json_gliP <- drawProteins::get_features("Q4WMJ7")
#make a copy
rel_json_gliP_trunc<- rel_json_gliP

#convert to DF
drawProteins::feature_to_dataframe(rel_json_gliP) -> rel_data_gliP
drawProteins::feature_to_dataframe(rel_json_gliP_trunc) -> rel_json_gliP_trunc


#shrink features to only those not in deletion
rel_json_gliP_trunc<- rel_json_gliP_trunc[rel_json_gliP_trunc$begin < 1278, ]
rel_json_gliP_trunc$order<- 2
rel_json_gliP_trunc[1,5]<- 1278
rel_json_gliP_trunc[1,4]<- 1278
rel_json_gliP_trunc$entryName<- "gliP \nClade 2"
rel_data_gliP$entryName<- "gliP \nAF293"

rel_data_gliP<- rbind(rel_json_gliP_trunc, rel_data_gliP)

#graph gilP
draw_canvas(rel_data_gliP) -> p
p <- draw_chains(p, rel_data_gliP,
                 fill = "grey",   # chain fill color
                 outline = "grey", # chain outline
                 label_size = 2)            
p <- draw_regions(p, rel_data_gliP) # add regions
p <- draw_domains(p, rel_data_gliP, label_domains = FALSE)

p <- p + theme_bw(base_size = 8) + # white background
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank()) +
  theme(axis.ticks = element_blank(),
        axis.text.y = element_blank()) +
  theme(panel.border = element_blank())

p<- p+scale_fill_manual(values = c("#E6F598",
                                   "#66C2A5",
                                   "#3288BD",
                                   "#3E5F83",
                                   "#D53E4F", 
                                   "#683D84",
                                   "#A03A7E"))

# add titles
rel_subtitle <- paste0("Carrier = peptidyl carrier protein (PCP) or T site \ngliP Clade 2 architecture= one ATC module + A: A-T-C-A \ngliP architecture= two ATC modules + T: A-T-C-A-T-C-T")

p <- p + labs(title = "Afu6g09660 / gliP",
              subtitle = rel_subtitle)
p
p + coord_cartesian(xlim = c(-275, max(rel_data_gliP$end)))
