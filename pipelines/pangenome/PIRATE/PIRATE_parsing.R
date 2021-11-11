#A fum pan genome analysis, first pass data investigation of PIRATE results
#started: 4.Oct.2020
#last updated: 27.Sep.2021

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
library(vegan)
packageVersion("vegan") 



setwd("~/Desktop/Project_Afum_pangenome_2/")
#with the double entered Af100-12_9 removed
PIRATE.gene_families<-as.data.frame(fread("PIRATE.gene_families.ordered_27Sep2021.tsv")) 

##basic stats:
#how many gene families are there?
n_gene_fams<- nrow(PIRATE.gene_families)
n_gene_fams
#15,476

#how many strains?

#names of cols to exclude
cols_to_exclude<- c("allele_name",
                    "gene_family",
                    "consensus_gene_name",
                    "consensus_product",
                    "threshold",
                    "alleles_at_maximum_threshold",
                    "number_genomes",
                    "average_dose",
                    "min_dose",
                    "max_dose",
                    "genomes_containing_fissions",
                    "genomes_containing_duplications",
                    "number_fission_loci",
                    "number_duplicated_loci",
                    "no_loci",
                    "products",
                    "gene_names",
                    "min_length(bp)",
                    "max_length(bp)",
                    "average_length(bp)", 
                    "cluster", 
                    "cluster_order")

strains_only<- PIRATE.gene_families[,!names(PIRATE.gene_families) %in% cols_to_exclude]
strain_names<- names(strains_only)

#length(strain_names)
ngenomes<- length(unique(strain_names)) 
ngenomes
#260


#print to cross ref for tree building
#write.table(strain_names, "strain_names_27Sep2021.txt", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#how many of the gene families are in every genome (of 261)
n_gene_fams_core_all<- sum(PIRATE.gene_families$number_genomes == 260)
n_gene_fams_core_all
#3584

#that's what percent out of the total?
(n_gene_fams_core_all*100)/n_gene_fams
#23.15844%

#present in 95% of genomes (n >=248)
cutoff<- round(.95*260)
n_gene_fams_core_w95per<- sum(PIRATE.gene_families$number_genomes >= cutoff)
n_gene_fams_core_w95per
#8600

#that's what percent out of the total?
(n_gene_fams_core_w95per*100)/n_gene_fams
#55.56991%

#how many of the gene families are singletons (accessory)
n_gene_fams_singletons<- sum(PIRATE.gene_families$number_genomes == 1)
n_gene_fams_singletons
#3258

#that's what percent out of the total?
(n_gene_fams_singletons*100)/n_gene_fams
#21.05195%

#how many accessory?)
n_accessory<- n_gene_fams - (n_gene_fams_singletons + n_gene_fams_core_w95per)
n_accessory
#3618

#that's what percent out of the total?
(n_accessory*100)/n_gene_fams
#13.91538%

#get average per genome 
n_accessory / ngenomes


##graph the distribution of gene presence in a gene family (distribution of core to accessory genes)
#plot
gene_fam_totals<-as.data.frame(PIRATE.gene_families$number_genomes)
colnames(gene_fam_totals) <- 'count'

#set groups
gene_fam_totals$group = 0                        
for (i in 1:nrow(gene_fam_totals)){
  if (gene_fam_totals$count[i] == 1) {
   gene_fam_totals$group[i] = "Singleton"
   } else if (gene_fam_totals$count[i] >= .95*261) {
   gene_fam_totals$group[i] = "Core"
   } else {
   gene_fam_totals$group[i] = "Accessory"
   }
}



#plot
#pallet: 
#singelton = "#316A6E",
#accessory = "#BA9141",
#core = "#6E572C")

p <- gene_fam_totals %>%
ggplot(aes(x = count)) +
  geom_bar(aes(color = group), fill = "white",
                 position = "identity") +
  ggtitle("Pangenome Distribution by gene family") +
  ylab("n gene families") + xlab("n genomes in family") +
theme_ipsum() +
  theme(
    plot.title = element_text(size=15), legend.title = element_blank()) +
  scale_color_manual(values = c("#BA9141", "#806633", "#316A6E")) 
p
#ggsave("dist_by_gene_fam.pdf",p, width=6, height=4, units="in")


###plot as doughnut chart 
#make df of totals
fam_dist_df <- data.frame(
  category=c("Singleton", "Accessory", "Core"),
  count=c(n_gene_fams_singletons, n_gene_fams -(n_gene_fams_singletons + n_gene_fams_core_w95per), n_gene_fams_core_w95per)
)

# Compute percentages
fam_dist_df$fraction <- fam_dist_df$count / sum(fam_dist_df$count)
# Compute the cumulative percentages (top of each rectangle)
fam_dist_df$ymax <- cumsum(fam_dist_df$fraction)
# Compute the bottom of each rectangle
fam_dist_df$ymin <- c(0, head(fam_dist_df$ymax, n=-1))
# Compute label position
fam_dist_df$labelPosition <- (fam_dist_df$ymax + fam_dist_df$ymin) / 2
# Compute a good label
fam_dist_df$label <- paste0(fam_dist_df$category, "\n", fam_dist_df$count)
#plot
p<- ggplot(fam_dist_df, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  geom_text( x=1.5, aes(y=labelPosition, label=label, color=category), size=4) + # x here controls label position (inner / outer)
  scale_fill_manual(values=c("#BA9141", "#806633", "#316A6E"))+
  scale_color_manual(values=c("#BA9141", "#806633", "#316A6E"))+
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")
#ggsave("donut.pdf",p, width=3.5, height=3.5, units="in")
p

##which strain has the highest/lowest number of singletons and accessory gene fams? 
gene_fam_by_strain<-as.data.frame(PIRATE.gene_families[,23:ncol(PIRATE.gene_families)])
ncol(gene_fam_by_strain)

##make binary (if gene = 1, if not = 0)
#fill in zeros
gene_fam_by_strain_zeros<- sapply(gene_fam_by_strain, gsub, pattern = "^\\s*$" , replacement = 0 )
#fill in ones
gene_fam_by_strain_ones<- as.data.frame(replace(gene_fam_by_strain_zeros, gene_fam_by_strain_zeros!="0", 1))
#change to numeric
gene_fam_by_strain_ones_num <- mutate_all(gene_fam_by_strain_ones, function(x) as.numeric(as.character(x)))

#subset to remove core genes from accessory and singletons 
all_accessory_1<- gene_fam_by_strain_ones_num[rowSums(gene_fam_by_strain_ones_num) > 1,]
in_95_percent<- .95 * 260
all_accessory<- all_accessory_1[rowSums(all_accessory_1) < in_95_percent,]

#subset to get only singletons 
singletons_only<- gene_fam_by_strain_ones_num[rowSums(gene_fam_by_strain_ones_num) == 1,]

#get average
ave_accessory<- colSums(all_accessory)
mean(ave_accessory) #974.3769

#get average
ave_singletons<- colSums(singletons_only)
mean(ave_singletons) #12.53077


#Which strain has the largest accessory genome (genes not in the core?) 
accessory_by_strain<- as.data.frame(colSums(all_accessory))
colnames(accessory_by_strain) <- "totals"
accessory_by_strain$strain<- row.names(accessory_by_strain)
#get max
accessory_by_strain[which.max(accessory_by_strain$totals),]
#get min
accessory_by_strain[which.min(accessory_by_strain$totals),]
#sort 
accessory_by_strain<- accessory_by_strain[order(accessory_by_strain$totals),]
#mean
mean(accessory_by_strain$totals)

#fix names so that they match 
name_map<-read.delim("clade_map_K3_20Jan2021.txt", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE)
row.names(accessory_by_strain) <- name_map$name_pop_genome[match(row.names(accessory_by_strain), name_map$name_Pan_genome)]
accessory_by_strain$pop_name<- row.names(accessory_by_strain)
#remove "DMC2" for graphing 
accessory_by_strain$pop_name<- sapply(accessory_by_strain$pop_name, gsub, pattern = "DMC2_", replacement = "")


#graph accessory genome size by strain
p<- accessory_by_strain %>%
  mutate(name = fct_reorder(pop_name, totals)) %>%
  ggplot( aes(x=name, y=totals))+
  geom_bar(stat="identity", fill="#BA9141", alpha=5, width=1, position = position_dodge(width=0.6)) +
  xlab("") + ylab("n accessory gene families") +
  ggtitle("accessory genome size by strain") +
  theme(text=element_text(size=9), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(size = 1.5, angle=90, hjust=1, vjust=0.5, margin=margin(-3,0,0,0)), legend.position = "none")
p


#  theme(text=element_text(size=9), 
#        axis.text.x = element_text(size = 2, angle=90, hjust=1), legend.position = "none")+
#  facet_zoom(ylim = c(min(accessory_by_strain$totals), max(accessory_by_strain$totals)), zoom.data = ifelse(a <= 6000,  FALSE))

#ggsave("accessory.pdf",p, width=6.9, height=3, units="in")



#for singletons 
#Which strain has the largest accessory genome (genes not in the core?) 
singletons_only_by_strain<- as.data.frame(colSums(singletons_only))
colnames(singletons_only_by_strain) <- "totals"
singletons_only_by_strain$strain<- row.names(singletons_only_by_strain)
#get max
singletons_only_by_strain[which.max(singletons_only_by_strain$totals),]
#get min
singletons_only_by_strain[which.min(singletons_only_by_strain$totals),]
#how many are zero?
no_singeltons<-data.frame(singletons_only_by_strain[singletons_only_by_strain$totals == 0,])
nrow(no_singeltons)
#how many are 1?
one_singeltons<-data.frame(singletons_only_by_strain[singletons_only_by_strain$totals == 1,])
nrow(one_singeltons)
#get average 
mean(singletons_only_by_strain$totals)
View(singletons_only_by_strain)
#sort 
singletons_only_by_strain<- singletons_only_by_strain[order(singletons_only_by_strain$totals),]
#fix names
row.names(singletons_only_by_strain) <- name_map$name_pop_genome[match(row.names(singletons_only_by_strain), name_map$name_Pan_genome)]
singletons_only_by_strain$pop_name<- row.names(singletons_only_by_strain)
#remove "DMC2" for graphing 
singletons_only_by_strain$pop_name<- sapply(singletons_only_by_strain$pop_name, gsub, pattern = "DMC2_", replacement = "")


#graph singleton genome size by strain
p<- singletons_only_by_strain %>%
  mutate(name = fct_reorder(pop_name, totals)) %>%
  ggplot( aes(x=name, y=totals)) +
  geom_bar(stat="identity", fill="#316A6E", alpha=2, width=1) +
  xlab("") + ylab("n singleton gene families") +
  ggtitle("singleton genome size by strain") +
  theme(text=element_text(size=9), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
  #theme(text=element_text(size=9), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(size = 2, angle=90, hjust=1), legend.position = "none")
p
#ggsave("singletons.pdf",p, width=6.9, height=3, units="in")

#print these for later
#write.table(singletons_only_by_strain, "n_singletons_by_strain.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(accessory_by_strain, "n_accessory_by_strain.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



###graph presence / absence matrix on to big tree

#tree with 260 strains
#tree_me <- read.tree("pop_for_pan2_.SNP.fasttree.tre")
tree_me <- read.tree("Afum_260_iq_tree_newick.tre")

#tree_me <- read.tree("pop_for_pan2_.SNP.mfa20Dec2020.iqtree.tre")

##use mapping file to rename the strains to match the way they appear in the tree

#name_map<-read.delim("clade_map_20Dec2020_3groups.tab", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE)
name_map<-read.delim("clade_map_K3_20Jan2021.txt", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE)

#load mat1-2 map
mat_map<-read.delim("mat_type_map_10Feb2021.txt", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE)
#accessory_by_strain_w_names<- accessory_by_strain
row.names(accessory_by_strain) <- name_map$name_pop_genome[match(row.names(accessory_by_strain), name_map$name_Pan_genome)]
row.names(singletons_only_by_strain) <- name_map$name_pop_genome[match(row.names(singletons_only_by_strain), name_map$name_Pan_genome)]

#attach the clade annotations - note- update this later when you designate clades
singletons_only_by_strain$clade<- name_map$clade[match(row.names(singletons_only_by_strain), name_map$name_pop_genome)]
accessory_by_strain$clade<- name_map$clade[match(row.names(accessory_by_strain), name_map$name_pop_genome)]

#get averages by clade
singeltons_clade1<- singletons_only_by_strain[singletons_only_by_strain$clade == "1",]
mean(singeltons_clade1$totals) # 14
singeltons_clade2<- singletons_only_by_strain[singletons_only_by_strain$clade == "2",]
mean(singeltons_clade2$totals) # 6
singeltons_clade3<- singletons_only_by_strain[singletons_only_by_strain$clade == "3",]
mean(singeltons_clade3$totals) # 16

accessory_clade1<- accessory_by_strain[accessory_by_strain$clade == "1",]
mean(accessory_clade1$totals) # 980
accessory_clade2<- accessory_by_strain[accessory_by_strain$clade == "2",]
mean(accessory_clade2$totals) # 973
accessory_clade3<- accessory_by_strain[accessory_by_strain$clade == "3",]
mean(accessory_clade3$totals) # 909

##is that significant?
#check variance 
res <- bartlett.test(totals ~ clade, data = singletons_only_by_strain)
res
var_sing <-fligner.test(totals ~ clade, data = singletons_only_by_strain)
var_sing
#variance is significantly different between groups - need to use non parametric test

#using perm
permKS(singletons_only_by_strain$clade ~ singletons_only_by_strain$totals, alternative="two.sided", method="exact.mc", 
       control=permControl(nmc=10^4-1))$p.value #not significant

permKS(accessory_by_strain$clade ~ accessory_by_strain$totals, alternative="two.sided", method="exact.mc", 
       control=permControl(nmc=10^4-1))$p.value #not significant

##plot tree
#match tree (for the couple sp. that are  miss-named)
tree_me$tip.label[tree_me$tip.label=="F18149"] <- "F18149-JCVI"
tree_me$tip.label[tree_me$tip.label=="Afu_343-P/11"] <- "Afu_343-P-11"
tree_me$tip.label[tree_me$tip.label=="AFIS1435CDC_6"] <- "AFIS1435_CDC_6"


#rename strains 
#Afum_DMC2_AF100_1_11
#Afum_DMC2_AF100_1_14
#Afum_DMC2_AF100_1_18
#Afum_DMC2_AF100_1_20_OE
#Afum_DMC2_AF100_1_24
#Afum_DMC2_AF100_1_3
#Afum_DMC2_AF100_1_4
#Afum_DMC2_AF100_1_8
#Afum_DMC2_AF100_12_42
#Afum_DMC2_AF100_12_9
#Afum_DMC2_AF100_2B
#Afum_DMC2_AF100_3B
#Afum_DMC2_AF100_4B
#Afum_DMC2_AF100_5B
#Afum_DMC2_AF100_6B
#Afum_DMC2_AF100_7B
#Afum_DMC2_AF100_8B
#Afum_DMC2_AF100_9B
#Afum_DMC2_AF1001_15


#remove the reference from the tree:
tree_me<- drop.tip(tree_me, "Af293-REF", trim.internal = TRUE, subtree = FALSE,
         root.edge = 0, rooted = is.rooted(tree_me), collapse.singles = TRUE,
         interactive = FALSE)

#split by clade
grA_me<- split(row.names(accessory_by_strain), accessory_by_strain$clade)

#spllit by mat type
grA_me_mat<- split(mat_map$name_pop_genome, mat_map$MAT_type)

#set colors by group for clade
tree_grA_me <- ggtree::groupOTU(tree_me, grA_me)
str(tree_grA_me)
levels(attributes(tree_grA_me)$group) 
levels(attributes(tree_grA_me)$group)[1] <- "1"
# Reorder factor levels if needed (this controls the order in the legend)
#attributes(tree_grA_me)$group <- factor(x = attributes(tree_grA_me)$group, 
#                                        levels = c("clade1", "clade2", "clade3", "clade4", "clade5", "clade6", "clade7", "clade8"))

attributes(tree_grA_me)$group <- factor(x = attributes(tree_grA_me)$group, 
                                        levels = c("1", "2", "3"))


#set colors
my_cols_me <- c(clade1 = "#56326E",
                clade2 ="#ED7F6F",
                clade3 = "#ABA778")

#my_cols_me <- c(clade1 = "#382147",
#                clade2 ="#7D7A70",
#                clade3 = "#ABA778",
#                clade4 = "#F2CC35")
                #clade5 = "#F2CC85",
                #clade7 = "#FFCA98",
                #clade6 = "#F7A583",
                #clade7 = "#ED7F6F",
                #clade8 = "#D4494E")
 

names(my_cols_me) <- levels(attributes(tree_grA_me)$group)
scales::show_col(my_cols_me); my_cols_me


###set colors by group for mat type 
tree_grA_me_mat <- ggtree::groupOTU(tree_me, grA_me_mat)
str(tree_grA_me_mat)
tree_grA_me_mat2 <- ggtree::groupOTU(tree_grA_me_mat, grA_me)
str(tree_grA_me_mat2)
levels(attributes(tree_grA_me_mat)$group) 
#levels(attributes(tree_grA_me_mat)$group)[1] <- "1"
attributes(tree_grA_me_mat)$group <- factor(x = attributes(tree_grA_me_mat)$group, 
                                        levels = c("MAT1","MAT2","unknown"))
my_cols_me_mat <- c(MAT1 = "#BEBDBD",
                    MAT2 ="#404041",
                    unknown = "#D4494E")

names(my_cols_me_mat) <- levels(attributes(tree_grA_me_mat)$group)
scales::show_col(my_cols_me_mat); my_cols_me_mat



#optional ultrametric tree 
#tree_grA_me_ultra<- force.ultrametric(tree_grA_me)

#root tree



#simple plot
tree_plot_me <- 
  ggtree(tr = tree_grA_me, 
         # color by group attribute, check str(tree_grA_me)
         mapping = aes(color = group), 
         layout  = 'circular', 
         branch.length = 'none', 
#  geom_treescale(x=3, y=NULL, color = "white") +
  # set line thickness
  size = .3) +
  # adjust coloring of main groups
  scale_color_manual(name = 'Clade', values = my_cols_me) + 
  theme(legend.title=element_text(size=9), # The title of legend 
      legend.text=element_text(size=7))
  #xlim(NA, NA)
#  guides(color = guide_legend(override.aes = list(size = 4))) 


# plot and ddd the tip labels
tree_plot_me + geom_tiplab(size = 1, align = TRUE, linesize = .1, linetype = 0) 

#ggsave(file="pan_genome_tree_cladogram.png",device="png")

#flip node
plotTree(tree_grA_me,node.numbers=T) #264?
plot(tree_grA_me, use.edge.length=TRUE, cex = .2, label.offset = 0)
nodelabels(cex=.5)
tree_to_map_on<- 
flip(tree_plot_me, 266, 325)+ geom_tiplab(size = 0, align = TRUE, linesize = .25, linetype = 3)



####map matting type
mat_map

#base plot
tree_plot_me_mat <- 
  ggtree(tr = tree_grA_me_mat, 
         # color by group attribute, check str(tree_grA_me)
         mapping = aes(color = group), 
         layout  = 'circular', 
  branch.length = 'none') + 
  geom_treescale(x=3, y=NULL, color = "white") +
  # set line thickness
  #  size = 1)
  # adjust coloring of main groups
  scale_color_manual(name = 'Clade', values = my_cols_me_mat) + 
  #scale_shape_manual("Clade", values = my_cols_me_mat, breaks=c("MAT1", "MAT2", "Unknown"), labels=c("O","R","N"))+
  theme(legend.title=element_text(size=9), # The title of legend 
        legend.text=element_text(size=7))
#guides(color = guide_legend(override.aes = list(size = 4))) 


tree_plot_me_mat + geom_tiplab(size = 1, align = TRUE, linesize = 0, linetype = 0) +
  geom_tippoint(aes(x=x+16, color=group), size=.50)

#combine annotations for clade and MAT type into one tree
pi<- tree_plot_me + geom_tiplab(size = .6, align = TRUE, linesize = .20, linetype = NA) 
p <- pi %<+% mat_map + geom_tippoint(aes(x=x+14,color=MAT_type), size=.30) + scale_color_manual(values=c("#56326E", "#ED7F6F","#ABA778", "#BEBDBD", "#404041", "#D4494E"))
  #theme(legend.title=element_text(size=10), # The title of legend 
  #      legend.text=element_text(size=7))+
  #guides(color = guide_legend(override.aes = list(linetype = c(1, 1, 1, 0, 0, 0), size = c(1,1,1,2,2,2), shape = c(NA, NA,NA, 16, 16, 16))))


plot(p)
#####


#shrink the the min val of each data set
accessory_by_strain$totals<- accessory_by_strain$totals - min(accessory_by_strain$totals)
singletons_only_by_strain$totals<- singletons_only_by_strain$totals - min(singletons_only_by_strain$totals)

##to visualize, re-scale both data sets on a scale of 0-1
#function
scale <- function(x){(x-min(x))/(max(x)-min(x))}
accessory_by_strain$totals<- scale(accessory_by_strain$totals)
singletons_only_by_strain$totals<- scale(singletons_only_by_strain$totals)

#join annotation data (accessory)
accessory_by_strain_input<- data.table(cbind((tip_lbs = row.names(accessory_by_strain)), clade = (accessory_by_strain$clade), val = (accessory_by_strain$totals)))
tree_plot_me <- tree_plot_me %<+% data.table(accessory_by_strain_input) 

#join annotation data (singleton)
singletons_only_by_strain_input<- data.table(cbind((tip_lbs = row.names(singletons_only_by_strain)), clade = (singletons_only_by_strain$clade), val2 = (singletons_only_by_strain$totals)))
tree_plot_me <- tree_plot_me %<+% data.table(singletons_only_by_strain_input) 

#process updated tree view data
tree_dt_me <- data.table(tree_plot_me$data)


# select only the tip labels and order by coord y
tree_dt_me <- tree_dt_me[isTip == TRUE][order(y)]
tree_dt_me

# Plot the tree with circular barplot and MAT type all together
my_factor <- 2
x_base_me <- max(tree_dt_me$x) + abs(min(as.numeric(tree_dt_me$val), na.rm = TRUE))*my_factor + 18
# Define variable to control the x coordinates of segments & labels
my_x_me <- x_base_me + max(as.numeric(tree_dt_me$val), na.rm = TRUE)*my_factor + 5.5

#for singletons
x_base_me2 <- max(tree_dt_me$x) + abs(min(as.numeric(tree_dt_me$val2), na.rm = TRUE))*my_factor + 25.5
# Define variable to control the x coordinates of segments & labels
my_x_me2 <- x_base_me + max(as.numeric(tree_dt_me$val2), na.rm = TRUE)*my_factor + 10.5

fill_in_value <- 0.2

tree_bars_me <- 
  #tree_labeled_me +
  p +
  #Add a disc to plot bars on top of it
  geom_rect(data = tree_dt_me,
            aes(xmin = x_base_me + min(as.numeric(val), na.rm = TRUE)*my_factor,
                ymin = y - 0.1,
                xmax = x_base_me + max(as.numeric(val))*(my_factor + 3),
                # Add also fill_in_value to `ymax` to force a complete circle
                ymax = max(y) + fill_in_value), 
            color = NA, # set NA so to avoid coloring borders
            fill = "#F0F1F2",
            alpha = 0.5) +
  # Add bars for each tip label
  geom_rect(data = tree_dt_me,
            aes(xmin = x_base_me,
                ymin = y - 0.1,
                #change the +3 here to increase the relative bar height 
                xmax = x_base_me + as.numeric(val)*(my_factor + 3),
                ymax = y + 0.1),
            fill = "#BA9141",
            show.legend = FALSE,
            # no borders? color for now 
            color = "#BA9141") +
  #add second set of rings for singletons
  #ring
  geom_rect(data = tree_dt_me,
            aes(xmin = x_base_me2 + min(as.numeric(val2), na.rm = TRUE)*my_factor,
                ymin = y - 0.1,
                #xmax = x_base_me2 + max(as.numeric(val2)*my_factor, na.rm = TRUE)*my_factor,
                xmax = x_base_me2 + max(as.numeric(val2))*(my_factor + 3),
                # Add also fill_in_value to `ymax` to force a complete circle
                ymax = max(y) + fill_in_value), 
            color = NA, # set NA so to avoid coloring borders
            fill = "#F0F1F2",
            alpha = 0.5) +
  #bars
  geom_rect(data = tree_dt_me,
            aes(xmin = x_base_me2,
                ymin = y - 0.1,
                xmax = x_base_me2 + as.numeric(val2)*(my_factor + 3) ,
                ymax = y + 0.1),
            #width = 1,
            fill = "#316A6E",
            show.legend = FALSE,
            # no borders 
            color = "#316A6E") +
#theme(legend.position = 'bottom',
  theme(legend.position = c(.5,.009), legend.direction = "horizontal",
      legend.background = element_rect(),
      legend.key = element_blank(), # removes the border
      legend.key.size = unit(.5, 'cm'), # sets overall area/size of the legend
      legend.text = element_text(size = 7), # text size
      title = element_text(size = 7),
      #legend.key.size = unit(1.5, "cm"),
      legend.key.width = unit(.5,"cm"),
      legend.key.height = unit(.2,"cm"))

#tree_bars_me + geom_tiplab(size = .6, offset = .03, align = TRUE, linesize = .04, linetype = 1)
tree_all_data<- tree_bars_me + geom_tiplab(size = .6, align = TRUE, linesize = .20, linetype = NA)+ guides(fill = guide_legend(ncol = 3))+
theme(legend.title=element_text(size=10), # The title of legend 
       legend.text=element_text(size=7))+
 guides(color = guide_legend(override.aes = list(linetype = c(1, 1, 1, 0, 0, 0), size = c(1,1,1,2,2,2), shape = c(NA, NA,NA, 16, 16, 16)), title = "", nrow = 3)) 

tree_all_data 

#ggsave(file="pan_genome_tree_test.pdf",device="pdf")



###are there unique gene families (or unique losses) by clade?
#View(accessory_by_strain)
#count data per family

#need to match annotation data to the all_accessory data 
#create 0/1 df 
PIRATE.gene_families$
gene_fam_by_strain_w_anno$cluster
gene_fam_by_strain_w_anno<-as.data.frame(PIRATE.gene_families[,c(1:4,21:ncol(PIRATE.gene_families))])
##make binary (if gene = 1, if not = 0)
#fill in zeros
gene_fam_by_strain_w_anno_zeros<- sapply(gene_fam_by_strain_w_anno, gsub, pattern = "^\\s*$" , replacement = 0 )
#fill in ones
gene_fam_by_strain_w_anno_ones<- as.data.frame(replace(gene_fam_by_strain_w_anno_zeros[,7:ncol(gene_fam_by_strain_w_anno_zeros)], gene_fam_by_strain_zeros!="0", 1))
#change to numeric
gene_fam_by_strain_w_anno_ones_num <- mutate_all(gene_fam_by_strain_w_anno_ones, function(x) as.numeric(as.character(x)))
#bind annotations
gene_fam_by_strain_w_anno_ones_num2<- cbind(gene_fam_by_strain_w_anno_zeros[,1:4], gene_fam_by_strain_w_anno_ones_num)


#subset to remove core genes and singletons 
#remove core
all_accessory_anno<- gene_fam_by_strain_w_anno_ones_num2[rowSums(gene_fam_by_strain_w_anno_ones_num2[,5:ncol(gene_fam_by_strain_w_anno_ones_num2)]) < 248,]
dim(all_accessory_anno)

#remove singletons
all_accessory_anno_no_sing<- all_accessory_anno[rowSums(all_accessory_anno[,5:ncol(all_accessory_anno)]) != 1,]
dim(all_accessory_anno_no_sing)

#split list of strain names by clade (as factors)
#accessory_by_strain_factor_list<- as.list(accessory_by_strain$strain, accessory_by_strain$clade)
accessory_by_strain_list <- setNames(accessory_by_strain$clade, accessory_by_strain$strain)
#accessory_by_strain_list_factor<- as.factor(accessory_by_strain_list)

#split the count data frame by the factor list
clade1_names<- accessory_by_strain_list[accessory_by_strain_list == "1"]
clade1<- cbind(all_accessory_anno_no_sing[,1:4], all_accessory_anno_no_sing[,colnames(all_accessory_anno_no_sing) %in% names(clade1_names)])

clade2_names<- accessory_by_strain_list[accessory_by_strain_list == "2"]
clade2<- cbind(all_accessory_anno_no_sing[,1:4], all_accessory_anno_no_sing[,colnames(all_accessory_anno_no_sing) %in% names(clade2_names)])

clade3_names<- accessory_by_strain_list[accessory_by_strain_list == "3"]
clade3<- cbind(all_accessory_anno_no_sing[,1:4], all_accessory_anno_no_sing[,colnames(all_accessory_anno_no_sing) %in% names(clade3_names)])


#are there gene families that are exclusive to clade 1?
#clade1_rowsums<- clade1[rowSums(clade1[,5:ncol(clade1)]) >0,]

#calculate rowsums

#QC make sure they match 
#length(names(clade7)) -5 
#length(clade7_names)  #yep

clade1$totals<-  rowSums(clade1[,5:ncol(clade1)])
clade2$totals<-  rowSums(clade2[,5:ncol(clade2)])
clade3$totals<-  rowSums(clade3[,5:ncol(clade3)])


#get fams exclusive to clade1
exclusive_to_clade1<- clade1[(clade1$totals > 0) & 
                               (clade2$totals == 0) &
                               (clade3$totals == 0),]

dim(exclusive_to_clade1) 

aveclade1<- nrow(exclusive_to_clade1) / length(clade1_names)
round(aveclade1) #ave is 5

#there are 1062 accessory gene fams exclusive to clade 1
length(clade1_names) #there are 200 strains in clade 1
sort(exclusive_to_clade1$totals) # top hit is distributed in 16 of the 201 isolates in clade1
length(clade1_names)*.95 
#get all present in more than 1/3 of the isolates 
exclusive_to_clade1_of_note<- exclusive_to_clade1[exclusive_to_clade1$totals > length(clade1_names)*.90,]
dim(exclusive_to_clade1_of_note)
#there is one calde defining gene fam gain (g005555 )
exclusive_to_clade1_of_note

#exclusive to clade2
exclusive_to_clade2<- clade2[(clade2$totals > 0) & 
                               (clade1$totals == 0) &
                               (clade3$totals == 0),]

aveclade2<- nrow(exclusive_to_clade2) / length(clade2_names)
round(aveclade2) #ave is 3

dim(exclusive_to_clade2) #there are 131 gene fams exclusive to clade2
length(clade2_names) #there are 46 isolates in clade 2
sort(exclusive_to_clade2$totals) #some are very wide spread, appearing in up to 43 of the 46 isolates
exclusive_to_clade2_of_note<- exclusive_to_clade2[exclusive_to_clade2$totals > length(clade2_names)*.95,]
dim(exclusive_to_clade2_of_note) #0 clade defining gains in Clade 2


#exclusive to clade3
exclusive_to_clade3<- clade3[(clade3$totals > 0) & 
                               (clade1$totals == 0) &
                               (clade2$totals == 0),]


dim(exclusive_to_clade3) 
aveclade3<- nrow(exclusive_to_clade3) / length(clade3_names)
round(aveclade3) #ave is 7

length(clade3_names)
sort(exclusive_to_clade3$totals)
#there are 102 accessory gene fams exclusive to clade3, some in all 15 isoaltes
exclusive_to_clade3_of_note<- exclusive_to_clade3[exclusive_to_clade3$totals > length(clade3_names)*.95,]
dim(exclusive_to_clade3_of_note) #14 are present in more than > 95% of all isolates


##graph distribution of accesory genes by clade 
#calculate percentage of genomes of all in clade that unique gene family is found in
exclusive_to_clade1$percent<- exclusive_to_clade1$totals / length(clade1_names) * 100
exclusive_to_clade2$percent<- exclusive_to_clade2$totals / length(clade2_names) * 100
exclusive_to_clade3$percent<- exclusive_to_clade3$totals / length(clade3_names) * 100

#format data in longform for graphing 
clade1_for_box<- data.frame(cbind(name = "Clade 1", value = round(exclusive_to_clade1$percent, digits = 2)))
clade2_for_box<- data.frame(cbind(name = "Clade 2", value = round(exclusive_to_clade2$percent, digits = 2)))
clade3_for_box<- data.frame(cbind(name = "Clade 3", value = round(exclusive_to_clade3$percent, digits = 2)))

box_df<- data.frame(rbind(clade1_for_box, clade2_for_box, clade3_for_box))

#plot
p<- box_df %>%
  ggplot( aes(x=name, y=as.numeric(value), fill=name, colour = factor(name))) +
geom_point(size=1, position=position_jitter(width=0.3, height=0))+
  scale_colour_manual(values=c("#56326E", "#ED7F6F","#ABA778"))+
  theme_bw() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11))+
  theme(axis.title.x=element_blank(), panel.background = element_blank())+
  scale_y_continuous(name="% of genomes in a clade", limits=c(0, 100))
  # + geom_hline(yintercept = 50)
p

#ggsave("distribution_of_accessory_genes_by_clade.pdf", p, width=20, height=20, 
#       device = "pdf", units = "cm")



#generate dfs of unique gene fams above the 50% cut off
#exclusive_to_clade1_above50<- exclusive_to_clade1[exclusive_to_clade1$percent > 50,]
#exclusive_to_clade2_above50<- exclusive_to_clade2[exclusive_to_clade2$percent > 50,]
#exclusive_to_clade3_above50<- exclusive_to_clade3[exclusive_to_clade3$percent > 50,]

#size of each search 
#nrow(exclusive_to_clade1_above50) #3
#nrow(exclusive_to_clade2_above50) #10
#nrow(exclusive_to_clade3_above50) #51

#annotations?
#exclusive_to_clade1_above50$consensus_product #none
#exclusive_to_clade2_above50$consensus_product #none
#exclusive_to_clade3_above50$consensus_product #none

#are these annotated in AF293?
#exclusive_to_clade1_above50$Afum_AF293 #only one is present in AF293 #g005563

#print these to pull out these gene fams for annotation in bash
#write.table(exclusive_to_clade1_above50, "exclusive_to_clade1_above50.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(exclusive_to_clade2_above50, "exclusive_to_clade2_above50.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(exclusive_to_clade3_above50, "exclusive_to_clade3_above50.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#Print all exclusive gene fams
#print these to pull out these gene fams for annotation in bash
#write.table(exclusive_to_clade1, "exclusive_to_clade1_all.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(exclusive_to_clade2, "exclusive_to_clade2_all.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(exclusive_to_clade3, "exclusive_to_clade3_all.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


##do the same for losses
#get fams lost exclusive to clade1
exclusive_to_clade1_lost<- clade1[(clade1$totals == 0) & 
                               (clade2$totals > 0) &
                               (clade3$totals > 0),]

dim(exclusive_to_clade1_lost) #there are 16 accessory gene fams exclusively lost in Clade 1

aveclade1_lost<- nrow(exclusive_to_clade1_lost) / length(clade1_names)
round(aveclade1_lost) #that's well less than one per genome. 
aveclade1_lost #0.07960199 to be exact

#all have been lost in all isolates of clade 1 - per the previous search perams. 
#how abundant are they in the other clades though? 
#GO BACK HERE AND ANSWER THIS - we want to know what's been lost in a clade but is present in the majority of other isolates in the other two clades. 

all_totals<- cbind(clade1_totals= clade1$totals, clade2_totals = clade2$totals, clade3_totals = clade3$totals)
row.names(all_totals)<- clade1$gene_family
all_totals<- data.frame(all_totals)

#add percent abundance in other clades 
all_totals$percent_abundance_in_Clade1<- all_totals$clade1_totals / length(clade1_names) * 100
all_totals$percent_abundance_in_Clade2<- all_totals$clade2_totals / length(clade2_names) * 100
all_totals$percent_abundance_in_Clade3<- all_totals$clade3_totals / length(clade3_names) * 100


#subset to losses in clade1
all_totals_lost_in_clade1<- all_totals[(all_totals$clade1_totals == 0) & 
                               (all_totals$clade2_totals > 0) &
                               (all_totals$clade3_totals > 0),]

dim(all_totals_lost_in_clade1) #16 of these

all_totals_lost_in_clade1_of_interest<- all_totals_lost_in_clade1[(all_totals_lost_in_clade1$percent_abundance_in_Clade2 > 90) |
                                                                    (all_totals_lost_in_clade1$percent_abundance_in_Clade3 > 90),]

dim(all_totals_lost_in_clade1_of_interest) #5 at greater than 90% 


#get the abundance in clades 2 and 3 for the genes identified in exclusive_to_clade1_lost$gene_family

#get gene fams lost exclusive to clade2
exclusive_to_clade2_lost<- clade2[(clade2$totals == 0) & 
                               (clade1$totals > 0) &
                               (clade3$totals > 0),]

dim(exclusive_to_clade2_lost) # there are 186 gene fams lost in Clade 2 that are present in Clades 1 and 3
aveclade2_lost<- nrow(exclusive_to_clade2_lost) / length(clade2_names)
round(aveclade2_lost) #that's an average of 4 genes per strain
aveclade3_lost<- nrow(exclusive_to_clade3_lost) / length(clade3_names)
round(aveclade3_lost) #that's a crazy 58 genes lost per strain.

#subset to losses in clade2 - total losses
all_totals_lost_in_clade2<- all_totals[(all_totals$clade1_totals > 0) & 
                                         (all_totals$clade2_totals == 0) &
                                         (all_totals$clade3_totals > 0),]

dim(all_totals_lost_in_clade2) #187 LOTS! subset to only presence/absence


all_totals_lost_in_clade2_of_interest<- all_totals_lost_in_clade2[(all_totals_lost_in_clade2$percent_abundance_in_Clade1 > 90) |
                                                                    (all_totals_lost_in_clade2$percent_abundance_in_Clade3 >90),]

dim(all_totals_lost_in_clade2_of_interest) #25 at 90%



#get gene fams lost exclusive to clade3
exclusive_to_clade3_lost<- clade3[(clade3$totals == 0) & 
                               (clade1$totals > 0) &
                               (clade2$totals > 0),]

dim(exclusive_to_clade3_lost) #875 lost in Clade 3 but present in Clades 1 and 2. 




#subset to losses in clade3
all_totals_lost_in_clade3<- all_totals[(all_totals$clade1_totals > 0) & 
                                         (all_totals$clade2_totals > 0) &
                                         (all_totals$clade3_totals == 0),]

dim(all_totals_lost_in_clade3) #wow there are a lot of these, lets subset to the ones at 100% presence/absence

all_totals_lost_in_clade3_of_interest<- all_totals_lost_in_clade3[(all_totals_lost_in_clade3$percent_abundance_in_Clade1 > 90) |
                                                                  (all_totals_lost_in_clade3$percent_abundance_in_Clade2 > 90),]

dim(all_totals_lost_in_clade3_of_interest) #115 - 10 with 100% present in the other two clades

###
#bind tables of interest, and attach annotations via. Jason's output
#first add row names
all_totals_lost_in_clade1_of_interest$gene_fam<- rownames(all_totals_lost_in_clade1_of_interest)
all_totals_lost_in_clade1_of_interest$greater_than_95<- (all_totals_lost_in_clade1_of_interest$percent_abundance_in_Clade2 > 95) & (all_totals_lost_in_clade1_of_interest$percent_abundance_in_Clade3 > 95)

all_totals_lost_in_clade2_of_interest$gene_fam<- rownames(all_totals_lost_in_clade2_of_interest)
all_totals_lost_in_clade2_of_interest$greater_than_95<- (all_totals_lost_in_clade2_of_interest$percent_abundance_in_Clade1 > 95) & (all_totals_lost_in_clade2_of_interest$percent_abundance_in_Clade3 > 95)

all_totals_lost_in_clade3_of_interest$gene_fam<- rownames(all_totals_lost_in_clade3_of_interest)
all_totals_lost_in_clade3_of_interest$greater_than_95<- (all_totals_lost_in_clade3_of_interest$percent_abundance_in_Clade2 > 95) & (all_totals_lost_in_clade3_of_interest$percent_abundance_in_Clade1 > 95)



#pivot long
lostinClade1.0<- all_totals_lost_in_clade1_of_interest %>%
  pivot_longer(
    cols = c(percent_abundance_in_Clade2, percent_abundance_in_Clade3))
lostinClade1.0$group<- "lost_in_Clade1"

lostinClade2.0<- all_totals_lost_in_clade2_of_interest %>%
  pivot_longer(
    cols = c(percent_abundance_in_Clade1, percent_abundance_in_Clade3))
lostinClade2.0$group<- "lost_in_Clade2"

lostinClade3.0<- all_totals_lost_in_clade3_of_interest %>%
  pivot_longer(
    cols = c(percent_abundance_in_Clade1, percent_abundance_in_Clade2))
lostinClade3.0$group<- "lost_in_Clade3"

#subset
lostinClade1.1<- data.frame(cbind(value = lostinClade1.0$value, group =lostinClade1.0$name, gene_fam=lostinClade1.0$gene_fam, greater_than_95 = lostinClade1.0$greater_than_95))
lostinClade1.1$value<- as.numeric(lostinClade1.1$value)
lostinClade1.1$group<- as.factor(lostinClade1.1$group)
lostinClade1.1$gene_fam<- as.factor(lostinClade1.1$gene_fam)
lostinClade1.1$greater_than_95<- as.factor(lostinClade1.1$greater_than_95)
lostinClade1.1$lost_in<- as.factor("lost_in_Clade1")

lostinClade2.1<- data.frame(cbind(value = lostinClade2.0$value, group =lostinClade2.0$name, gene_fam=lostinClade2.0$gene_fam, greater_than_95 = lostinClade2.0$greater_than_95))
lostinClade2.1$value<- as.numeric(lostinClade2.1$value)
lostinClade2.1$group<- as.factor(lostinClade2.1$group)
lostinClade2.1$gene_fam<- as.factor(lostinClade2.1$gene_fam)
lostinClade2.1$greater_than_95<- as.factor(lostinClade2.1$greater_than_95)
lostinClade2.1$lost_in<- as.factor("lost_in_Clade2")

lostinClade3.1<- data.frame(cbind(value = lostinClade3.0$value, group =lostinClade3.0$name, gene_fam=lostinClade3.0$gene_fam), greater_than_95 = lostinClade3.0$greater_than_95)
lostinClade3.1$value<- as.numeric(lostinClade3.1$value)
lostinClade3.1$group<- as.factor(lostinClade3.1$group)
lostinClade3.1$gene_fam<- as.factor(lostinClade3.1$gene_fam)
lostinClade3.1$greater_than_95<- as.factor(lostinClade3.1$greater_than_95)
lostinClade3.1$lost_in<- as.factor("lost_in_Clade3")

#get total of both above %95
sum(lostinClade1.1$greater_than_95 =="TRUE") /2 #0
sum(lostinClade2.1$greater_than_95 =="TRUE") /2 #1
sum(lostinClade3.1$greater_than_95 =="TRUE") /2 #37


##ggplot
##ggplot way to color lines
clade1_loss<- ggplot(lostinClade1.1, aes(x = group, y = value)) +
  #geom_boxplot(aes(fill = group), alpha = 0.2, col = "grey") +
  geom_line(aes(group = gene_fam, col = greater_than_95), alpha=0.4) +
  scale_colour_manual(values = c("grey", "#ED7F6F", "#ABA778", "#E1AF00"))+
  geom_point(aes(col = group)) +
  ylab("% of genomes in a clade")+
  xlab("Lost in Clade 1")+
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

clade1_loss


clade2_loss<- ggplot(lostinClade2.1, aes(x = group, y = value)) +
  #geom_boxplot(aes(fill = group), alpha = 0.2, col = "grey") +
  geom_line(aes(group = gene_fam, col = greater_than_95), alpha=0.4) +
  scale_colour_manual(values = c("grey", "#56326E", "#ABA778", "#E1AF00"))+
  geom_point(aes(col = group)) +
  ylab("% of genomes in a clade")+
  xlab("Lost in Clade 2")+
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
clade2_loss

clade3_loss<- ggplot(lostinClade3.1, aes(x = group, y = value)) +
  #geom_boxplot(aes(fill = group), alpha = 0.2, col = "grey") +
  geom_line(aes(group = gene_fam, col = greater_than_95), alpha=0.4) +
  scale_colour_manual(values = c("grey", "#56326E", "#ED7F6F", "#E1AF00"))+
  geom_point(aes(col = group)) +
  ylab("% of genomes in a clade")+
  xlab("Lost in Clade 3")+
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
clade3_loss

#plot together

all<- ggarrange(clade1_loss, clade2_loss, clade3_loss + rremove("x.text"), 
                labels = c("A", "B", "C"),
                ncol = 3, nrow = 1,
                heights = c(2, 2, 2))
all

#ggsave("distribution_of_losses_by_clade.pdf", all, width=20, height=20, 
#       device = "pdf", units = "cm")





##############################
#write to file for analysis in bash
just_singeltons<- gene_fam_by_strain_w_anno_ones_num2[rowSums(gene_fam_by_strain_w_anno_ones_num2[,5:ncol(gene_fam_by_strain_w_anno_ones_num2)]) == 1,]
just_core<- gene_fam_by_strain_w_anno_ones_num2[rowSums(gene_fam_by_strain_w_anno_ones_num2[,5:ncol(gene_fam_by_strain_w_anno_ones_num2)]) >= 248,]
just_acessory_1<- gene_fam_by_strain_w_anno_ones_num2[rowSums(gene_fam_by_strain_w_anno_ones_num2[,5:ncol(gene_fam_by_strain_w_anno_ones_num2)]) > 1,]
just_acessory<- just_acessory_1[rowSums(just_acessory_1[,5:ncol(just_acessory_1)]) < 248,]

#write.table(just_singeltons, "acessory_gene_numbers.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(just_core, "core_gene_numbers.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(just_acessory, "singleton_gene_numbers.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##############################


##graph distribution of accesory genes LOST by clade
#calculate percentage of genomes of all in clade that unique gene family is found in
exclusive_to_clade1$percent<- exclusive_to_clade1$totals / length(clade1_names) * 100
exclusive_to_clade2$percent<- exclusive_to_clade2$totals / length(clade2_names) * 100
exclusive_to_clade3$percent<- exclusive_to_clade3$totals / length(clade3_names) * 100

#format data in longform for graphing 
clade1_for_box<- data.frame(cbind(name = "Clade 1", value = round(exclusive_to_clade1$percent, digits = 2)))
clade2_for_box<- data.frame(cbind(name = "Clade 2", value = round(exclusive_to_clade2$percent, digits = 2)))
clade3_for_box<- data.frame(cbind(name = "Clade 3", value = round(exclusive_to_clade3$percent, digits = 2)))

box_df<- data.frame(rbind(clade1_for_box, clade2_for_box, clade3_for_box))

#plot
p<- box_df %>%
  ggplot( aes(x=name, y=as.numeric(value), fill=name, colour = factor(name))) +
  geom_point(size=1, position=position_jitter(width=0.3, height=0))+
  scale_colour_manual(values=c("#56326E", "#ED7F6F","#ABA778"))+
  theme_bw() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11))+
  theme(axis.title.x=element_blank(), panel.background = element_blank())+
  scale_y_continuous(name="% of genomes in a clade", limits=c(0, 100))+
  geom_hline(yintercept = 50)
p

#ggsave("distribution_of_accessory_genes_by_clade.pdf", p, width=20, height=20, 
#       device = "pdf", units = "cm")



#generate dfs of unique gene fams above the 50% cut off
exclusive_to_clade1_above50<- exclusive_to_clade1[exclusive_to_clade1$percent > 50,]
exclusive_to_clade2_above50<- exclusive_to_clade2[exclusive_to_clade2$percent > 50,]
exclusive_to_clade3_above50<- exclusive_to_clade3[exclusive_to_clade3$percent > 50,]

#size of each search 
nrow(exclusive_to_clade1_above50) #3
nrow(exclusive_to_clade2_above50) #10
nrow(exclusive_to_clade3_above50) #51

#annotations?
exclusive_to_clade1_above50$consensus_product #none
exclusive_to_clade2_above50$consensus_product #none
exclusive_to_clade3_above50$consensus_product #none

#are these annotated in AF293?
exclusive_to_clade1_above50$Afum_AF293 #only one is present in AF293 #g005563


##graph these results
##subset the clade-trees from the larger tree
#fixnames calde1
clade1_names_fixed <- name_map$name_pop_genome[match(names(clade1_names), name_map$name_Pan_genome)]
#subset tree to just clade1
clade1_tree<- keep.tip(tree_grA_me, clade1_names_fixed)
#fixnames calde2
clade2_names_fixed <- name_map$name_pop_genome[match(names(clade2_names), name_map$name_Pan_genome)]
#subset tree to just clade2
clade2_tree<- keep.tip(tree_grA_me, clade2_names_fixed)
#fixnames clade3
clade3_names_fixed <- name_map$name_pop_genome[match(names(clade3_names), name_map$name_Pan_genome)]
#subset tree to just clade1
clade3_tree<- keep.tip(tree_grA_me, clade3_names_fixed)


#clade 1 plot
clade1_tree_plot <- 
  ggtree(tr = clade1_tree, 
         # color by group attribute, check str(tree_grA_me)
         #mapping = aes(color = group), 
         layout  = 'rectangular') +
         #branch.length = "none") + 
  geom_treescale(x=3, y=NULL, color = NA) +
  scale_color_manual(name = 'Clade', values = my_cols_me) + 
  theme(legend.title=element_text(size=9), # The title of legend 
        legend.text=element_text(size=7)) +
  guides(color = guide_legend(override.aes = list(size = 4))) 

clade1_tree_plot + geom_tiplab(size = 1, align = TRUE, linesize = .25, linetype = 3)

##plot presence / absence of gene fams 

#mutate cols into rows
row.names(exclusive_to_clade1_of_note)<- exclusive_to_clade1_of_note$gene_family
clade_1_binary<- data.frame(t(exclusive_to_clade1_of_note[,5:(ncol(exclusive_to_clade1_of_note) -1)]))
#turn into presence absence 
clade_1_binary_presence<- data.frame(sapply(clade_1_binary, gsub, pattern = "1", replacement = "present"))
clade_1_binary_presence_absence<- sapply(clade_1_binary_presence, gsub, pattern = "0", replacement = "absent")
#View(clade_1_binary_presence_absence)
row.names(clade_1_binary_presence_absence) <- row.names(clade_1_binary)
#fix rownames
clade_1_binary_presence_absence_fixed <- name_map$name_pop_genome[match(rownames(clade_1_binary_presence_absence), name_map$name_Pan_genome)]
row.names(clade_1_binary_presence_absence) <- clade_1_binary_presence_absence_fixed
str(clade_1_binary)
#plot
clade1_tree_plot <-  gheatmap(clade1_tree_plot, 
                              clade_1_binary_presence_absence, 
                              offset=0.01, width=0.6, low="white", high="black", 
                              colnames = F, color="white") +
  scale_fill_manual(values=c("white", "#78B7C5")) +
  ggtitle("Clade 1 specific gene families")


clade1_tree_plot + geom_tiplab(size = .7, align = TRUE, linesize = .25, offset = .65, linetype = 0)
#ggsave(file="Clade_1_specific_families.pdf",device="pdf")



#clade 2
##plot presence / absence of gene fams 
clade2_tree_plot <- 
  ggtree(tr = clade2_tree, 
         # color by group attribute, check str(tree_grA_me)
         #mapping = aes(color = group), 
         layout  = 'rectangular') + 
         #branch.length = "none") + 
  geom_treescale(x=3, y=NULL, color = NA) +
  scale_color_manual(name = 'Clade', values = my_cols_me) + 
  theme(legend.title=element_text(size=9), # The title of legend 
        legend.text=element_text(size=7)) +
  guides(color = guide_legend(override.aes = list(size = 4))) 

#mutate cols into rows
row.names(exclusive_to_clade2_of_note)<- exclusive_to_clade2_of_note$gene_family
clade_2_binary<- data.frame(t(exclusive_to_clade2_of_note[,5:(ncol(exclusive_to_clade2_of_note) -1)]))
#turn into presence absence 
clade_2_binary_presence<- data.frame(sapply(clade_2_binary, gsub, pattern = "1", replacement = "present"))
clade_2_binary_presence_absence<- sapply(clade_2_binary_presence, gsub, pattern = "0", replacement = "absent")
row.names(clade_2_binary_presence_absence) <- row.names(clade_2_binary)
#fix rownames
clade_2_binary_presence_absence_fixed <- name_map$name_pop_genome[match(rownames(clade_2_binary_presence_absence), name_map$name_Pan_genome)]
row.names(clade_2_binary_presence_absence) <- clade_2_binary_presence_absence_fixed


#plot
clade2_tree_plot <-  gheatmap(clade2_tree_plot, 
                              clade_2_binary_presence_absence, 
                              width=1, low="white", high="black", 
                              colnames = F, color="white") +
  scale_fill_manual(values=c("white", "#78B7C5")) +
  ggtitle("Clade 2 specific gene families")

clade2_tree_plot + geom_tiplab(size = 1, align = TRUE, linesize = .25, offset = .6, linetype = 0)
#ggsave(file="Clade_2_specific_families.pdf")


#clade 3
##plot presence / absence of gene fams 
clade3_tree_plot <- 
  ggtree(tr = clade3_tree, 
         # color by group attribute, check str(tree_grA_me)
         #mapping = aes(color = group), 
         layout  = 'rectangular') +
         #branch.length = "none") + 
  #geom_treescale(x=3, y=NULL, color = NA) +
  #scale_color_manual(name = 'Clade', values = my_cols_me) + 
  theme(legend.title=element_text(size=9), # The title of legend 
        legend.text=element_text(size=7))
  #guides(color = guide_legend(override.aes = list(size = 4))) 


#mutate cols into rows
row.names(exclusive_to_clade3_of_note)<- exclusive_to_clade3_of_note$gene_family
clade_3_binary<- data.frame(t(exclusive_to_clade3_of_note[,5:(ncol(exclusive_to_clade3_of_note) -1)]))
#turn into presence absence 
clade_3_binary_presence<- data.frame(sapply(clade_3_binary, gsub, pattern = "1", replacement = "present"))
clade_3_binary_presence_absence<- sapply(clade_3_binary_presence, gsub, pattern = "0", replacement = "absent")
row.names(clade_3_binary_presence_absence) <- row.names(clade_3_binary)
#fix rownames
clade_3_binary_presence_absence_fixed <- name_map$name_pop_genome[match(rownames(clade_3_binary_presence_absence), name_map$name_Pan_genome)]
row.names(clade_3_binary_presence_absence) <- clade_3_binary_presence_absence_fixed


#plot
clade3_tree_plot <-  gheatmap(clade3_tree_plot, 
                              clade_3_binary_presence_absence, 
                              width=1, low="white", high="black", 
                              colnames = F, color="white") +
  scale_fill_manual(values=c("white", "#78B7C5")) +
  ggtitle("Clade 3 specific gene families")

#clade3_tree_plot
clade3_tree_plot + geom_tiplab(size = 1, align = TRUE, linesize = .25, offset = .06, linetype = 0)

#ggsave(file="Clade_3_specific_families.pdf")


##look at gene family losses
#get losses exclusive to clade 1:
exclusive_to_clade1_loss<- clade1[(clade1$totals == 0) & 
                               (clade2$totals > 0) &
                               (clade3$totals > 0),]

dim(exclusive_to_clade1_loss) #there are 16 gene fams lost in Clade 1, that are present in both Clade 2 and Clade 3

#lost in clade 1 and 2, put present in 3
exclusive_to_clade1_and_2<- clade1[(clade1$totals == 0) & 
                                    (clade2$totals == 0) &
                                    (clade3$totals > 0),]

dim(exclusive_to_clade1_and_2) #there are 102 lost in clade 1 and 2 but not 3.

#lost in clade 1 and 3, put present in 2
exclusive_to_clade1_and_3<- clade1[(clade1$totals == 0) & 
                                     (clade2$totals > 0) &
                                     (clade3$totals == 0),]

dim(exclusive_to_clade1_and_3) #there are 131 lost in clade 1 and 3 but not 2.

length(clade1_names)
tail(sort(exclusive_to_clade1_loss$totals)) # top hit is distributed in 201 of the of 201 isolates in clade1

#get all present in more than 1/3 of the isolates 
exclusive_to_clade1_of_note_loss<- exclusive_to_clade1_loss[exclusive_to_clade1_loss$totals > length(clade1_names)*.33,]
dim(exclusive_to_clade1_of_note_loss)
#there are 0 losses exclusive to clade 1 that are in all isolates of clade 1


##Clade2
##look at gene family losses
#get losses exclusive to clade 2:
exclusive_to_clade2_loss<- clade2[(clade1$totals > 0) & 
                                    (clade2$totals == 0) &
                                    (clade3$totals > 0),]

dim(exclusive_to_clade2_loss) #there are 186 gene fams lost in Clade 2, that are present in both Clade 1 and Clade 3


#lost in clade 2 and 1, but present in 3
exclusive_to_clade2_and_1_loss<- clade2[(clade1$totals == 0) & 
                                     (clade2$totals == 0) &
                                     (clade3$totals > 0),]

dim(exclusive_to_clade2_and_1_loss) #there are 102 lost in clade 1 and 2 but not 3.

#lost in clade 2 and 3, but present in 1
exclusive_to_clade2_and_3_loss<- clade2[(clade1$totals > 0) & 
                                     (clade2$totals == 0) &
                                     (clade3$totals == 0),]

dim(exclusive_to_clade2_and_3_loss) #there are 1063 lost in clade 2 and 3 but not 1.

length(clade2_names)
#think about this more...
#tail(sort(exclusive_to_clade2_loss$totals)) # top hit is distributed in 201 of the of 46 isolates in clade2

#get all present in more than 1/3 of the isolates 
exclusive_to_clade2_of_note_loss<- exclusive_to_clade2_loss[exclusive_to_clade2_loss$totals > length(clade2_names)*.33,]
dim(exclusive_to_clade2_of_note_loss)
#there are 0 losses exclusive to clade 1 that are in all isolates of clade 1



##Clade3
##look at gene family losses
#get losses exclusive to clade 3:
exclusive_to_clade3_loss<- clade2[(clade1$totals > 0) & 
                                    (clade2$totals > 0) &
                                    (clade3$totals == 0),]

dim(exclusive_to_clade3_loss) #there are 875 gene fams lost in Clade 3, that are present in both Clade 1 and Clade 2

exclusive_to_clade3_of_note_loss<- exclusive_to_clade3_loss[exclusive_to_clade3_loss$totals > length(clade2_names)*.33,]
















###check quality by BUSCO scores 
# get names of the lowest 5% of accessory genomes
accessory_by_strain_bottom_5percent<- accessory_by_strain[accessory_by_strain$totals < quantile(accessory_by_strain$totals , 0.05 ) , ]
dim(accessory_by_strain_top_5percent)
dim(accessory_by_strain)


#get BUSCO scores
assembly_stats<-as.data.frame(fread("asm_stats.tsv"))

#subset to only problem species 
length(accessory_by_strain_bottom_5percent$strain)

highlight<- c("B9781_CDC-19",
              "F17999",
              "F18454",
              "AZTEC_19_231",    
              "AZTEC_1_235",
              "AZTEC_20_243",
              "AZTEC_16_237",
              "DMC2_AF100-1_20_C",
              "AF100-12_18",
              "AF100-12_21",
              "AF100-12_35",      
              "JCM_10253",
              "M128",
              "08-19-02-46",    
              "117535A-11",
              "C56A1-11")


assembly_stats_subset<- assembly_stats[assembly_stats$SampleID %in% highlight,]
dim(assembly_stats_subset)
length(highlight)

length(assembly_stats_subset$SampleID)

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(assembly_stats$BUSCO_Complete,
  col = ifelse(assembly_stats$SampleID %in% highlight, "#F21A00", "black"),
  ylab = "BUSCO score (%)",
  xlab = "genome number (arbitrary)")
legend("topright", inset=c(-.6,0), legend=c("bottom 5% acessory #"), pch=1,col =c("#F21A00"),pt.cex = 1, cex = .5)


#sort by Busco score 
assembly_stats_sorted<- assembly_stats[order(assembly_stats$BUSCO_Complete),]

#View(assembly_stats_sorted)

###view how depth is influencing coverage and quality
depth_df<-as.data.frame(fread("all_isolates_by_depth.txt"))

#get depth for the isolates in the pan genome study
dim(depth_df)

depth_df_less_than_ten<- depth_df[depth_df$depth<10,]
dim(depth_df_less_than_ten) #there are 40 of them. 
#are any of these genomes in the pan genome? 
depth_df_less_than_ten$strain
#View(depth_df_less_than_ten)

seqs_to_remove<-depth_df_less_than_ten$depth %in% strain_names

test<- intersect(depth_df_less_than_ten$strain, strain_names)

strain_names_match<- sapply(strain_names, gsub, pattern = "Afum_" , replacement = "")

name_map_subset<- name_map[strain_names %in% name_map$name_Pan_genome,]
dim(name_map)
dim(name_map_subset)
depth_df_subset<- depth_df[depth_df$strain %in% strain_names_match,]
dim(depth_df_subset)


#View(strain_names_match)

#print and fix these damn names from the mapping file
length(strain_names)
dim(name_map)

#get list of genomes that are not in the 
strain_names_match<- 


#Print BUSCO and depth scores for supplemental
#BUSCO
assembly_stats<-as.data.frame(fread("asm_stats.tsv"))
dim(assembly_stats)
#depth
depth_df<-as.data.frame(fread("all_isolates_by_depth.txt"))
dim(depth_df)
#strains to subset
tree_grA_me$tip.label


#View(assembly_stats)
#subset
BUSCO_subset<- assembly_stats[assembly_stats$SampleID %in% tree_grA_me$tip.label,]
dim(BUSCO_subset)


#247 - should be 266. Missing some - because some are named differently 
Depth_subset<- depth_df[depth_df$strain %in% tree_grA_me$tip.label,]
dim(Depth_subset)
#258 - should be 266. Missing some - because some are named differently

#write.table(BUSCO_subset, "temp_BUSCO.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(Depth_subset, "temp_Depth.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)





###where do the genes fall on each chr?
#graph where the core genes fall on the chromosome 

#graph where the accessory genes fall on the chromosome

#how many genes fall at chromosome ends? (in the last X BPs). 


#combine accessory and singletons and plot by clade? 

#is genome open or closed?
#colors: "#BA9141" = accessory, "#806633" = core, "#316A6E" = singleton
sp <- specaccum(t(all_accessory), "random", permutations=100)
summary(sp)
plot(sp, ci.type="poly", col="#BA9141", lwd=2, ci.lty=0, ci.col="#BA914199", xlab="Genomes", ylab="Gene Families", ylim = c(0,7000))
#boxplot(sp, col="yellow", add=TRUE)

#plot singletons
#sp2<- specaccum(t(singletons_only), "random", permutations=100)
#plot(sp2, ci.type="poly", col="#316A6E", lwd=2, ci.lty=0, ci.col="#316A6E99", xlab="Genomes", ylab="Gene Families")
#boxplot(sp2, col="yellow", add=TRUE)

#plot core
#core_only<- gene_fam_by_strain_ones_num[rowSums(gene_fam_by_strain_ones_num) > 248,]
#sp3<- specaccum(t(core_only), "random", permutations=100)
#plot(sp3, ci.type="poly", col="#806633", lwd=2, ci.lty=0, ci.col="#80663399", xlab="Genomes", ylab="Gene Families")
#boxplot(sp3, col="yellow", add=TRUE)
par(new=TRUE)
#plot combined  singletons
singletons_and_accessory_only<- gene_fam_by_strain_ones_num[rowSums(gene_fam_by_strain_ones_num) < 248,]
sp4<- specaccum(t(singletons_and_accessory_only), "random", permutations=100)
plot(sp4, ci.type="poly", col="#316A6E", lwd=2, ci.lty=0, ci.col="#316A6E99", xlab="Genomes", ylab="Gene Families", ylim = c(0,7000))
#boxplot(sp2, col="yellow", add=TRUE)
