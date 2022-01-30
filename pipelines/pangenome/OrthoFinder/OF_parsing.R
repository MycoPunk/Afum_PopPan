#A fum pan genome analysis, using proteins clustered by Ortho Finder results
#last updated: 30.Jan.2022

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

setwd("~/Desktop/Project_Afum_pangenome_3/")
OF.gene_counts<-as.data.frame(fread("Orthogroups.GeneCount.tsv")) 
OF.gene_families<-as.data.frame(fread("Orthogroups.tsv")) 
OF.unassigned<- as.data.frame(fread("Orthogroups_UnassignedGenes.tsv")) #this is where your singletons live


#combine OF families with singletons
OF.gene_families_all<- rbind(OF.gene_families, OF.unassigned)

##basic stats:
#how many gene families are there?
n_gene_fams<- nrow(OF.gene_families_all)
n_gene_fams
#15,309 (vs. 15,476 w. PIRATE)

#fix names
names(OF.gene_families_all)<- sapply(names(OF.gene_families_all), gsub, pattern = "Aspergillus_fumigatus_", replacement = "" )
names(OF.gene_families_all)<- sapply(names(OF.gene_families_all), gsub, pattern = ".proteins", replacement = "" )

#how many strains?
#names of cols to exclude
cols_to_exclude<- "Orthogroup"

strains_only<- OF.gene_families_all[,!names(OF.gene_families_all) %in% cols_to_exclude]
strain_names<- names(strains_only)

#length(strain_names)
ngenomes<- length(unique(strain_names)) 
ngenomes
#260

#add number_genomes column to get totals
#first replace blank calls with NAs
strains_only<- apply(strains_only, 2, function(x) gsub("^$|^ $", NA, x))
OF.gene_families_all$number_genomes<- rowSums(!is.na(strains_only)) #cols that are not blank
#replace blanks with NAs for all
OF.gene_families_all<- data.frame(apply(OF.gene_families_all, 2, function(x) gsub("^$|^ $", NA, x)))

#how many of the gene families are in every genome (of 260)
n_gene_fams_core_all<- sum(OF.gene_families_all$number_genomes == 260)
n_gene_fams_core_all
#3,902 (vs. 3,584 in PIRATE)

#that's what percent out of the total?
(n_gene_fams_core_all*100)/n_gene_fams
#25.48827 vs. 23.15844 in PIRATE)

#present in 95% of genomes
cutoff<- round(.95*260)
n_gene_fams_core_w95per<- sum(OF.gene_families_all$number_genomes >= cutoff)
n_gene_fams_core_w95per
#8,866 (vs. 8,600 in PIRATE)

#that's what percent out of the total?
(n_gene_fams_core_w95per*100)/n_gene_fams
#57.91365% (vs. 55.56991% in PIRATE)

#how many of the gene families are singletons (accessory)
n_gene_fams_singletons<- sum(as.numeric(OF.gene_families_all$number_genomes) == 1)
n_gene_fams_singletons
#2,109 (vs. 3,258 in PIRATE)

#that's what percent out of the total?
(n_gene_fams_singletons*100)/n_gene_fams
#13.77621%

#how many accessory?)
n_accessory<- n_gene_fams - (n_gene_fams_singletons + n_gene_fams_core_w95per)
n_accessory
#4,334 accessory (vs. 3,618 in PIRATE)

#that's what percent out of the total?
(n_accessory*100)/n_gene_fams
#28.31014% (vs. 13.91538% in PIRATE)

#get average per genome 
n_accessory / ngenomes

##graph the distribution of gene presence in a gene family (distribution of core to accessory genes)
#plot
gene_fam_totals<-as.data.frame(OF.gene_families_all$number_genomes)
colnames(gene_fam_totals) <- 'count'

#set groups
gene_fam_totals$group = 0                        
for (i in 1:nrow(gene_fam_totals)){
  if (as.numeric(gene_fam_totals$count[i]) == 1) {
    gene_fam_totals$group[i] = "Singleton"
  } else if (gene_fam_totals$count[i] >= cutoff) {
    gene_fam_totals$group[i] = "Core"
  } else {
    gene_fam_totals$group[i] = "Accessory"
  }
}

#fix class
gene_fam_totals$count<- as.numeric(gene_fam_totals$count)


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
#ggsave("dist_by_gene_fam_OF.pdf",p, width=6, height=4, units="in")



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
#ggsave("donut_OF.pdf",p, width=3.5, height=3.5, units="in")
p

##which strain has the highest/lowest number of singletons and accessory gene fams? 
gene_fam_by_strain<-as.data.frame(OF.gene_families_all[,2:261])

##make binary (if gene = 1, if not = 0)
#fill in zeros
#gene_fam_by_strain_zeros<- sapply(gene_fam_by_strain, gsub, pattern = "^\\s*$" , replacement = 0 )

#replace all NAs with 0s
gene_fam_by_strain[is.na(gene_fam_by_strain)] <- 0

#fill in ones
gene_fam_by_strain_ones<- as.data.frame(replace(gene_fam_by_strain, gene_fam_by_strain!="0", 1))
#change to numeric
gene_fam_by_strain_ones_num <- mutate_all(gene_fam_by_strain_ones, function(x) as.numeric(as.character(x)))
#subset to remove core genes from accessory and singletons 
all_accessory_1<- gene_fam_by_strain_ones_num[rowSums(gene_fam_by_strain_ones_num) > 1,]

all_accessory<- all_accessory_1[rowSums(all_accessory_1) < cutoff,]

#subset to get only singletons 
singletons_only<- gene_fam_by_strain_ones_num[rowSums(gene_fam_by_strain_ones_num) == 1,]


#get average accessory per genome
ave_accessory<- colSums(all_accessory)
mean(ave_accessory) #1163.338, vs. 974.3769 with pirate
#variance
accessory_variance<- var(ave_accessory)
accessory_variance #14,087.88
#sd
accesory_sd<- sd(ave_accessory)
accesory_sd #0.4420526

#get average singletons
ave_singletons<- colSums(singletons_only)
mean(ave_singletons) #8.112, vs. 12.5308 with pirate
#variance
singletons_variance<- var(ave_singletons)
singletons_variance #166.9643
#sd
singletons_sd<- sd(ave_singletons)
singletons_sd #12.92147


#fix names (remove X's that were auto-added)
names(all_accessory)<- sapply(names(all_accessory), gsub, pattern = "X", replacement = "" )
names(singletons_only)<- sapply(names(singletons_only), gsub, pattern = "X", replacement = "" )

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
#name_map<-read.delim("clade_map_K3_20Jan2021.txt", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE, check.names = TRUE)
name_map<-read.delim("clade_map_K3_3Jan2022.txt", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE, check.names = TRUE)


row.names(accessory_by_strain) <- name_map$name_to_use_in_paper[match(row.names(accessory_by_strain), name_map$OF_name)]
accessory_by_strain$pop_name<- row.names(accessory_by_strain)
#remove "DMC2" for graphing 
#accessory_by_strain$pop_name<- sapply(accessory_by_strain$pop_name, gsub, pattern = "DMC2_", replacement = "")


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

ggsave("accessory.pdf",p, width=6.9, height=3, units="in")



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
#sort 
singletons_only_by_strain<- singletons_only_by_strain[order(singletons_only_by_strain$totals),]


#fix names
row.names(singletons_only_by_strain) <- name_map$name_to_use_in_paper[match(row.names(singletons_only_by_strain), name_map$OF_name)]
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
#write.table(singletons_only_by_strain, "n_singletons_by_strain_OF.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(accessory_by_strain, "n_accessory_by_strain_OF.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



###graph presence / absence matrix on to big tree

#tree with 260 strains
tree_me <- read.tree("Afum_260_iq_tree_newick_v2.tre")
#remove the reference from the tree:
tree_me<- drop.tip(tree_me, "Af293-REF", trim.internal = TRUE, subtree = FALSE,
                root.edge = 0, rooted = is.rooted(tree_me), collapse.singles = TRUE,
                interactive = FALSE)
#root tree
tree_me<- root(tree_me, node = 505, resolve.root = TRUE,
            interactive = FALSE, edgelabel = FALSE)

#attach the clade annotations 
singletons_only_by_strain$clade<- name_map$clade[match(row.names(singletons_only_by_strain), name_map$name_to_use_in_paper)]
accessory_by_strain$clade<- name_map$clade[match(row.names(accessory_by_strain), name_map$name_to_use_in_paper)]

##use mapping file to rename the tree branches to match the way appear in by singleton and accessory files
#rename tips to represent clades
#function to rename tips
rename.tips <- function(phy, old_names, new_names) {
  mpos <- match(old_names,phy$tip.label)
  phy$tip.label[mpos] <- new_names
  return(phy)
}

tree_me_rename<- rename.tips(tree_me, old_names = name_map$name_pop_genome_new, new_names = name_map$name_to_use_in_paper)
tree_me_rename$tip.label

#get averages by clade
singeltons_clade1<- singletons_only_by_strain[singletons_only_by_strain$clade == "1",]
mean(singeltons_clade1$totals) # 9 (vs. 14 using PIRATE)
singeltons_clade2<- singletons_only_by_strain[singletons_only_by_strain$clade == "2",]
mean(singeltons_clade2$totals) # 4 (vs. 6 using PIRATE)
singeltons_clade3<- singletons_only_by_strain[singletons_only_by_strain$clade == "3",]
mean(singeltons_clade3$totals) # 8 (vs. 16 using PIRATE)

accessory_clade1<- accessory_by_strain[accessory_by_strain$clade == "1",]
mean(accessory_clade1$totals) # 1171 (vs. 982 using PIRATE)
accessory_clade2<- accessory_by_strain[accessory_by_strain$clade == "2",]
mean(accessory_clade2$totals) # 1160 (vs. 982 using PIRATE)
accessory_clade3<- accessory_by_strain[accessory_by_strain$clade == "3",]
mean(accessory_clade3$totals) # 1075 (vs. 974 using PIRATE)

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
#split by clade
grA_me<- split(row.names(accessory_by_strain), accessory_by_strain$clade)

#split by mat type
grA_me_mat<- split(name_map$name_to_use_in_paper, name_map$MAT_type)

#set colors by group for clade
tree_grA_me <- ggtree::groupOTU(tree_me_rename, grA_me)
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

names(my_cols_me) <- levels(attributes(tree_grA_me)$group)
scales::show_col(my_cols_me); my_cols_me


###set colors by group for mat type 
tree_grA_me_mat <- ggtree::groupOTU(tree_me_rename, grA_me_mat)
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
  

#simple plot
tree_plot_me <- 
  ggtree(tr = tree_grA_me, 
         # color by group attribute, check str(tree_grA_me)
         mapping = aes(color = group), 
         layout  = 'circular', 
        # layout  = 'rectangular', 
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

####map matting type
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

name_map_sm<- name_map[,c("name_to_use_in_paper", "MAT_type")]

#combine annotations for clade and MAT type into one tree
pi<- tree_plot_me + geom_tiplab(size = .6, align = TRUE, linesize = .20, linetype = NA) 
p <- pi %<+% name_map_sm + geom_tippoint(aes(x=x+14,color=MAT_type), size=.30) + scale_color_manual(values=c("#56326E", "#ED7F6F","#ABA778", "#BEBDBD", "#404041", "#D4494E"))
plot(p)


#theme(legend.title=element_text(size=10), # The title of legend 
#      legend.text=element_text(size=7))+
#guides(color = guide_legend(override.aes = list(linetype = c(1, 1, 1, 0, 0, 0), size = c(1,1,1,2,2,2), shape = c(NA, NA,NA, 16, 16, 16))))

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

#ggsave(file="pan_genome_tree_w_bars.pdf",device="pdf")

#add MAT1-2-4 to this tree?
OF.gene_families<-as.data.frame(fread("Orthogroups.tsv")) 
OF.unassigned<- as.data.frame(fread("Orthogroups_UnassignedGenes.tsv")) #this is where your singletons live
OFtoAfu<- as.data.frame(fread("OF_blastP_results_filtered.txt", header = FALSE)) #this is the OG to Afu designations
colnames(OFtoAfu)<- c("Afu_name", "OG_name")
#combine OF families with singletons
OF.gene_families_all<- rbind(OF.gene_families, OF.unassigned)
#fix names
names(OF.gene_families_all)<- sapply(names(OF.gene_families_all), gsub, pattern = "Aspergillus_fumigatus_", replacement = "" )
names(OF.gene_families_all)<- sapply(names(OF.gene_families_all), gsub, pattern = ".proteins", replacement = "" )
#names of cols to exclude
cols_to_exclude<- "Orthogroup"
strains_only<- OF.gene_families_all[,!names(OF.gene_families_all) %in% cols_to_exclude]
strain_names<- names(strains_only)
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
#subset large df to only Afu3g06160 (MAT1-2-4)
MAT1.2.4<- gene_fam_by_strain_w_anno_ones_num2[gene_fam_by_strain_w_anno_ones_num2$Af293_reference_gene_name %in% "Afu3g06160",]
#subset to only counts
all_vars<- MAT1.2.4[,2:261]
#change to presence / absence
all_vars_presence<- data.frame(lapply(all_vars, gsub, pattern = "1", replacement = "present"))
all_vars_presence_absence<- data.frame(lapply(all_vars_presence, gsub, pattern = "0", replacement = "absent"))
#fix col names X's introduced in last step because R likes to make life difficult
colnames(all_vars_presence_absence)<- colnames(all_vars)
#add rowames
rownames(all_vars_presence_absence)<- MAT1.2.4$Af293_reference_gene_name
#transpose for graphing
all_vars_presence_absence_t<- data.frame(t(all_vars_presence_absence))
#fix names to match OG names
row.names(all_vars_presence_absence_t) <- name_map$name_to_use_in_paper[match(row.names(all_vars_presence_absence_t), name_map$OF_name)]

#get counts
rownames(all_vars)<- MAT1.2.4$Af293_reference_gene_name
rowSums(all_vars) #there are 116 strains w/ MAT-1-2-4


#plot tree with all annotation
tree_plot <-  gheatmap(tree_bars_me, 
                       all_vars_presence_absence_t,
                       offset=15, width=0.02, low="white", high="black", 
                       colnames = F, 
                       #colnames_angle = 45,
                       #colnames_position = "top",
                       #font.size = 1,
                       #colnames_offset_y = -6,
                       #colnames_offset_x = .05,
                       color="white") +
  scale_fill_manual(values=c("white", "black")) 

tree_all_data<- tree_plot + geom_tiplab(size = .09, align = TRUE, linesize = .20, linetype = NA)+ guides(fill = guide_legend(ncol = 3))+
  theme(legend.title=element_text(size=10), # The title of legend 
        legend.text=element_text(size=7))+
  guides(color = guide_legend(override.aes = list(linetype = c(1, 1, 1, 0, 0, 0), size = c(1,1,1,2,2,2), shape = c(NA, NA,NA, 16, 16, 16)), title = "", nrow = 3)) 

tree_all_data 

ggsave("pan_genome_tree_w_bars_and_mat1-2-4.pdf", tree_all_data, width=22, height=22, 
       device = "pdf", units = "cm")


sum(name_map$MAT_type == "MAT1") #144
sum(name_map$MAT_type == "MAT2") #105
sum(name_map$MAT_type == "unknown") #11
######




###are there unique gene families (or unique losses) by clade?

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



#subset to remove core genes and singletons 
#remove core
all_accessory_anno<- gene_fam_by_strain_w_anno_ones_num2[rowSums(gene_fam_by_strain_w_anno_ones_num2[,2:ncol(gene_fam_by_strain_w_anno_ones_num2)]) < cutoff,]
dim(all_accessory_anno)

#remove singletons
all_accessory_anno_no_sing<- all_accessory_anno[rowSums(all_accessory_anno[,2:ncol(all_accessory_anno)]) != 1,]
dim(all_accessory_anno_no_sing)

#split list of strain names by clade
accessory_by_strain_list <- setNames(accessory_by_strain$clade, accessory_by_strain$strain)

##split the count data frame by the factor list
clade1_names<- accessory_by_strain_list[accessory_by_strain_list == "1"]
clade1<- cbind(OF_gene_fam = all_accessory_anno_no_sing[,1], all_accessory_anno_no_sing[,colnames(all_accessory_anno_no_sing) %in% names(clade1_names)])

clade2_names<- accessory_by_strain_list[accessory_by_strain_list == "2"]
clade2<- cbind(OF_gene_fam = all_accessory_anno_no_sing[,1], all_accessory_anno_no_sing[,colnames(all_accessory_anno_no_sing) %in% names(clade2_names)])

clade3_names<- accessory_by_strain_list[accessory_by_strain_list == "3"]
clade3<- cbind(OF_gene_fam = all_accessory_anno_no_sing[,1], all_accessory_anno_no_sing[,colnames(all_accessory_anno_no_sing) %in% names(clade3_names)])

dim(clade1) #note this is +1 in ea row as the OG is the first col
dim(clade2)
dim(clade3)

#are there gene families that are exclusive to clade 1?
clade1$totals<-  rowSums(clade1[,2:ncol(clade1)])
clade2$totals<-  rowSums(clade2[,2:ncol(clade2)])
clade3$totals<-  rowSums(clade3[,2:ncol(clade3)])


#get fams exclusive to clade1
exclusive_to_clade1<- clade1[(clade1$totals > 0) & 
                               (clade2$totals == 0) &
                               (clade3$totals == 0),]

nrow(exclusive_to_clade1)

aveclade1<- nrow(exclusive_to_clade1) / length(clade1_names) #there are 1256 accessory gene fams exclusive to clade 1 (1062  using PIRATE)
round(aveclade1) #ave is 6

tail(sort(exclusive_to_clade1$totals)) # most abundant gene fam is present in in 142 of the 200 isolates in clade1


length(clade1_names)*.90 
#get all present in more than 90% of the isolates in that clade
exclusive_to_clade1_of_note<- exclusive_to_clade1[exclusive_to_clade1$totals > length(clade1_names)*.90,]
dim(exclusive_to_clade1_of_note)
#there are no gene fams in =>95% of the strians in clade 1
exclusive_to_clade1_of_note

#exclusive to clade2
exclusive_to_clade2<- clade2[(clade2$totals > 0) & 
                               (clade1$totals == 0) &
                               (clade3$totals == 0),]

dim(exclusive_to_clade2) #there are 95 gene fams exclusive to clade2
aveclade2<- nrow(exclusive_to_clade2) / length(clade2_names)
round(aveclade2) #ave is 2
length(clade2_names) #there are 45 isolates in clade 2
tail(sort(exclusive_to_clade2$totals)) #most abundant appear in 42 of the of the 45 isolates (two gene fams)
exclusive_to_clade2_of_note<- exclusive_to_clade2[exclusive_to_clade2$totals > length(clade2_names)*.90,]
dim(exclusive_to_clade2_of_note) #2 clade-defining gains in Clade 2: OG0010568 & OG0010571


#exclusive to clade3
exclusive_to_clade3<- clade3[(clade3$totals > 0) & 
                               (clade1$totals == 0) &
                               (clade2$totals == 0),]


dim(exclusive_to_clade3) 
aveclade3<- nrow(exclusive_to_clade3) / length(clade3_names)
round(aveclade3) #ave is 8

length(clade3_names)
tail(sort(exclusive_to_clade3$totals))
#there are 115 accessory gene fams exclusive to clade3, some in all 15 isolates
exclusive_to_clade3_of_note<- exclusive_to_clade3[exclusive_to_clade3$totals > length(clade3_names)*.90,]
dim(exclusive_to_clade3_of_note) #23 are present in more than > 90% of all isolates

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



#attach annotations to clade exclusive gene fams and print
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

#rename cols so that they match 
colnames(anno_all)[1] <- "OF_gene_fam"


#add anno col
exclusive_to_clade1$interpro_annotation<- anno_all$InterPro_annotations_description[match(exclusive_to_clade1$OF_gene_fam, anno_all$OF_gene_fam)]
exclusive_to_clade2$interpro_annotation<- anno_all$InterPro_annotations_description[match(exclusive_to_clade2$OF_gene_fam, anno_all$OF_gene_fam)]
exclusive_to_clade3$interpro_annotation<- anno_all$InterPro_annotations_description[match(exclusive_to_clade3$OF_gene_fam, anno_all$OF_gene_fam)]

exclusive_to_clade1$GO_annotation<- anno_all$GO_annotations[match(exclusive_to_clade1$OF_gene_fam, anno_all$OF_gene_fam)]
exclusive_to_clade2$GO_annotation<- anno_all$GO_annotations[match(exclusive_to_clade2$OF_gene_fam, anno_all$OF_gene_fam)]
exclusive_to_clade3$GO_annotation<- anno_all$GO_annotations[match(exclusive_to_clade3$OF_gene_fam, anno_all$OF_gene_fam)]


#if cell is empty, replace with NA
exclusive_to_clade1<- as.data.frame(apply(exclusive_to_clade1, 2, function(x) gsub("^$|^ $|^-$", NA, x)))
exclusive_to_clade2<- as.data.frame(apply(exclusive_to_clade2, 2, function(x) gsub("^$|^ $|^-$", NA, x)))
exclusive_to_clade3<- as.data.frame(apply(exclusive_to_clade3, 2, function(x) gsub("^$|^ $|^-$", NA, x)))

##validate - OG exclusive to C2 and C3 should primarily NOT have Afu blast hits if they are real (because Af293 is in C1)
#load data from OG to Afu BLAST search
OFtoAfu<- as.data.frame(fread("OF_blastP_results_filtered.txt", header = FALSE)) #this is the OG to Afu designations

#assign col names
fmt6_names<- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
colnames(OFtoAfu)<- fmt6_names

#intersections between these two sets?
intersect(exclusive_to_clade2$OF_gene_fam, OFtoAfu$sseqid) #none
length(exclusive_to_clade2$OF_gene_fam) #none out of 95
intersect(exclusive_to_clade3$OF_gene_fam, OFtoAfu$sseqid) #2 "OG0011685" "OG0013027"
length(exclusive_to_clade3$OF_gene_fam) #two out of 115 in that set



#Attach Afu annotations for Clade 1
exclusive_to_clade1$Af293_reference_gene_name<- OFtoAfu$qseqid[match(exclusive_to_clade1$OF_gene_fam, OFtoAfu$sseqid)]

#how many clade specific gene families in C1 are annotated with an Af293 gene name?
exclusive_to_clade1[exclusive_to_clade1=='NA'] <- NA
how_many <- exclusive_to_clade1[!is.na(exclusive_to_clade1$Af293_reference_gene_name),]
dim(how_many) #28


#Print all exclusive gene fams for sup table
#write.table(exclusive_to_clade1, "exclusive_to_clade1_all.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(exclusive_to_clade2, "exclusive_to_clade2_all.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(exclusive_to_clade3, "exclusive_to_clade3_all.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


##do the same for losses
#get fams lost exclusive to clade1
exclusive_to_clade1_lost<- clade1[(clade1$totals == 0) & 
                                    (clade2$totals > 0) &
                                    (clade3$totals > 0),]

dim(exclusive_to_clade1_lost) #there are 24 accessory gene fams exclusively lost in Clade 1

aveclade1_lost<- nrow(exclusive_to_clade1_lost) / length(clade1_names)
round(aveclade1_lost) #that's well less than one per genome. 
aveclade1_lost #0.12 to be exact

#all have been lost in all isolates of clade 1 - per the previous search perams. 
#how abundant are they in the other clades though? 
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

dim(all_totals_lost_in_clade1) #24 of these
all_totals_lost_in_clade1_of_interest<- all_totals_lost_in_clade1[(all_totals_lost_in_clade1$percent_abundance_in_Clade2 > 90) |
                                                                    (all_totals_lost_in_clade1$percent_abundance_in_Clade3 > 90),]

dim(all_totals_lost_in_clade1_of_interest) #5 clade defining

#get the abundance in clades 2 and 3 for the genes identified in exclusive_to_clade1_lost$gene_family

#get gene fams lost exclusive to clade2
exclusive_to_clade2_lost<- clade2[(clade2$totals == 0) & 
                                    (clade1$totals > 0) &
                                    (clade3$totals > 0),]

dim(exclusive_to_clade2_lost) # there are 270 gene fams lost in Clade 2 that are present in Clades 1 and 3
aveclade2_lost<- nrow(exclusive_to_clade2_lost) / length(clade2_names)
round(aveclade2_lost) #that's an average of 6 genes per strain


#subset to losses in clade2 - total losses
all_totals_lost_in_clade2<- all_totals[(all_totals$clade1_totals > 0) & 
                                         (all_totals$clade2_totals == 0) &
                                         (all_totals$clade3_totals > 0),]

dim(all_totals_lost_in_clade2) #270 lost exclusively in Clade2


all_totals_lost_in_clade2_of_interest<- all_totals_lost_in_clade2[(all_totals_lost_in_clade2$percent_abundance_in_Clade1 > 90) |
                                                                    (all_totals_lost_in_clade2$percent_abundance_in_Clade3 >90),]

dim(all_totals_lost_in_clade2_of_interest) #25 are Clade defining losses 


#get gene fams lost exclusive to clade3
exclusive_to_clade3_lost<- clade3[(clade3$totals == 0) & 
                                    (clade1$totals > 0) &
                                    (clade2$totals > 0),]

dim(exclusive_to_clade3_lost) #991 lost in Clade 3 but present in Clades 1 and 2. 
aveclade3_lost<- nrow(exclusive_to_clade3_lost) / length(clade3_names)
round(aveclade3_lost) #that's an average of 66 genes per strain

#subset to losses in clade3
all_totals_lost_in_clade3<- all_totals[(all_totals$clade1_totals > 0) & 
                                         (all_totals$clade2_totals > 0) &
                                         (all_totals$clade3_totals == 0),]

all_totals_lost_in_clade3_of_interest<- all_totals_lost_in_clade3[(all_totals_lost_in_clade3$percent_abundance_in_Clade1 > 90) |
                                                                    (all_totals_lost_in_clade3$percent_abundance_in_Clade2 > 90),]

dim(all_totals_lost_in_clade3_of_interest) #125 clade defining losses


###
#bind tables of interest, and attach annotations via. Jason's output
#first add row names
all_totals_lost_in_clade1_of_interest$gene_fam<- rownames(all_totals_lost_in_clade1_of_interest)
all_totals_lost_in_clade1_of_interest$greater_than_90<- (all_totals_lost_in_clade1_of_interest$percent_abundance_in_Clade2 > 90) | (all_totals_lost_in_clade1_of_interest$percent_abundance_in_Clade3 > 90)

all_totals_lost_in_clade2_of_interest$gene_fam<- rownames(all_totals_lost_in_clade2_of_interest)
all_totals_lost_in_clade2_of_interest$greater_than_90<- (all_totals_lost_in_clade2_of_interest$percent_abundance_in_Clade1 > 90) | (all_totals_lost_in_clade2_of_interest$percent_abundance_in_Clade3 > 90)

all_totals_lost_in_clade3_of_interest$gene_fam<- rownames(all_totals_lost_in_clade3_of_interest)
all_totals_lost_in_clade3_of_interest$greater_than_90<- (all_totals_lost_in_clade3_of_interest$percent_abundance_in_Clade2 > 90) | (all_totals_lost_in_clade3_of_interest$percent_abundance_in_Clade1 > 90)



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
lostinClade1.1<- data.frame(cbind(value = lostinClade1.0$value, group =lostinClade1.0$name, gene_fam=lostinClade1.0$gene_fam, greater_than_90 = lostinClade1.0$greater_than_90))
lostinClade1.1$value<- as.numeric(lostinClade1.1$value)
lostinClade1.1$group<- as.factor(lostinClade1.1$group)
lostinClade1.1$gene_fam<- as.factor(lostinClade1.1$gene_fam)
lostinClade1.1$greater_than_90<- as.factor(lostinClade1.1$greater_than_90)
lostinClade1.1$lost_in<- as.factor("lost_in_Clade1")

lostinClade2.1<- data.frame(cbind(value = lostinClade2.0$value, group =lostinClade2.0$name, gene_fam=lostinClade2.0$gene_fam, greater_than_90 = lostinClade2.0$greater_than_90))
lostinClade2.1$value<- as.numeric(lostinClade2.1$value)
lostinClade2.1$group<- as.factor(lostinClade2.1$group)
lostinClade2.1$gene_fam<- as.factor(lostinClade2.1$gene_fam)
lostinClade2.1$greater_than_90<- as.factor(lostinClade2.1$greater_than_90)
lostinClade2.1$lost_in<- as.factor("lost_in_Clade2")

lostinClade3.1<- data.frame(cbind(value = lostinClade3.0$value, group =lostinClade3.0$name, gene_fam=lostinClade3.0$gene_fam), greater_than_90 = lostinClade3.0$greater_than_90)
lostinClade3.1$value<- as.numeric(lostinClade3.1$value)
lostinClade3.1$group<- as.factor(lostinClade3.1$group)
lostinClade3.1$gene_fam<- as.factor(lostinClade3.1$gene_fam)
lostinClade3.1$greater_than_90<- as.factor(lostinClade3.1$greater_than_90)
lostinClade3.1$lost_in<- as.factor("lost_in_Clade3")

#get total of both above %90
sum(lostinClade1.1$greater_than_90 =="TRUE") /2 #5
sum(lostinClade2.1$greater_than_90 =="TRUE") /2 #25
sum(lostinClade3.1$greater_than_90 =="TRUE") /2 #125


##ggplot
##ggplot way to color lines
clade1_loss<- ggplot(lostinClade1.1, aes(x = group, y = value)) +
  #geom_boxplot(aes(fill = group), alpha = 0.2, col = "grey") +
  geom_line(aes(group = gene_fam, col = greater_than_90), alpha=0.4) +
  scale_colour_manual(values = c("#ED7F6F", "#ABA778","grey"))+
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
  geom_line(aes(group = gene_fam, col = greater_than_90), alpha=0.4) +
  scale_colour_manual(values = c("#56326E", "#ABA778", "grey"))+
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
  geom_line(aes(group = gene_fam, col = greater_than_90), alpha=0.4) +
  #scale_colour_manual(values = c("grey", "#56326E", "#ED7F6F", "#E1AF00"))+
  scale_colour_manual(values = c("#56326E", "#ED7F6F", "grey"))+
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


###get MAT-type stats
Mat_df<- data.frame(table(name_map$clade, name_map$MAT_type))

Mat_df_wide<- pivot_wider(Mat_df, id_cols = Var1)
Mat_df_wide

Mat_df_wide<- Mat_df %>%
  pivot_wider(names_from = Var2, values_from = Freq, values_fill = 0)

C1_mat_ratio<- reduce.fraction(c(Mat_df_wide$MAT1[1], Mat_df_wide$MAT2[1]))
C2_mat_ratio<-reduce.fraction(c(Mat_df_wide$MAT1[2], Mat_df_wide$MAT2[2]))
C3_mat_ratio<-reduce.fraction(c(Mat_df_wide$MAT1[3], Mat_df_wide$MAT2[3]))
C1_mat_ratio
C2_mat_ratio
C3_mat_ratio

#total MAT1
total_MAT1<- sum(as.numeric(Mat_df_wide$MAT1))
total_MAT1

#total MAT2
total_MAT2<-sum(as.numeric(Mat_df_wide$MAT2))
total_MAT2

#ratio
all_mat_ratio<- reduce.fraction(c(total_MAT1, total_MAT2))
all_mat_ratio #can't be reduced 

