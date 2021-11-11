#this scrip looks at the abundance of CAZYs genes in the A.fumigatus pan genome. 

#load modules
library(phylosignal)
library(adephylo)
library(phylobase)
library(phytools)
library(caper)
library(geiger)
library(pheatmap)
library(tidyr)
library(dplyr)
library(phangorn) 
library(tidyverse)
library(RColorBrewer)
library(viridis)
require(data.table)
library(FSA)
library(rcompanion)

library(dplyr) #note- this version causes a conflict error with ggtree., down grade to v1.0.5. 
packageVersion("dplyr") #‘1.0.7’
#remove.packages("dplyr")
#devtools::install_github("hadley/dplyr@v1.0.5")
#remotes::install_github("YuLab-SMU/ggtree", force = TRUE)
#remotes::install_github("YuLab-SMU/tidytree")
library(ggtree)
library(tidyverse)
library(ggplot2)
library(ape)
library(car)
library(caret)
library(resample)


#set seed for reproducibility
set.seed(666)

setwd("/CAZY")

#read in the files
all_cazy_list<- list.files(path ="~/CAZY", pattern = "*.CAZY.txt", full.names=T)
all_cazy<-lapply(all_cazy_list, fread, sep = c("\t"), header = FALSE, colClasses = 'character', stringsAsFactors=FALSE)

##clean up the files to only relivent info 
#remove "CAZy:" prefix
all_cazy[]<-lapply(rapply(all_cazy, function(x) 
  gsub("CAZy:", "", x), how = "list"), 
  as.data.frame)

#name the cols
col_names<-c("junk1", "junk2", "CAZy","strain")
all_cazy <-lapply(all_cazy, setNames, col_names)

#remove the first two cols
all_cazy_clean<- lapply(all_cazy, function(x) x[!(names(x) %in% c("junk1", "junk2"))])

#clean up
rm(all_cazy)

#bind all dfs together 
big_df<-bind_rows(all_cazy_clean)

#add "Afum_" prefix to strain names
big_df$strain = paste0('Afum_', big_df$strain)

#call table
big_df_table<- data.frame(table(big_df))

#clean up
rm(all_cazy_clean)
rm(big_df)

#make sure you got them all - there should be 260
length(unique(big_df_table$strain)) #good
unique(big_df_table$strain)
#how many cazy fams were identified?
length(unique(big_df_table$CAZy))
#basic stats 
stats_df<-data.frame(totals_CAZy =rowsum(big_df_table$Freq, group = big_df_table$strain))
mean(stats_df$totals_CAZy) #460.1038
sd(stats_df$totals_CAZy) #7.675412
min(stats_df$totals_CAZy) #422
max(stats_df$totals_CAZy) #478


##split dfs by clade association 
#read in clade associations
name_map<-read.delim("~/Desktop/Project_Afum_pangenome_2/clade_map_K3_20Jan2021.txt", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE)

#check that you got them all, and not names are messed up - there should be 260
setdiff(big_df_table$strain,name_map$name_Pan_genome) #good
setdiff(name_map$name_Pan_genome, big_df_table$strain) #good

clade_1_names<- name_map[name_map$clade == "1",]
clade_2_names<- name_map[name_map$clade == "2",]
clade_3_names<- name_map[name_map$clade == "3",]

clade_1<- big_df_table[big_df_table$strain %in% clade_1_names$name_Pan_genome,]
clade_2<- big_df_table[big_df_table$strain %in% clade_2_names$name_Pan_genome,]
clade_3<- big_df_table[big_df_table$strain %in% clade_3_names$name_Pan_genome,]

#read in tree
tree_me <- read.tree("~/Desktop/Project_Afum_pangenome_2/Afum_260_iq_tree_newick.tre")


#remove the reference from the tree:
tree_me<- drop.tip(tree_me, "Af293-REF", trim.internal = TRUE, subtree = FALSE,
                   root.edge = 0, rooted = is.rooted(tree_me), collapse.singles = TRUE,
                   interactive = FALSE)


#match tree (for the couple sp. that are miss-named)
tree_me$tip.label[tree_me$tip.label=="F18149"] <- "F18149-JCVI"
tree_me$tip.label[tree_me$tip.label=="Afu_343-P/11"] <- "Afu_343-P-11"
tree_me$tip.label[tree_me$tip.label=="AFIS1435CDC_6"] <- "AFIS1435_CDC_6"

#root tree based on small outgroup tree (at node 266)
tree_me<- root(tree_me, node = 266, resolve.root = FALSE,
               interactive = FALSE, edgelabel = FALSE)

#split by clade
grA_me<- split(name_map$name_pop_genome, name_map$clade)


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
names(my_cols_me) <- levels(attributes(tree_grA_me)$group)
scales::show_col(my_cols_me); my_cols_me
#simple plot
tree_plot_me <- 
  ggtree(tr = tree_grA_me, 
         # color by group attribute, check str(tree_grA_me)
         mapping = aes(color = group), 
         #layout  = 'circular', 
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


# plot and add the tip labels
tree_plot_me + geom_tiplab(size = 1, align = TRUE, linesize = .1, linetype = 0)

#map the CAZYs to the tree
#format data

CAZYs_to_map<- data.frame(rbind(clade_1, clade_2, clade_3))
CAZYs_to_map_wide<- pivot_wider(CAZYs_to_map, id_cols = strain, values_from = Freq, names_from = CAZy)
#replace strain names with mapping file 

merged_temp<- data.frame(merge(CAZYs_to_map_wide, name_map, by.x = "strain", by.y = "name_Pan_genome"))
rownames(merged_temp)<- merged_temp$name_pop_genome

#drop cols you don't need
#CAZYs_to_map_clean<- subset(merged_temp, select=-c(strain,name_pop_genome, clade))
CAZYs_to_map_clean<- subset(merged_temp, select=-c(strain,name_pop_genome, OF_name, MAT_type))


#### tree to plot on
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

#map and rotate the nodes to be in order of the clades 
#tree_to_map_on<- flip(tree_plot_me, 387, 264) %>% flip(309, 295) + geom_tiplab(size = 0, align = TRUE, linesize = .25, linetype = 3)
#tree_to_map_on
####
#asthetics
tree_plot_me<- tree_plot_me +geom_tiplab(size = 0, align = TRUE, linesize = .25, linetype = 3)
tree_plot_me


#plot
CAZY_tree_plot <-  gheatmap(tree_plot_me, 
                            CAZYs_to_map_clean[1:(length(CAZYs_to_map_clean)-1)],
                            low = "white",
                            high = "black",
                              offset=0.02, width=2, 
                              colnames = T, 
                              colnames_angle = 45,
                              colnames_position = "top",
                              font.size = 1,
                              colnames_offset_y = -6,
                              colnames_offset_x = .05,
                              color="white")
  #scale_fill_viridis_c(option="A", name="continuous\nvalue") 

p<- CAZY_tree_plot + geom_tiplab(size = .8, align = TRUE, linesize = .25, offset = 2.75, linetype = 0)
p


#remove cols with zero variance 
colVars(CAZYs_to_map_clean)


#which cols have 0 variance?
which(apply(CAZYs_to_map_clean, 2, var) == 0)

#remove cols with 0 variance
CAZYs_to_map_clean_no_zero <- CAZYs_to_map_clean[ - as.numeric(which(apply(CAZYs_to_map_clean, 2, var) == 0))]
ncol(CAZYs_to_map_clean)
ncol(CAZYs_to_map_clean_no_zero)
colVars(CAZYs_to_map_clean_no_zero) #looks good

#####
##check variance and normalcy assumptions 
#leveneTest from the car package, look for homogenity of variance.
x <- "clade"
CAZYs_to_map_clean_no_zero<- CAZYs_to_map_clean_no_zero[c(x, setdiff(names(CAZYs_to_map_clean_no_zero), x))]

levene_results<- lapply(subset(CAZYs_to_map_clean_no_zero, select = -clade), leveneTest, group = as.factor(CAZYs_to_map_clean_no_zero$clade))
#extract p-vals
extract1<- lapply(levene_results, `[[`, c("Pr(>F)")) 
#pval- less than 0.05? if so, reject= variance is not equal. 
extract2<- data.frame(sapply(extract1, "[[", 1) < 0.05)
table(extract2) #26 are less than 0.05, so not all categories have equal variance. = need to use non parametric test.


#test normalcy of residuals with shapiro.test (if < 0.05, data is normal)
shapiro_results<- apply(CAZYs_to_map_clean_no_zero[,2:ncol(CAZYs_to_map_clean_no_zero)],2,shapiro.test)
#get p-values
shapiro_results2 <- sapply(shapiro_results, `[`, c("statistic","p.value"))
shapiro_results3<- t(shapiro_results2)
shapiro_results3 #data is not normal - need to run non-parametric test
#####

##run non-parametric (kruskal walliace) test as a loop on each column. 
#fix data type
CAZYs_to_map_clean_no_zero<- as.data.frame(CAZYs_to_map_clean_no_zero)
#loop
all_kw_p_vals<- data.frame(p.val =apply(CAZYs_to_map_clean_no_zero[,-1], 2, function(x) kruskal.test(x,CAZYs_to_map_clean_no_zero[,1])$p.value))
#bonf adjustment for multiple comp. 
all_kw_p_vals$p.adjust<- p.adjust(all_kw_p_vals$p.val, method = "bonferroni", 
                               n = nrow(all_kw_p_vals))

p_vals_kw_sig<- all_kw_p_vals[all_kw_p_vals$p.adjust < 0.001,]


######
names_of_sig_results<- rownames(p_vals_kw_sig)


######

##conduct posthoc tests. 

#subset only significant results
df.subset_and_clade <- CAZYs_to_map_clean_no_zero[, c("clade",names_of_sig_results)]

#as loop
DT_loop_vals<- apply(df.subset_and_clade[,-1], 2, function(x) dunnTest(x,df.subset_and_clade[,1], method="bonferroni"))
length(DT_loop_vals)

### Compact letter display loop
#instantiate empty data frame
CAZy_table<- data.frame(matrix(NA, nrow = 3, ncol = 0))
                        
for (i in seq_along(DT_loop_vals)){
  
  name <- (names(DT_loop_vals)[i])
  PT = DT_loop_vals[[i]]$res
  x<- cldList(comparison = PT$Comparison,
              p.value    = PT$P.adj,
              threshold  = 0.05)
  row_to_add<-x$Letter
  CAZy_table[[name]] <- row_to_add
}
CAZy_table


#######
#make means table 
C1_mean_temp2<- df.subset_and_clade[df.subset_and_clade$clade == 1,]
C2_mean_temp2<- df.subset_and_clade[df.subset_and_clade$clade == 2,]
C3_mean_temp2<- df.subset_and_clade[df.subset_and_clade$clade == 3,]


clade1_means<- round(colMeans(C1_mean_temp2),1)
clade2_means<- round(colMeans(C2_mean_temp2),1)
clade3_means<- round(colMeans(C3_mean_temp2),1)

df_of_CAZY_means<- as.data.frame(rbind("Clade 1" = clade1_means, "Clade 2" = clade2_means, "Clade 3" = clade3_means))
df_of_CAZY_means<- df_of_CAZY_means[-c(1)]
rownames(CAZy_table)<- c("Clade 1", "Clade 2", "Clade 3")

#######

#combine post-hoc and mean tables for printing 
mypaste <- function(x,y) paste(x, " (", y, ")", sep="")
table_for_print<- data.frame(mapply(mypaste, df_of_CAZY_means, CAZy_table))
rownames(table_for_print)<- c("Clade 1", "Clade 2", "Clade 3")

#write.table(table_for_print, "CAZyme_table.txt", sep="\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

####### 



#get names of significant results. 
names_of_sig_results<- rownames(p_vals_kw_sig)
#subset from the main df
#sig_df<- CAZYs_to_map_clean_no_zero[names(CAZYs_to_map_clean_no_zero) %in% names_of_sig_results,]
df.subset <- CAZYs_to_map_clean_no_zero[, names_of_sig_results]

#Scale
test<- apply(df.subset, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))

#map these onto tree. 
CAZY_tree_plot <-  gheatmap(tree_plot_me, 
                            test,
                            #low = "white",
                            #high = "black",
                            offset=0.02, width=2, 
                            colnames = T, 
                            colnames_angle = 45,
                            colnames_position = "top",
                            font.size = 2,
                            colnames_offset_y = -6,
                            colnames_offset_x = .05,
                            color="white") +
#scale_fill_viridis_c(option="F", name="continuous\nvalue") 
scale_fill_viridis_c(option="F", name="normalized\nabundance") 
p<- CAZY_tree_plot + geom_tiplab(size = .8, align = TRUE, linesize = .25, offset = 2.75, linetype = 0)
p
ggsave(file="CAZY_tree.pdf",device="pdf", p, width=8, height=8, units="in")


for_print<- data.frame(names_of_sig_results)
#write.table(for_print, "CAZy_anno.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#big supplemental table with all counts
#fix names
good_names <- name_map$name_pop_genome[match(CAZYs_to_map_wide$strain, name_map$name_Pan_genome)]
CAZYs_to_map_wide_fixed<- cbind(good_names, CAZYs_to_map_wide)
#write.table(CAZYs_to_map_wide_fixed, "Afum_Table_S4.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
