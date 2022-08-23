#generates a TreeMix plot for the A.fumigatus pan genome dataset

#######

sessionInfo()
#set wd
setwd("~/Desktop/Project_Afum_pangenome_3/TreeMix_w_OG_optimized")

#plot
source("~/Desktop/R/treemix-1.13/src/plotting_funcs.R")
plot_tree(stem="TreeMix")
