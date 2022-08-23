#This script determines the optimum number of migration events using bootstrapped replicates of TreeMix 

#######
#install.packages('OptM')
library('OptM')

sessionInfo()
folder<- "~/Desktop/Project_Afum_pangenome_3/TreeMix_results_OG"

test.optM = optM(folder)
test.optM

plot_optM(test.optM, method = "Evanno")
