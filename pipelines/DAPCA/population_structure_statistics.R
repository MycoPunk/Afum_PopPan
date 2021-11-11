##this script runs population stats on a vcf file and population defined previously using DAPCA
##Last updated 10.Nov.2021

##set wd
setwd("")

##load libraries
library("poppr")
library("ape") # To visualize the tree using the "nj" function
library("magrittr")
library("adegenet")
library("reshape")
library('vcfR')
library('parallel')
library('hierfstat')

sessionInfo() #hierfstat_0.5-7 = cran version #hierfstat_0.5-9 = github version

##set seed for reproducibility
set.seed(666)

##read in data
vcf<- read.vcfR("Pop_for_pan_260.All.SNP.combined_selected.vcf.gz") #W/o TE's and w 260 strains
Afum_grp<-read.delim("grp_temp.txt", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE)
#name_map<-read.delim("clade_map_K3_20Jan2021.txt", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE, check.names = TRUE)
#Afum_grp<- data.frame(cbind(strain = name_map$name_pop_genome, pop = name_map$clade))


##Convert to genID object
genind_ob<- vcfR2genind(vcf, pop = Afum_grp$pop, ploidy = 1)

#set pop as strata
strata(genind_ob) #NULL
strata(genind_ob)<-Afum_grp
strata(genind_ob)

##Create hierfstat object
my_genind2 <- genind2hierfstat(genind_ob) 

##transform into dosage format changing all 1s to 2s for all cols except the population col
my_genind2[,-1][my_genind2[,-1]==1] <- 2

## check for clones, and if so run clone correction
#Get number of MLG in GO set
mlg(genind_ob)
# Number of Individuals:  260 
# Number of MLG:  260 
sum(strata(genind_ob) ==1)
sum(strata(genind_ob) ==2)
sum(strata(genind_ob) ==3)

#run clone correction
genind_ob_cc<- clonecorrect(genind_ob, strata = ~pop, combine = FALSE)
mlg(genind_ob_cc)

sum(strata(genind_ob_cc) ==1)
sum(strata(genind_ob_cc) ==2)
sum(strata(genind_ob_cc) ==3)
# Number of Individuals:  260 
# Number of MLG:  260 
#= No clones in the dataset

#double using the explicit output generated from calling clonecorrect on a genclone object
genind_clone_ob<- as.genclone(genind_ob, mlgclass = TRUE)
mlg(genind_clone_ob)
genind_clone_ob_cc<- clonecorrect(genind_clone_ob, strata = ~pop, combine = FALSE)
genind_clone_ob_cc
#nope - there are no clones in this set. 

##get Fst values for ea pop and overall
fst.dosage(my_genind2[,-1],pop=my_genind2$pop)
#2         1         3       All 
#0.8045510 0.7775002 0.8718761 0.8179758 

##get confidence interval for ea. population 

#Fst_ea<- betas(my_genind_dip_w_pop,nboot=101,lim=c(0.025,0.975),betaijT=FALSE)
#Fst_ea


##diplodize the haploid data by changing 0 to 11 and 2 to 22
#make copy
#predip<- my_genind2
#change 0 to 11
#predip[,-1][predip[,-1]==0] <- 11
#change 2 to 22
#predip[,-1][predip[,-1]==2] <- 22


#using dosage format to calculate betas
#dosage_betas<- beta.dosage(my_genind2[,-1],inb=FALSE,Mb=FALSE)
#dosage_betas

##Calculate between population Fst values and confidence intervals 

#pairwise Fst CIs - note - you don't want to do this if your population sizes are unequal. You want to use dosage format.
#ci_pairwiseFst_boot1000_dipT<- boot.ppfst(predip, nboot=1000, quant=c(0.025,0.975), diploid=TRUE)

#ci_pairwiseFst$ll

#ci_pairwiseFst$ul

#to calculate pairwise fst note - you don't want to do this if your population sizes are unequal. You want to use dosage format.
#pairWC<- genet.dist(predip, diploid=F, method = "WC84") #pairwise.WCfst



#in dosage format
pairWC_dos <- fs.dosage(my_genind2[,-1], pop = my_genind2[,1])
pairWC_dos$Fst2x2
#2         1         3
#2        NA 0.5768640 0.8709889
#1 0.5768640        NA 0.8599834
#3 0.8709889 0.8599834        NA



##AMOVA to calculate pairwise sequence divergence
#define strata
strata(genind_ob)<- as.data.frame(genind_ob$pop)

#run AMOVA
amova_of_genid2<- poppr.amova(genind_ob,
  hier = ~genind_ob.pop,
  clonecorrect = FALSE,
  within = FALSE,
  filter = FALSE,
  threshold = 0,
  #algorithm = "farthest_neighbor",
  threads = 2,
  missing = "loci",
  method = "ade4",
  #nperm = 0
)

amova_of_genid2

#test for significance in population differentiation
genid2signif<- ade4::randtest(amova_of_genid2, nrepet = 999, alter ="two-sided")
#difference is significant



###calculate total nucleotide diversity 

#in poppr
montab <- mlg.table(genind_ob, strata = ~genind_ob.pop)
monstat <- diversity_stats(montab)
ci_monstat_T<- diversity_ci(genind_ob, n = 1000, rarefy = T, 
                            raw = FALSE, ci = 95, plot = TRUE,
                            center = TRUE, n.rare = 2)
ci_monstat_T

#produces no confidence intervals
ci_monstat_F<- diversity_ci(genind_ob, n = 10000, rarefy = FALSE, 
                            raw = TRUE, ci = 95, n.boot = 15, plot = TRUE, center = TRUE)
ci_monstat_F
