##this script runs population stats on a vcf file and population defined previously using DAPCA
##Last updated 15.Apr.2021

##set wd
#setwd("")

##load libraries
library("poppr")
library("ape") 
library("magrittr")
library("adegenet")
library("reshape")
library('vcfR')
library('parallel')
library('hierfstat')

#sessionInfo()

##set seed for reproducibility
set.seed(666)

##read in data
vcf<- read.vcfR("262_strains.selected.SNP.NO_TEs.vcf.gz") #W/o TE's and w 262 strains
Afum_grp<-read.delim("grp_temp.txt", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE)

##Convert to genID object
genind_ob<- vcfR2genind(vcf, pop = Afum_grp$pop, ploidy = 1)

##Create hierfstat object
my_genind2 <- genind2hierfstat(genind_ob) 

##transform into dosage format changeing all 1s to 2s for all cols except the population col
my_genind2[,-1][my_genind2[,-1]==1] <- 2

## check for clones, and if so run clone correction
#Get number of MLG in GO set
mlg(genind_ob)
# Number of Individuals:  262 
# Number of MLG:  262 
sum(strata(genind_ob, ~genind_ob.pop) ==1)
sum(strata(genind_ob, ~genind_ob.pop) ==2)
sum(strata(genind_ob, ~genind_ob.pop) ==3)

#run clone correction
genind_ob_cc<- clonecorrect(genind_ob,  strata = ~genind_ob.pop, combine = TRUE)
mlg(genind_ob_cc)
sum(strata(genind_ob_cc, ~genind_ob.pop) ==1)
sum(strata(genind_ob_cc, ~genind_ob.pop) ==2)
sum(strata(genind_ob_cc, ~genind_ob.pop) ==3)
# Number of Individuals:  262 
# Number of MLG:  262 
#= No clones in the dataset


##get Fst values for ea pop and overall
fst.dosage(my_genind2[,-1],pop=my_genind2$pop)
#2         1         3       All 
#0.8045356 0.7774566 0.8721277 0.8180399 

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
#1         2         3
#1 NA 0.8386362 0.8154802
#2 NA        NA 0.5492461
#3 NA        NA        NA

#ci_pairwiseFst$ul
#1         2         3
#1 NA 0.8436003 0.8207094
#2 NA        NA 0.5609904
#3 NA        NA        NA

#to calculate pairwise fst note - you don't want to do this if your population sizes are unequal. You want to use dosage format.
#pairWC<- genet.dist(predip, diploid=F, method = "WC84") #pairwise.WCfst
#2         1
#1 0.5561156          
#3 0.8408365 0.8181832


#in dosage format
pairWC_dos <- fs.dosage(my_genind2[,-1], pop = my_genind2[,1])
pairWC_dos$Fst2x2
#2         1         3
#2        NA 0.5772499 0.8710694
#1 0.5772499        NA 0.8600231
#3 0.8710694 0.8600231        NA



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

#test for significance in population differentiation
genid2signif<- ade4::randtest(amova_of_genid2, nrepet = 999, alter ="two-sided")
genid2signif



###calculate total nucleotide diversity 
#convert to DNAbin format 

#bin_ob<- vcfR2DNAbin(vcf)
#nuc.div(bin_ob)

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
