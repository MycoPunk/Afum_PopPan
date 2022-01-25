##this script runs population stats on a vcf file and populations defined previously using DAPC
#note, unless you're going to subset the data, you probably want to run this on the cluster w/ highmem.
##Last updated 25.Jan.2022

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

sessionInfo() #hierfstat_0.5-7 = cran version #hierfstat_0.5-9 = github version

##set seed for reproducibility
set.seed(666)

##read in data
vcf<- read.vcfR("Pop_for_pan_260.v2.All.SNP.combined_selected.NO.TE.vcf.gz") #W/o TE's and w 260 strains
Afum_grp<-read.delim("clade_map_K3_3Jan2022.txt", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE, check.names = TRUE)

#check that all sample names are retained and in the right order
setdiff(colnames(vcf@gt)[-1], Afum_grp$name_pop_genome_new)# character(0) = TRUE - all the names are correct and all the samples are in there
all(colnames(vcf@gt)[-1] == Afum_grp$name_pop_genome_new) # FALSE - samples are not in the correct order

#put samples in the right order (to match vcf order)
to_be_ordered_strata<- as.data.frame(cbind(Afum_grp$name_pop_genome_new, Afum_grp$clade))
sort_on<- as.data.frame(colnames(vcf@gt)[-1])
ordered_strata<- to_be_ordered_strata[order(match(to_be_ordered_strata[,1],sort_on[,1])),]

#check that it's in order now
all(colnames(vcf@gt)[-1] == ordered_strata$V1) # TRUE - pop now ordered to match vcf. 

##Convert to genID object
genind_ob<- vcfR2genind(vcf, ploidy = 1)

#set pop as strata
colnames(ordered_strata)<- c("strain", "pop")
strata(genind_ob) #NULL
strata(genind_ob)<-ordered_strata
strata(genind_ob) # now set to clade

#check pop specifically 
head(pop(genind_ob)) #nope
setPop(genind_ob) <- ~pop
head(pop(genind_ob)) #now set to the pop strata

#clean up large vcf to free memory
rm(vcf)

#create hierfstat object
my_genind2 <- genind2hierfstat(genind_ob) 

#transform into dosage format changing all 1s to 2s for all cols except the population col
my_genind2[,-1][my_genind2[,-1]==1] <- 2

## check for clones, and if so run clone correction
#Get number of MLG in vcf
mlg(genind_ob)
# Number of Individuals:  260 
# Number of MLG:  260 
sum(strata(genind_ob) ==1)
sum(strata(genind_ob) ==2)
sum(strata(genind_ob) ==3)
#looks like there are no clones - as many MLGs as there are individuals.

#this will run clone correction if needed (it is not here so I'm hashing it for this project)
#genind_ob_cc<- clonecorrect(genind_ob, strata = ~pop, combine = FALSE)
#mlg(genind_ob_cc)

#sum(strata(genind_ob_cc) ==1)
#sum(strata(genind_ob_cc) ==2)
#sum(strata(genind_ob_cc) ==3)


#double check using the explicit output generated from calling clonecorrect on a genclone object
#genind_clone_ob<- as.genclone(genind_ob, mlgclass = TRUE)
#mlg(genind_clone_ob)
#genind_clone_ob_cc<- clonecorrect(genind_clone_ob, strata = ~pop, combine = FALSE)#
#genind_clone_ob_cc 

##get Fst values for ea pop and overall
fst.dosage(my_genind2[,-1],pop=my_genind2$pop)
#        2         1         3       All 
#0.8252919 0.7669710 0.8794938 0.8239189 

##get pairwise Fst values
pairWC_dos <- fs.dosage(my_genind2[,-1], pop = my_genind2[,1])
pairWC_dos$Fst2x2
#          2         1         3
#2        NA 0.5308424 0.8873124
#1 0.5308424        NA 0.8592142
#3 0.8873124 0.8592142        NA


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

#$componentsofcovariance
#                               Sigma         %
#Variations  Between samples 45009.85  69.79054
#Variations  Within samples  19482.91  30.20946
#Total variations            64492.76 100.00000

#$statphi
#                        Phi
#Phi-samples-total 0.6979054


#test for significance in population differentiation
genid2signif<- ade4::randtest(amova_of_genid2, nrepet = 999)
genid2signif
#difference is significant
