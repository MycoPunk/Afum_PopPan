#this script identifies the number of clusters present in the A fum dataset. 
#using DAPC in the the adegnet package

#set wd
#setwd("~")

#load libraries
library("poppr")
library("ape") # To visualize the tree using the "nj" function
library("magrittr")
library("adegenet")
library("reshape")
library('vcfR')
library('parallel')

#set seed for reproducibility
set.seed(666)

#read in vcf
vcf_me<- read.vcfR("Pop_for_pan_260.All.SNP.combined_selected.vcf.gz") #W/o TE's and w 260 strains

#convert vcf to genlight object
gl_Afum <- vcfR2genlight(vcf_me)
rm(vcf_me)

#run find clusters -note this takes about 90 minutes
grp1<-find.clusters(gl_Afum, max.n.clust=15) #note n-clusters needs to be less than the n of individuals
#Choose the number PCs to retain (>=1): 
#  200 #for the selection of K, there is no benefit to reducing the number of PCs
#Choose the number of clusters (must be >=2): 
#  3 #3 is where the 'elbow' in the curve is

###run DAPC

#first, determine the number of PC's to retain in the analysis 
#run dapc with default settings
dapc1 <- dapc(gl_Afum, grp1$grp, n.da=100, n.pca=15)
temp <- optim.a.score(dapc1)
temp #a-score recommends 3 PCs, with 15 possible at K = 3

#re-run dapc with pca identified above
dapc2 <- dapc(gl_Afum, grp1$grp, n.da=100, n.pca=3)
#dapc2 <- dapc(gl_Afum, grp1$grp, n.da=100, n.pca=6)

#plot the k=3 population groups 
#colors for this project: "#382147","#7D7A70","#ABA778","#F2CC35", "#F2CC85", "#FFCA98", "#F7A583","#ED7F6F","#D4494E"

myCol <- c(clade3 = "#ABA778",
           clade2 ="#ED7F6F",
           clade1 = "#56326E")

#scales::show_col(myCol); myCol

scatter(dapc2, posi.da="bottomright", bg="white",
        pch=20, 
        #cell=0, 
        cstar=0, 
        col=myCol, 
        scree.pca=TRUE,
        posi.pca="bottomleft",
        solid=.4,
        cex = 1.5,
        clab=0,
        leg=TRUE,
        txt.leg=paste("Clade",3:1))


#in a single dimension 
scatter(dapc2,1,1, col=myCol, bg="white",
        scree.da=FALSE, legend=TRUE, solid=.4)

##extract group membership
summary(dapc2)
#how many strains / groups?
dim(dapc2$posterior)
##  260   3
#assignplot(dapc2, subset=1:150,cex.lab=.10,pch=1)
#round(head(dapc2$posterior),3)

#get groups
groups<- data.frame(clade = dapc2$grp)

#look for the most admixed individuals
admixed <- which(apply(dapc2$posterior,1, function(e) all(e<.80))) #those having no more than 0.80 probability of membership to any group
length(admixed) #there are no individuals likely to be admixed between the three populations

##print clade assignments
write.table(groups, "DAPC_groupings_K3_20Jan2021.txt", sep="\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
name_map<-read.delim("namemappingfile_12Oct2020.tab", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE)

#clades<-read.delim("DAPCA_groupings.txt", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE)
clades<- data.frame(cbind(name_pop_genome= rownames(groups), groups))
#match pan_genome_names w/ name_Pan_genome

#change oddly formatted strains
clades$name_pop_genome[clades$name_pop_genome=="AFIS1435CDC_6"] <- "AFIS1435_CDC_6"
clades$name_pop_genome[clades$name_pop_genome=="Afu_343-P/11"] <- "Afu_343-P-11"

#match names
clades$name_Pan_genome <- name_map$name_Pan_genome[match(clades$name_pop_genome, name_map$name_pop_genome)]

#looks good - print it
write.table(clades, "clade_map_K3_20Jan2021.tab", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#write posteriors
write.table(dapc2$posterior, "DAPC_posteriors_K3_20Jan2021.tab", sep="\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
