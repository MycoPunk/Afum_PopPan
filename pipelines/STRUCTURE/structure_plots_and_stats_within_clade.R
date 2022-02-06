#this script creates STRUCTURE plots from STRUCTURE meanQ outputs, looks for introgression, 
#this script is to look at C1 and C2 respectively at K=5 for each, based on DAPCA results

#load packages
library(pophelper)
library(grid)
library(gridExtra)
library(svglite)
library(tidyr)
library(multcomp)


#####for C1
setwd("~/Desktop/Project_Afum_pangenome_3/MeanQ_C1")
#first get the group assignments (note strain_names_in_order.csv is just the second col of the faststructure output Afum_260.nosex)
#inds <- read.table("strain_names_in_order.csv",header=FALSE,stringsAsFactors=F)
inds <- read.table("../strain_names_in_order_C1.csv",header=FALSE,stringsAsFactors=F)
inds$V1

#load the first itteration of the meanQ files for K3:5
ffiles <- list.files(path=".",pattern="*.meanQ",full.names=T)
flist <- readQ(files=ffiles)

if(length(unique(sapply(flist,nrow)))==1) flist <- lapply(flist,"rownames<-",inds$V1)
# show row names of all runs and all samples
lapply(flist, rownames)

# view head of first converted file
head(flist[[1]])
#view file attributes like this
attributes(flist[[1]])

#tabulateQ
tr1 <- tabulateQ(qlist=flist)
tabulateQ(flist)
tabulateQ(flist, writetable=TRUE, sorttable=TRUE)

#summarize Q
sr1 <- summariseQ(tr1)
summariseQ(tr1, writetable=TRUE)

#plot
p1 <- plotQ(alignK(sortQ(flist)[1:3]),imgoutput="join",returnplot=T,exportplot=F,quiet=T,basesize=11,
            clustercol=c("#73726A",
                         "#A6A498",
                         "#D9D8D7",
                         "#F2F2F2",
                         "#000000"),
            showlegend=T, legendlab=c("sub-clade 1.1","sub-clade 1.2","sub-clade 1.3","sub-clade 1.4","sub-clade 1.5"),
            legendkeysize=10,legendtextsize=10,legendmargin=c(2,2,2,0),
            showindlab=F,useindlab=T,showyaxis=T, indlabsize = 4,
            sortind="all", sharedindlab = F, panelspacer=0.4,
            barbordersize=1)
p2<-grid.arrange(p1$plot[[1]])

ggsave(file="subSTRUCTURE_C1.pdf",device="pdf", p2, width=12, height=8, units="in")

#get K=5
flist2<- as.data.frame(flist[3])

admixed<- flist2[apply(flist2[1:5], 1, function(x) all(x < .85)), ] 
nrow(admixed)
#here there are 87 of them! 
admixed

#plot admixed individuals only
admixed_names<- rownames(admixed)
admixed_names
#subset input 
only_admixed<- sapply(flist,function(x) x[match(c(admixed_names),rownames(x)),])

#plot
p1 <- plotQ(alignK(sortQ(only_admixed)[3]),returnplot=T,exportplot=F,quiet=T,basesize=20,barsize = .5,
            clustercol=c("#73726A",
                         "#A6A498",
                         "#D9D8D7",
                         "#F2F2F2",
                         "#000000"),
            showlegend=T, legendlab=c("sub-clade 1.1","sub-clade 1.2","sub-clade 1.3","sub-clade 1.4","sub-clade 1.5"),
            legendkeysize=10,legendtextsize=10,legendmargin=c(2,2,2,0),
            showindlab=F,useindlab=T,showyaxis=T, indlabsize = 10,
            sortind="all", sharedindlab = F, panelspacer=0.4,
            barbordersize=0)
grid.arrange(p1$plot[[1]])


#get max cluster assignments for K=3
#rename
colnames(flist2)<- c("1", "2", "3", "4", "5")
flist2$clade_assignment<- colnames(flist2)[max.col(flist2,ties.method="first")]

#get n in each clade 
table(flist2$clade_assignment) #note this is the auto naming for groups, and they don't match up perfictly with the DAPC results
#1  2  3  4  5 
#68 54  5 56 17 

#NOTE- the total numbers of isolates in each group are not consistent with DAPCA, and groups are not equivelent 

#format for printing to match the DAPCA sup table. 
colnames(flist2)<- c("prob. sub-clade1.1", "prob. sub-clade1.2", "prob. sub-clade1.3", "prob. sub-clade1.4", "prob. sub-clade1.5", "clade_assignment")
#reorder to last col is first
print_df<- flist2[,c(ncol(flist2),1:(ncol(flist2)-1))]
print_df$strain<- row.names(print_df)
print_df<- print_df[,c(ncol(print_df),1:(ncol(print_df)-1))]

#read in mapping file to rename strains to match paper naming conventions
Afum_grp<-read.delim("../clade_map_K3_3Jan2022.txt", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE, check.names = TRUE)
print_df$strain<- Afum_grp$name_to_use_in_paper[match(print_df$strain, Afum_grp$name_pop_genome_new)]

#write posteriors with assignments for supplemental table
write.table(print_df, "posteriors_fastSTRUCTURE_C1.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


##### for C2
gc()
setwd("~/Desktop/Project_Afum_pangenome_3/MeanQ_C2")
inds <- read.table("../strain_names_in_order_C2.csv",header=FALSE,stringsAsFactors=F)
inds$V1

#load the first itteration of the meanQ files for K2:5
ffiles <- list.files(path=".",pattern="*.meanQ",full.names=T)
flist <- readQ(files=ffiles)

if(length(unique(sapply(flist,nrow)))==1) flist <- lapply(flist,"rownames<-",inds$V1)
# show row names of all runs and all samples
lapply(flist, rownames)

# view head of first converted file
head(flist[[1]])
#view file attributes like this
attributes(flist[[1]])

#tabulateQ
tr1 <- tabulateQ(qlist=flist)
tabulateQ(flist)
tabulateQ(flist, writetable=TRUE, sorttable=TRUE)

#summarize Q
sr1 <- summariseQ(tr1)
summariseQ(tr1, writetable=TRUE)

#plot
p1 <- plotQ(alignK(sortQ(flist)[1:3]),imgoutput="join",returnplot=T,exportplot=F,quiet=T,basesize=11,
            clustercol=c("#73726A",
                         "#A6A498",
                         "#D9D8D7",
                         "#F2F2F2",
                         "#000000"),
            showlegend=T, legendlab=c("sub-clade 1.1","sub-clade 1.2","sub-clade 1.3","sub-clade 1.4","sub-clade 1.5"),
            legendkeysize=10,legendtextsize=10,legendmargin=c(2,2,2,0),
            showindlab=F,useindlab=T,showyaxis=T, indlabsize = 4,
            sortind="all", sharedindlab = F, panelspacer=0.4,
            barbordersize=1)
p2<-grid.arrange(p1$plot[[1]])
ggsave(file="subSTRUCTURE_C2.pdf",device="pdf", p2, width=12, height=8, units="in")

#get K=5
flist2<- as.data.frame(flist[3])

admixed<- flist2[apply(flist2[1:5], 1, function(x) all(x < .85)), ] 
nrow(admixed)
#here there are 7 of them 
admixed
#08-12-12-13
#08-31-08-91
#08-36-03-25  
#10-01-02-27   
#12-7505446
#B5859      
#F13619

#plot admixed individuals only
admixed_names<- rownames(admixed)
admixed_names
#subset input 
only_admixed<- sapply(flist,function(x) x[match(c(admixed_names),rownames(x)),])

#plot
p1 <- plotQ(alignK(sortQ(only_admixed)[3]),returnplot=T,exportplot=F,quiet=T,basesize=20,barsize = .5,
            clustercol=c("#73726A",
                         "#A6A498",
                         "#D9D8D7",
                         "#F2F2F2",
                         "#000000"),
            showlegend=T, legendlab=c("sub-clade 1.1","sub-clade 1.2","sub-clade 1.3","sub-clade 1.4","sub-clade 1.5"),
            legendkeysize=10,legendtextsize=10,legendmargin=c(2,2,2,0),
            showindlab=T,useindlab=T,showyaxis=T, indlabsize = 10,
            sortind="all", sharedindlab = F, panelspacer=0.4,
            barbordersize=0)
grid.arrange(p1$plot[[1]])


#get max cluster assignments for K=3
#rename
colnames(flist2)<- c("1", "2", "3", "4", "5")
flist2$clade_assignment<- colnames(flist2)[max.col(flist2,ties.method="first")]

#get n in each clade 
table(flist2$clade_assignment) #note this is the auto naming for groups
#1  2  3  4  5 
#9  2 28  5  1  

#NOTE: again, the total numbers of isolates in each group are NOT consistent with DAPCA 

#format for printing to match the DAPCA sup table. 
colnames(flist2)<- c("prob. sub-clade2.1", "prob. sub-clade2.2", "prob. sub-clade2.3", "prob. sub-clade2.4", "prob. sub-clade2.5", "clade_assignment")
#reorder to last col is first
print_df<- flist2[,c(ncol(flist2),1:(ncol(flist2)-1))]
print_df$strain<- row.names(print_df)
print_df<- print_df[,c(ncol(print_df),1:(ncol(print_df)-1))]

#change strain names to match paper naming conventions
print_df$strain<- Afum_grp$name_to_use_in_paper[match(print_df$strain, Afum_grp$name_pop_genome_new)]

#write posteriors with assignments for supplemental table
write.table(print_df, "posteriors_fastSTRUCTURE_C2.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
