#this script creates STRUCTURE plots from STRUCTURE meanQ outputs, looks for introgression, and does stats and makes a figure to look at the marginal likelihood values for K

#load packages
library(pophelper)
library(grid)
library(gridExtra)
library(svglite)
library(tidyr)
library(multcomp)

setwd("")

#first get the group assignments (note strain_names_in_order.csv is just the second col of the faststructure output Afum_260.nosex)
inds <- read.table("strain_names_in_order.csv",header=FALSE,stringsAsFactors=F)
inds$V1

#load the first itteration of the meanQ files for K2:5
ffiles <- list.files(path=".",pattern="*.meanQ",full.names=T)
flist <- readQ(files=ffiles)

#length(flist[[1]])
#rownames(flist[[1]]) <- (inds$V1)
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
p1 <- plotQ(alignK(sortQ(flist)[1:4]),imgoutput="join",returnplot=T,exportplot=F,quiet=T,basesize=11,
            clustercol=c("#56326E",
                         "#ABA778",
                         "#ED7F6F",
                         "#F2CC85", 
                         "#D4494E",
                         "#7D7A70"), 
                         #"#7E7A70"),
            showlegend=T, legendlab=c("Clade 1","Clade 2","Clade 3","Clade 4","Clade 5"),
            legendkeysize=10,legendtextsize=10,legendmargin=c(2,2,2,0),
            showindlab=F,useindlab=T,showyaxis=T, indlabsize = 4,
            sortind="all", sharedindlab = F, panelspacer=0.4,
            barbordersize=1)
grid.arrange(p1$plot[[1]])

#my_cols_me<-c("#56326E","#ED7F6F","#ABA778","#F2CC85", "#D4494E","#7D7A70")
#scales::show_col(my_cols_me); my_cols_me

#get K=3
flist2<- as.data.frame(flist[2])

admixed<- flist2[apply(flist2[1:3], 1, function(x) all(x < .85)), ] 
admixed
nrow(admixed)
#here there are 8 of them 

#Afum_260.3.Cluster1 Afum_260.3.Cluster2 Afum_260.3.Cluster3
#10-01-02-27            0.152166            0.000000            0.847833
#12-7505220             0.298457            0.000000            0.701542
#AF293                  0.812018            0.152990            0.034992
#CF098                  0.812468            0.153031            0.034501
#CM7632                 0.805103            0.157081            0.037816
#F18304                 0.809874            0.153508            0.036618
#F7763                  0.822984            0.160289            0.016727
#RSF2S8                 0.416723            0.000000            0.583277


#plot admixed individuals only
admixed_names<- rownames(admixed)
admixed_names
#subset input 
only_admixed<- sapply(flist,function(x) x[match(c(admixed_names),rownames(x)),])

#plot
p1 <- plotQ(alignK(sortQ(only_admixed)[2]),returnplot=T,exportplot=F,quiet=T,basesize=20,barsize = .5,
            clustercol=c("#56326E",
                         "#ABA778",
                         "#ED7F6F",
                         "#F2CC85", 
                         "#D4494E",
                         "#7D7A70", 
                         "#7E7A70"),
            showlegend=T, legendlab=c("Clade 1","Clade 2","Clade 3"),
            legendkeysize=10,legendtextsize=10,legendmargin=c(2,2,2,0),
            showindlab=T,useindlab=T,showyaxis=T, indlabsize = 10,
            sortind="all", sharedindlab = F, panelspacer=0.4,
            barbordersize=0)
grid.arrange(p1$plot[[1]])


#make plots from chooseK.py
choose_k <- read.csv("choose_k_results.csv",header=FALSE, row.names = 1, stringsAsFactors=F, check.names = F)

choose_k_sub<- choose_k[,4:ncol(choose_k)]
#transpose
choose_k_t<- data.frame(t(choose_k_sub))

#make long
choose_k_long <- choose_k_t %>% pivot_longer(cols = 2:31, names_to = "Run_n", values_to = "ML")

choose_k_long$Run_n<- gsub("X", "", choose_k_long$Run_n)
choose_k_long$itteration<- gsub("ML K =", "", choose_k_long$itteration)

#plot
p<- ggplot(data = choose_k_long, mapping = aes(x = as.numeric(itteration), y = as.numeric(ML), group = Run_n, color = Run_n)) + 
  geom_line(color="#D4494E", size=1, alpha=.2, linetype=1) +
  theme_minimal() +
  xlab("K") + ylab("Marginal Likelihood")+
  scale_x_discrete(limits=c("1","2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"))
p
#ggsave(plot = p, width = 6, height = 2, dpi = 300, filename = "Marginal_likelyhood.pdf")


#get max cluster assignments for K=3
#rename
colnames(flist2)<- c("1", "2", "3")
flist2$clade_assignment<- colnames(flist2)[max.col(flist2,ties.method="first")]

#get n in each clade 
table(flist2$clade_assignment) #note this is the auto naming for groups "1","2",and "3" where "1" is Clade1, but "2" is actually Clade3 and "3" is actually Clade2 - the colors in the above plots reflect this
#1   2   3 
#200  15  45 

#the total numbers of isolates in each group are consistent with DAPCA 

#format for printing to match the DAPCA sup table. 
colnames(flist2)<- c("prob_clade1", "prob_clade1", "prob_clade1", "clade_assignment")
#reorder to last col is first
print_df<- flist2[,c(ncol(flist2),1:(ncol(flist2)-1))]
print_df$strain<- row.names(print_df)
print_df<- print_df[,c(ncol(print_df),1:(ncol(print_df)-1))]

#write posteriors with assignments for supplemental table
write.table(print_df, "posteriors_fastSTRUCTURE.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


###Run stats on what point is no longer a significant increase in K 
#create differences df - a dataframe of differences between that col n (k) and col n - 1 (k - 1)
choose_k_sub_input<- choose_k_sub[2:nrow(choose_k_sub),]
str(choose_k_sub_input)

#make numeric
choose_k_sub_input<-as.data.frame(lapply(choose_k_sub_input,as.numeric))
class(choose_k_sub_input$V5)

K_diffs<-  data.frame(t(apply(choose_k_sub_input,1,diff)))

#rename confusing col headers
colnames(K_diffs)<- seq(1:14)

#make long
choose_k_long_diffs <- K_diffs %>% pivot_longer(cols = 1:14, names_to = "K", values_to = "ML")

#run anova
#aov(K_diffs_clean)

#change to factor
choose_k_long_diffs$K = as.factor(choose_k_long_diffs$K)
#run linear model and anova
k_diffs_lm <- lm(ML ~ K, data = choose_k_long_diffs)
k_diffs_anova<- anova(k_diffs_lm)
summary(k_diffs_anova)

#posthoc test
# running glht()
post.hoc <- glht(k_diffs_lm, linfct = mcp(K = 'Tukey'))

# displaying the result table with summary()
summary(post.hoc, test = adjusted("bonferroni"))
#does not increase significantly after k = 4

#difference between K vals
K_diffs_clean<- K_diffs[,-1]
abs(mean(K_diffs_clean[,1])) #mean difference between K1 and K2
abs(mean(K_diffs_clean[,2])) #mean difference between K2 and K3
abs(mean(K_diffs_clean[,3])) #mean difference between K3 and K4
abs(mean(K_diffs_clean[,4])) #mean difference between K4 and K5

