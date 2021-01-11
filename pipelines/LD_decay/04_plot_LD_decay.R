#this script creates a plot for LD decay (r2) estimated using PLINK from filtered SNP data. It was written for R version 4.0.1

library(dplyr)
library(stringr)
library(ggplot2)
library(purrr)
library(data.table)

R.Version()

#list files to read in
dfr_1 <-list.files(path = ".",
             pattern = "Clade1_samp*", 
             full.names = T) %>% 
  map_df(~fread(.))
colnames(dfr_1) <- c("dist","rsq")

dfr_2 <-list.files(path = ".",
                   pattern = "Clade2_samp*", 
                   full.names = T) %>% 
  map_df(~fread(.))
colnames(dfr_2) <- c("dist","rsq")

dfr_3 <-list.files(path = ".",
                   pattern = "Clade3_samp*", 
                   full.names = T) %>% 
  map_df(~fread(.))
colnames(dfr_3) <- c("dist","rsq")

#get median value for each distance
rsq_means_Clade1<- aggregate(rsq~dist, data=dfr_1, FUN=function(x) c(mean=mean(x), count=length(x)))
rsq_means_Clade2<- aggregate(rsq~dist, data=dfr_2, FUN=function(x) c(mean=mean(x), count=length(x)))
rsq_means_Clade3<- aggregate(rsq~dist, data=dfr_3, FUN=function(x) c(mean=mean(x), count=length(x)))


#convert and back to get it into the correct format (ugh)
rsq_means_Clade1<- as.matrix(rsq_means_Clade1)
rsq_means_Clade1<- as.data.frame(rsq_means_Clade1)
rsq_means_Clade2<- as.matrix(rsq_means_Clade2)
rsq_means_Clade2<- as.data.frame(rsq_means_Clade2)
rsq_means_Clade3<- as.matrix(rsq_means_Clade3)
rsq_means_Clade3<- as.data.frame(rsq_means_Clade3)

#add group annotation 
rsq_means_Clade1$Clade<- "Clade1"
rsq_means_Clade2$Clade<- "Clade2"
rsq_means_Clade3$Clade<- "Clade3"

#bind all
rsq_means<- rbind(rsq_means_Clade1, rsq_means_Clade2, rsq_means_Clade3)

#print here from supercomputer 
#write.table(rsq_means, "rsq_means.tab", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


#get group median
gm_1<- median(rsq_means_Clade1$rsq.mean)
gm_2<- median(rsq_means_Clade2$rsq.mean)
gm_3<- median(rsq_means_Clade3$rsq.mean)
gm<- data.frame(Clade = c("Clade1", "Clade2", "Clade3"), gm =c(gm_1,gm_2,gm_3))

#find x values closest to these Y values 
C1_y_index<- which.min(abs(rsq_means_Clade1$rsq.mean - gm_1))
gm_1y<- rsq_means_Clade1[C1_y_index, 1]
C2_y_index<- which.min(abs(rsq_means_Clade2$rsq.mean - gm_2))
gm_2y<- rsq_means_Clade1[C2_y_index, 1]
C3_y_index<- which.min(abs(rsq_means_Clade3$rsq.mean - gm_3))
gm_3y<- rsq_means_Clade3[C3_y_index, 1]
gm$gm_y<- c(gm_1y, gm_2y, gm_3y)


#set pallet
myCol <- c(Clade1 = "#56326E",
           Clade2 = "#ED7F6F",
           Clade3 = "#ABA778")

#plot
sp<- ggplot(rsq_means, aes(x=dist,y=rsq.mean))+
  geom_point(shape=20,size=0.15,alpha=0.7, aes(color = Clade)) +
  geom_smooth(aes(color = Clade), show.legend = FALSE) +
  geom_smooth(aes(color = Clade), se=FALSE) +
  labs(x="Pairwise Distance (bp)",y=expression(LD~(r^{2})))+
  theme_bw(base_size=14)+ 
  theme(panel.border=element_blank(),
        axis.ticks=element_blank())
sp + scale_color_manual(values=myCol) + 
  geom_segment(aes(x = gm[1,3], y = .04, xend = gm[1,3], yend = 0), colour = "#56326E",
                arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
  geom_segment(aes(x = gm[2,3], y = .04, xend = gm[2,3], yend = 0), colour = "#ED7F6F",
             arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
  geom_segment(aes(x = gm[3,3], y = .04, xend = gm[3,3], yend = 0), colour = "#ABA778",
               arrow = arrow(length = unit(0.18, "cm"), type = "closed"))

ggsave("LD_decay_plot.pdf", plot = last_plot())





