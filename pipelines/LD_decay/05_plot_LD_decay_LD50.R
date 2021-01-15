library(dplyr)
library(stringr)
library(ggplot2)
library(purrr)
library(data.table)

setwd("~/bigdata/pop_genomics/LD_decay/LD_out")

#get median value for each distance
rsq_means_Clade1 <- read.delim("rsq_means_Clade1.tab",sep="",header=T,check.names=F,stringsAsFactors=F)
rsq_means_Clade2 <- read.delim("rsq_means_Clade2.tab",sep="",header=T,check.names=F,stringsAsFactors=F)
rsq_means_Clade3 <- read.delim("rsq_means_Clade3.tab",sep="",header=T,check.names=F,stringsAsFactors=F)



#sliding window size = n 
n <- 50
#get average over sliding windows 
rsq_means_Clade1_windows<- aggregate(rsq_means_Clade1,list(rep(1:(nrow(rsq_means_Clade1)%/%n+1),each=n,len=nrow(rsq_means_Clade1))),mean)[-1];
rsq_means_Clade2_windows<- aggregate(rsq_means_Clade2,list(rep(1:(nrow(rsq_means_Clade2)%/%n+1),each=n,len=nrow(rsq_means_Clade2))),mean)[-1];
rsq_means_Clade3_windows<- aggregate(rsq_means_Clade3,list(rep(1:(nrow(rsq_means_Clade3)%/%n+1),each=n,len=nrow(rsq_means_Clade3))),mean)[-1];

#add group annotation 
rsq_means_Clade1_windows$Clade<- "Clade1"
rsq_means_Clade2_windows$Clade<- "Clade2"
rsq_means_Clade3_windows$Clade<- "Clade3"

colnames(rsq_means_Clade1_windows) <- c("dist","rsq.mean", "rsq.count", "Clade")
colnames(rsq_means_Clade2_windows) <- c("dist","rsq.mean", "rsq.count", "Clade")
colnames(rsq_means_Clade3_windows) <- c("dist","rsq.mean", "rsq.count", "Clade")
#bind
rsq_means<- data.frame(rbind(rsq_means_Clade1_windows, rsq_means_Clade2_windows, rsq_means_Clade3_windows))


#set pallet
myCol <- c(Clade1 = "#56326E",
           Clade2 = "#ED7F6F",
           Clade3 = "#ABA778")


#get group mean
gm_1<- mean(rsq_means_Clade1$rsq.mean)
gm_2<- mean(rsq_means_Clade2$rsq.mean)
gm_3<- mean(rsq_means_Clade3$rsq.mean)
gm<- data.frame(Clade = c("Clade1", "Clade2", "Clade3"), gm =c(gm_1,gm_2,gm_3))

#find x values closest to these Y values 
C1_y_index<- which.min(abs(rsq_means_Clade1$rsq.mean - gm_1))
gm_1y<- rsq_means_Clade1[C1_y_index, 1]

C2_y_index<- which.min(abs(rsq_means_Clade2$rsq.mean - gm_2))
gm_2y<- rsq_means_Clade2[C2_y_index, 1]

C3_y_index<- which.min(abs(rsq_means_Clade3$rsq.mean - gm_3))
gm_3y<- rsq_means_Clade3[C3_y_index, 1]

#add these new intersection values
gm$gm_y<- c(gm_1y, gm_2y, gm_3y)


#plot
sp<- ggplot(rsq_means, aes(x=dist,y=rsq.mean))+
  geom_line(aes(color = Clade), size=.1, alpha=0.9) +  
  #geom_smooth(aes(color = Clade), show.legend = FALSE, se=FALSE) +
  #geom_smooth(aes(color = Clade), se=FALSE) +
  labs(x="Pairwise Distance (bp)",y=expression(LD~(r^{2})))+
  theme_bw(base_size=14)+ 
  theme(panel.border=element_blank(),
        axis.ticks=element_blank())
sp + scale_color_manual(values=myCol) + 
  #add lines to indicate the point at which each Clade is half decayed
  geom_segment(aes(x = gm[1,3], y = .04, xend = gm[1,3], yend = 0), colour = "#56326E",
               arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
  geom_segment(aes(x = gm[2,3], y = .04, xend = gm[2,3], yend = 0), colour = "#ED7F6F",
               arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
  geom_segment(aes(x = gm[3,3], y = .04, xend = gm[3,3], yend = 0), colour = "#ABA778",
               arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
  scale_x_log10()

setwd("~/bigdata/pop_genomics/LD_decay/plots")
ggsave("LD_decay_wLD50.pdf", plot = last_plot())
