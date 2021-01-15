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

#shrink
rsq_means_Clade1_sm<- rsq_means_Clade1[1:25000,]
rsq_means_Clade2_sm<- rsq_means_Clade2[1:25000,]
rsq_means_Clade3_sm<- rsq_means_Clade3[1:25000,]

#sliding window size = n 
n <- 10
#get average over sliding windows 
rsq_means_Clade1_windows<- aggregate(rsq_means_Clade1_sm,list(rep(1:(nrow(rsq_means_Clade1_sm)%/%n+1),each=n,len=nrow(rsq_means_Clade1_sm))),mean)[-1];
rsq_means_Clade2_windows<- aggregate(rsq_means_Clade2_sm,list(rep(1:(nrow(rsq_means_Clade2_sm)%/%n+1),each=n,len=nrow(rsq_means_Clade2_sm))),mean)[-1];
rsq_means_Clade3_windows<- aggregate(rsq_means_Clade3_sm,list(rep(1:(nrow(rsq_means_Clade3_sm)%/%n+1),each=n,len=nrow(rsq_means_Clade3_sm))),mean)[-1];

#add group annotation 
rsq_means_Clade1_windows$Clade<- "Clade1"
rsq_means_Clade2_windows$Clade<- "Clade2"
rsq_means_Clade3_windows$Clade<- "Clade3"

colnames(rsq_means_Clade1_windows) <- c("dist","rsq.mean", "rsq.count", "Clade")
colnames(rsq_means_Clade2_windows) <- c("dist","rsq.mean", "rsq.count", "Clade")
colnames(rsq_means_Clade3_windows) <- c("dist","rsq.mean", "rsq.count", "Clade")
#bind
rsq_means_sm<- data.frame(rbind(rsq_means_Clade1_windows, rsq_means_Clade2_windows, rsq_means_Clade3_windows))


#set pallet
myCol <- c(Clade1 = "#56326E",
           Clade2 = "#ED7F6F",
           Clade3 = "#ABA778")
#plot
sp<- ggplot(rsq_means_sm, aes(x=dist,y=rsq.mean))+
  #geom_point(shape=20,size=0.15,alpha=0.7, aes(color = Clade)) +
  geom_line(aes(color = Clade), size=.1, alpha=0.9) +
  geom_smooth(aes(color = Clade), show.legend = FALSE) +
  geom_smooth(aes(color = Clade), se=FALSE) +
  labs(x="Pairwise Distance (bp)",y=expression(LD~(r^{2})))+
  theme_bw(base_size=14)+ 
  theme(panel.border=element_blank(),
        axis.ticks=element_blank())
sp + scale_color_manual(values=myCol)

setwd("~/bigdata/pop_genomics/LD_decay/plots")
ggsave("LD_decay_zoomed_in.pdf", plot = last_plot())
