#This script creates a supplemental figure to look at the influence of sample size on LD decay estimates.

library(dplyr)
library(stringr)
library(ggplot2)
library(purrr)
library(data.table)
library(ggpubr)

setwd("~/bigdata/pop_genomics/LD_decay_downsample_test/LD_out")

#set seed for reproducibility
set.seed(666)

#get median value for each distance
rsq_means_n_5 <- read.delim("~/bigdata/pop_genomics/LD_decay_downsample_test/LD_out/rsq_means_n5.tab",sep="",header=T,check.names=F,stringsAsFactors=F)
rsq_means_n_10 <- read.delim("~/bigdata/pop_genomics/LD_decay_downsample_test/LD_out/rsq_means_n10.tab",sep="",header=T,check.names=F,stringsAsFactors=F)
rsq_means_n_20 <- read.delim("~/bigdata/pop_genomics/LD_decay_downsample_test/LD_out/rsq_means_n20.tab",sep="",header=T,check.names=F,stringsAsFactors=F)
rsq_means_n_50<- read.delim("~/bigdata/pop_genomics/LD_decay_downsample_test/LD_out/rsq_means_n50.tab",sep="",header=T,check.names=F,stringsAsFactors=F)
rsq_means_n_100<- read.delim("~/bigdata/pop_genomics/LD_decay_downsample_test/LD_out/rsq_means_n100.tab",sep="",header=T,check.names=F,stringsAsFactors=F)


#add group annotation 
rsq_means_n_5$group<- "n_5" 
rsq_means_n_10$group<- "n_10" 
rsq_means_n_20$group<- "n_20"
rsq_means_n_50$group<- "n_50"
rsq_means_n_100$group<- "n_100"


#get intersection of R2 half way between R2max and R2 mean 
#get y w/o averaging
LD50_n_5_y<- min(rsq_means_n_5$rsq.mean) + (max(rsq_means_n_5$rsq.mean) - min(rsq_means_n_5$rsq.mean)) /2
LD50_n_10_y<- min(rsq_means_n_10$rsq.mean) + (max(rsq_means_n_10$rsq.mean) - min(rsq_means_n_10$rsq.mean)) /2
LD50_n_20_y<- min(rsq_means_n_20$rsq.mean) + (max(rsq_means_n_20$rsq.mean) - min(rsq_means_n_20$rsq.mean)) /2
LD50_n_50_y<- min(rsq_means_n_50$rsq.mean) + (max(rsq_means_n_50$rsq.mean) - min(rsq_means_n_50$rsq.mean)) /2
LD50_n_100_y<- min(rsq_means_n_100$rsq.mean) + (max(rsq_means_n_100$rsq.mean) - min(rsq_means_n_100$rsq.mean)) /2

#calculate x vals 
LD50_n_5_x <- approx(x = rsq_means_n_5$rsq.mean, y = rsq_means_n_5$dist, xout = LD50_n_5_y)$y
LD50_n_5_x

LD50_n_10_x <- approx(x = rsq_means_n_10$rsq.mean, y = rsq_means_n_10$dist, xout = LD50_n_10_y)$y
LD50_n_10_x

LD50_n_20_x <- approx(x = rsq_means_n_20$rsq.mean, y = rsq_means_n_20$dist, xout = LD50_n_20_y)$y
LD50_n_20_x

LD50_n_50_x <- approx(x = rsq_means_n_50$rsq.mean, y = rsq_means_n_50$dist, xout = LD50_n_50_y)$y
LD50_n_50_x

LD50_n_100_x <- approx(x = rsq_means_n_100$rsq.mean, y = rsq_means_n_100$dist, xout = LD50_n_100_y)$y
LD50_n_100_x



#concat LD50y, LD50x, and group designations
gm<- data.frame(group = c("n_5", "n_10", "n_20", "n_50", "n_100"), 
                LD50y =c(LD50_n_5_y, 
                         LD50_n_10_y,
                         LD50_n_20_y,
                         LD50_n_50_y,
                         LD50_n_100_y),
                LD50x =c(LD50_n_5_x, 
                         LD50_n_10_x,
                         LD50_n_20_x,
                         LD50_n_50_x,
                         LD50_n_100_x))

#set pallet
myCol <- c(n_5 = "dark grey",
           n_10 = "#F7A583",
           n_20 = "#D4494E",
           n_50 = "#F2CC35",
           n_100 = "#5A51F5")

#bind
rsq_means_grand<- data.frame(rbind(rsq_means_n_5,
                                   rsq_means_n_10,
                                   rsq_means_n_20,
                                   rsq_means_n_50,
                                   rsq_means_n_100))


nrow(rsq_means_grand) #too many data points for easy plot rendering
#subset for rendering
random_input<- rsq_means_grand[sample(nrow(rsq_means_grand), 10000),]
nrow(random_input)

##plot all
sp<- ggplot(data=random_input,aes(x=dist,y=rsq.mean))+
  geom_line(data=random_input, aes(color=group),size=0.15,alpha=0.9)+
  #geom_smooth(aes(color = group), show.legend = FALSE, size=.8) +
  #geom_smooth(aes(color = group), se=FALSE, size=.8) +
  labs(x="Distance (BP)",y=expression(LD~(r^{2})))+
  theme_bw() + scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8))
all_print<- sp + scale_color_manual(values=myCol)
#ggsave("LD_decay_wLD50_w.o.subsettting.pdf", plot = last_plot(), width = 8, height = 4)


##plot zoomed in
##plot zoomed in
#subset to only the first 20,000 bp - but 
rsq_means_n_5_zoom<- rsq_means_n_5[1:5000,]
rsq_means_n_10_zoom<- rsq_means_n_10[1:5000,]
rsq_means_n_20_zoom<- rsq_means_n_20[1:5000,]
rsq_means_n_50_zoom<- rsq_means_n_50[1:5000,]
rsq_means_n_100_zoom<- rsq_means_n_100[1:5000,]

#retain the first positions for detailed graphing 
rsq_means_n_5_zoom_first<- rsq_means_n_5_zoom[1:20,]
rsq_means_n_10_zoom_first<- rsq_means_n_10_zoom[1:20,]
rsq_means_n_20_zoom_first<- rsq_means_n_20_zoom[1:20,]
rsq_means_n_50_zoom_first<- rsq_means_n_50_zoom[1:20,]
rsq_means_n_100_zoom_first<- rsq_means_n_100_zoom[1:20,]

#retain the last positions for subsetting (so the graph is managibly small)
rsq_means_n_5_zoom_last<- rsq_means_n_5_zoom[21:nrow(rsq_means_n_5_zoom),]
rsq_means_n_10_zoom_last<- rsq_means_n_10_zoom[21:nrow(rsq_means_n_10_zoom),]
rsq_means_n_20_zoom_last<- rsq_means_n_20_zoom[21:nrow(rsq_means_n_20_zoom),]
rsq_means_n_50_zoom_last<- rsq_means_n_50_zoom[21:nrow(rsq_means_n_50_zoom),]
rsq_means_n_100_zoom_last<- rsq_means_n_100_zoom[21:nrow(rsq_means_n_100_zoom),]

#subset the lasts
rsq_means_n_5_zoom_last_subet<- rsq_means_n_5_zoom_last[sample(nrow(rsq_means_n_5_zoom_last), 4000),]
rsq_means_n_10_zoom_last_subet<- rsq_means_n_10_zoom_last[sample(nrow(rsq_means_n_10_zoom_last), 4000),]
rsq_means_n_20_zoom_last_subet<- rsq_means_n_20_zoom_last[sample(nrow(rsq_means_n_20_zoom_last), 4000),]
rsq_means_n_50_zoom_last_subet<- rsq_means_n_50_zoom_last[sample(nrow(rsq_means_n_50_zoom_last), 4000),]
rsq_means_n_100_zoom_last_subet<- rsq_means_n_100_zoom_last[sample(nrow(rsq_means_n_100_zoom_last), 4000),]


#bind first and last
rsq_means_n_5_zoom_zoom_all<- rbind(rsq_means_n_5_zoom_first, rsq_means_n_5_zoom_last_subet)
rsq_means_n_10_zoom_zoom_all<- rbind(rsq_means_n_10_zoom_first, rsq_means_n_10_zoom_last_subet)
rsq_means_n_20_zoom_zoom_all<- rbind(rsq_means_n_20_zoom_first, rsq_means_n_20_zoom_last_subet)
rsq_means_n_50_zoom_zoom_all<- rbind(rsq_means_n_50_zoom_first, rsq_means_n_50_zoom_last_subet)
rsq_means_n_100_zoom_zoom_all<- rbind(rsq_means_n_100_zoom_first, rsq_means_n_100_zoom_last_subet)

#bind
rsq_means_grand_zoom<- data.frame(rbind(rsq_means_n_5_zoom_zoom_all,
                                        rsq_means_n_10_zoom_zoom_all,
                                        rsq_means_n_20_zoom_zoom_all,
                                        rsq_means_n_50_zoom_zoom_all,
                                        rsq_means_n_100_zoom_zoom_all))

#plot
sp_zoom<- ggplot(data=rsq_means_grand_zoom,aes(x=dist,y=rsq.mean))+
  geom_line(data=rsq_means_grand_zoom, aes(color=group),size=0.15,alpha=0.9)+
  #geom_smooth(aes(color = group), show.legend = FALSE) +
  #geom_smooth(aes(color = group), se=FALSE) +
  labs(x="Distance (log BP)",y=expression(LD~(r^{2})))+
  theme_bw()
zoom_print<- sp_zoom + scale_color_manual(values=myCol) + 
  #add lines to indicate the point at which each group is half decayed
  geom_segment(aes(x = gm[1,3], y = .08, xend = gm[1,3], yend = 0), colour = "dark grey",
               arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
  geom_segment(aes(x = gm[2,3], y = .08, xend = gm[2,3], yend = 0), colour = "#F7A583",
               arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
  geom_segment(aes(x = gm[3,3], y = .08, xend = gm[3,3], yend = 0), colour = "#D4494E",
               arrow = arrow(length = unit(0.15, "cm"), type = "closed")) + 
  geom_segment(aes(x = gm[4,3], y = .08, xend = gm[4,3], yend = 0), colour = "#F2CC35",
               arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
  geom_segment(aes(x = gm[5,3], y = .08, xend = gm[5,3], yend = 0), colour = "#5A51F5",
               arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
   scale_x_log10() +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8))


#plot together
p<-ggarrange(all_print, zoom_print,
             #labels = c("", ""),
             ncol = 1, nrow = 2)
p
ggsave("LD50_downsample_control.pdf",p, width=9, height=7, units="in")

