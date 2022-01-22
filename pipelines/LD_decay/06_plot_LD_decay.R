library(dplyr)
library(stringr)
library(ggplot2)
library(purrr)
library(data.table)
library(ggpubr)

setwd("~/bigdata/pop_genomics/LD_decay_3/LD_out")

#set seed for reproducibility
set.seed(666)

#get median value for each distance
rsq_means_all <- read.delim("~/bigdata/pop_genomics/LD_decay_3/LD_out_nosubset/rsq_means_all_samples_all_clades.tab",sep="",header=T,check.names=F,stringsAsFactors=F)
rsq_means_Clade1_no_sample <- read.delim("~/bigdata/pop_genomics/LD_decay_4/LD_out_save_no_subset/rsq_means_all_samples_Clade_1_all_strains.tab",sep="",header=T,check.names=F,stringsAsFactors=F)
rsq_means_Clade2_no_sample <- read.delim("~/bigdata/pop_genomics/LD_decay_4/LD_out_save_no_subset/rsq_means_all_samples_Clade_2_all_strains.tab",sep="",header=T,check.names=F,stringsAsFactors=F)
rsq_means_Clade3_no_sample <- read.delim("~/bigdata/pop_genomics/LD_decay_4/LD_out_save_no_subset/rsq_means_all_samples_Clade_3_all_strains.tab",sep="",header=T,check.names=F,stringsAsFactors=F)

#add group annotation 
rsq_means_all$Clade<- "All"
rsq_means_Clade1_no_sample$Clade<- "Clade1"
rsq_means_Clade2_no_sample$Clade<- "Clade2"
rsq_means_Clade3_no_sample$Clade<- "Clade3"

#get intersection of R2 half way between R2max and R2 mean 
#get y w/o averaging
LD50_all_y<- min(rsq_means_all$rsq.mean) + (max(rsq_means_all$rsq.mean) - min(rsq_means_all$rsq.mean)) /2
LD50_C1_y<- min(rsq_means_Clade1_no_sample$rsq.mean) + (max(rsq_means_Clade1_no_sample$rsq.mean) - min(rsq_means_Clade1_no_sample$rsq.mean)) /2
LD50_C2_y<- min(rsq_means_Clade2_no_sample$rsq.mean) + (max(rsq_means_Clade2_no_sample$rsq.mean) - min(rsq_means_Clade2_no_sample$rsq.mean)) /2
LD50_C3_y<- min(rsq_means_Clade3_no_sample$rsq.mean) + (max(rsq_means_Clade3_no_sample$rsq.mean) - min(rsq_means_Clade3_no_sample$rsq.mean)) /2

#calculate x vals 
LD50_all_x <- approx(x = rsq_means_all$rsq.mean, y = rsq_means_all$dist, xout = LD50_all_y)$y
LD50_all_x

LD50_C1_x <- approx(x = rsq_means_Clade1_no_sample$rsq.mean, y = rsq_means_Clade1_no_sample$dist, xout = LD50_C1_y)$y
LD50_C1_x

LD50_C2_x <- approx(x = rsq_means_Clade2_no_sample$rsq.mean, y = rsq_means_Clade2_no_sample$dist, xout = LD50_C2_y)$y
LD50_C2_x

LD50_C3_x <- approx(x = rsq_means_Clade3_no_sample$rsq.mean, y = rsq_means_Clade3_no_sample$dist, xout = LD50_C3_y)$y
LD50_C3_x

#concat LD50y, LD50x, and group designations
gm<- data.frame(Clade = c("All", "Clade1", "Clade2", "Clade3"), 
                LD50y =c(LD50_all_y, 
                         LD50_C1_y,
                         LD50_C2_y,
                         LD50_C3_y),
                LD50x =c(LD50_all_x, 
                         LD50_C1_x,
                         LD50_C2_x,
                         LD50_C3_x))


#set pallet
myCol <- c(All = "dark grey",
           Clade1 = "#56326E",
           Clade2 = "#ED7F6F",
           Clade3 = "#ABA778")

#bind
rsq_means_grand<- data.frame(rbind(rsq_means_Clade1_no_sample,
                             rsq_means_Clade2_no_sample,
                             rsq_means_Clade3_no_sample,
                             rsq_means_all))


nrow(rsq_means_grand) #too many data points for easy plot rendering
#subset for east rendering
random_input<- rsq_means_grand[sample(nrow(rsq_means_grand), 20000),]
nrow(random_input)

##plot all
sp<- ggplot(data=random_input,aes(x=dist,y=rsq.mean))+
  geom_line(data=random_input, aes(color=Clade),size=0.1,alpha=0.9)+
  geom_smooth(aes(color = Clade), show.legend = FALSE, size=.8) +
  geom_smooth(aes(color = Clade), se=FALSE, size=.8) +
  labs(x="Distance (BP)",y=expression(LD~(r^{2})))+
  theme_bw() + scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8))
all_print<- sp + scale_color_manual(values=myCol)
#ggsave("LD_decay_wLD50_w.o.subsettting.pdf", plot = last_plot(), width = 8, height = 4)


##plot zoomed in
#subset to only the first 10,000 bp
rsq_means_all_zoom<- rsq_means_all[1:20000,]
rsq_means_C1_zoom<- rsq_means_Clade1_no_sample[1:20000,]
rsq_means_C2_zoom<- rsq_means_Clade2_no_sample[1:20000,]
rsq_means_C3_zoom<- rsq_means_Clade3_no_sample[1:20000,]


rsq_means_grand_zoom<- data.frame(rbind(rsq_means_C1_zoom,
                                        rsq_means_C2_zoom,
                                        rsq_means_C3_zoom,
                                        rsq_means_all_zoom))


#plot
sp_zoom<- ggplot(data=rsq_means_grand_zoom,aes(x=dist,y=rsq.mean))+
  geom_line(data=rsq_means_grand_zoom, aes(color=Clade),size=0.1,alpha=0.9)+
  #geom_smooth(aes(color = Clade), show.legend = FALSE) +
  #geom_smooth(aes(color = Clade), se=FALSE) +
  labs(x="Distance (BP) log scale",y=expression(LD~(r^{2})))+
  theme_bw()
zoom_print<- sp_zoom + scale_color_manual(values=myCol) + 
  #add lines to indicate the point at which each Clade is half decayed
  geom_segment(aes(x = gm[1,3], y = .08, xend = gm[1,3], yend = 0), colour = "dark grey",
               arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
  geom_segment(aes(x = gm[2,3], y = .08, xend = gm[2,3], yend = 0), colour = "#56326E",
               arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
  geom_segment(aes(x = gm[3,3], y = .08, xend = gm[3,3], yend = 0), colour = "#ED7F6F",
               arrow = arrow(length = unit(0.15, "cm"), type = "closed")) + 
  geom_segment(aes(x = gm[4,3], y = .08, xend = gm[4,3], yend = 0), colour = "#ABA778",
               arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
   scale_x_log10() +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8))

#plot together
p<-ggarrange(all_print, zoom_print,
             #labels = c("", ""),
             ncol = 1, nrow = 2)
#p
#ggsave("LD50_wo_subset.pdf",p, width=9, height=7, units="in")


##do the same, but for normalized n=12 samples per group 
#get median value for each distance
rsq_means_all_n12 <- read.delim("~/bigdata/pop_genomics/LD_decay_4/LD_out/rsq_means_n12_Clade_all.tab",sep="",header=T,check.names=F,stringsAsFactors=F)
rsq_means_Clade1_n12 <- read.delim("~/bigdata/pop_genomics/LD_decay_4/LD_out/rsq_means_n12_Clade_1.tab",sep="",header=T,check.names=F,stringsAsFactors=F)
rsq_means_Clade2_n12 <- read.delim("~/bigdata/pop_genomics/LD_decay_4/LD_out/rsq_means_n12_Clade_2.tab",sep="",header=T,check.names=F,stringsAsFactors=F)
rsq_means_Clade3_n12 <- read.delim("~/bigdata/pop_genomics/LD_decay_4/LD_out/rsq_means_n12_Clade_3.tab",sep="",header=T,check.names=F,stringsAsFactors=F)

#add group annotation 
rsq_means_all_n12$Clade<- "All"
rsq_means_Clade1_n12$Clade<- "Clade1"
rsq_means_Clade2_n12$Clade<- "Clade2"
rsq_means_Clade3_n12$Clade<- "Clade3"

#get intersection of R2 half way between R2max and R2 mean 
#get y w/o averaging
LD50_all_y2<- min(rsq_means_all_n12$rsq.mean) + (max(rsq_means_all_n12$rsq.mean) - min(rsq_means_all_n12$rsq.mean)) /2
LD50_C1_y2<- min(rsq_means_Clade1_n12$rsq.mean) + (max(rsq_means_Clade1_n12$rsq.mean) - min(rsq_means_Clade1_n12$rsq.mean)) /2
LD50_C2_y2<- min(rsq_means_Clade2_n12$rsq.mean) + (max(rsq_means_Clade2_n12$rsq.mean) - min(rsq_means_Clade2_n12$rsq.mean)) /2
LD50_C3_y2<- min(rsq_means_Clade3_n12$rsq.mean) + (max(rsq_means_Clade3_n12$rsq.mean) - min(rsq_means_Clade3_n12$rsq.mean)) /2

#calculate x vals 
LD50_all_x2 <- approx(x = rsq_means_all_n12$rsq.mean, y = rsq_means_all_n12$dist, xout = LD50_all_y2)$y
LD50_all_x2

LD50_C1_x2 <- approx(x = rsq_means_Clade1_n12$rsq.mean, y = rsq_means_Clade1_n12$dist, xout = LD50_C1_y2)$y
LD50_C1_x2

LD50_C2_x2 <- approx(x = rsq_means_Clade2_n12$rsq.mean, y = rsq_means_Clade2_n12$dist, xout = LD50_C2_y2)$y
LD50_C2_x2

LD50_C3_x2 <- approx(x = rsq_means_Clade3_n12$rsq.mean, y = rsq_means_Clade3_n12$dist, xout = LD50_C3_y2)$y
LD50_C3_x2

#concat LD50y, LD50x, and group designations
gm2<- data.frame(Clade = c("All", "Clade1", "Clade2", "Clade3"), 
                LD50y =c(LD50_all_y2, 
                         LD50_C1_y2,
                         LD50_C2_y2,
                         LD50_C3_y2),
                LD50x =c(LD50_all_x2, 
                         LD50_C1_x2,
                         LD50_C2_x2,
                         LD50_C3_x2))


#bind
rsq_means_grand<- data.frame(rbind(rsq_means_Clade1_n12,
                                   rsq_means_Clade2_n12,
                                   rsq_means_Clade3_n12,
                                   rsq_means_all_n12))


nrow(rsq_means_grand) #too many data points for easy plot rendering
#subset for east rendering
random_input<- rsq_means_grand[sample(nrow(rsq_means_grand), 20000),]
nrow(random_input)

##plot all
sp2<- ggplot(data=random_input,aes(x=dist,y=rsq.mean))+
  geom_line(data=random_input, aes(color=Clade),size=0.1,alpha=0.9)+
  geom_smooth(aes(color = Clade), show.legend = FALSE, size=.8) +
  geom_smooth(aes(color = Clade), se=FALSE, size=.8) +
  labs(x="Distance (BP)",y=expression(LD~(r^{2})))+
  theme_bw() + scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8))
all_print2<- sp2 + scale_color_manual(values=myCol)
#ggsave("LD_decay_wLD50_at_n12.pdf", plot = last_plot(), width = 8, height = 4)


##plot zoomed in
#subset to only the first 10,000 bp
rsq_means_all_n12_zoom<- rsq_means_all_n12[1:20000,]
rsq_means_C1_zoom<- rsq_means_Clade1_n12[1:20000,]
rsq_means_C2_zoom<- rsq_means_Clade2_n12[1:20000,]
rsq_means_C3_zoom<- rsq_means_Clade3_n12[1:20000,]


rsq_means_grand_zoom<- data.frame(rbind(rsq_means_C1_zoom,
                                        rsq_means_C2_zoom,
                                        rsq_means_C3_zoom,
                                        rsq_means_all_n12_zoom))


#plot
sp_zoom2<- ggplot(data=rsq_means_grand_zoom,aes(x=dist,y=rsq.mean))+
  geom_line(data=rsq_means_grand_zoom, aes(color=Clade),size=0.1,alpha=0.9)+
#  geom_smooth(aes(color = Clade), show.legend = FALSE) +
#  geom_smooth(aes(color = Clade), se=FALSE) +
  labs(x="Distance (BP) log scale",y=expression(LD~(r^{2})))+
  theme_bw()
zoom_print2<- sp_zoom2 + scale_color_manual(values=myCol) + 
  #add lines to indicate the point at which each Clade is half decayed
  geom_segment(aes(x = gm2[1,3], y = .08, xend = gm2[1,3], yend = 0), colour = "dark grey",
               arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
  geom_segment(aes(x = gm2[2,3], y = .08, xend = gm2[2,3], yend = 0), colour = "#56326E",
               arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
  geom_segment(aes(x = gm2[3,3], y = .08, xend = gm2[3,3], yend = 0), colour = "#ED7F6F",
               arrow = arrow(length = unit(0.15, "cm"), type = "closed")) + 
  geom_segment(aes(x = gm2[4,3], y = .08, xend = gm2[4,3], yend = 0), colour = "#ABA778",
               arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
  scale_x_log10() +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8))

#plot together
p2<-ggarrange(all_print, zoom_print,
              all_print2, zoom_print2,
              #labels = c("", ""),
              ncol = 2, nrow = 4)
#p2
ggsave("LD_decay_at_n12_20000pts.pdf",p2, width=9, height=7, units="in")
