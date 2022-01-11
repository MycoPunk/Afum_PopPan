#this script creates figures of depth profiles over MAT loci for MAT-1 and MAT-2 for the 9 strains that have significant alignments to both idomorphs. 

#set packages 
library(data.table)
library(tidyverse)
library(ggforce)
library(ggimage)
library(ggplot2)
library(seqinr)
library(ShortRead)


setwd("~/Desktop/Project_Afum_pangenome_3/Raw_MAT_alignments")
fl <- list.files(path="MAT1-1/fastas",pattern="*.fa",full.names=T)
fl
#read in MAT1 fastas
MAT1_08_36_03_25<- data.frame(read.fasta(file = "MAT1-1/fastas/08-36-03-25.MAT1.consensus.fa", 
                       seqtype = "DNA",as.string = FALSE, set.attributes = FALSE, forceDNAtolower = TRUE))
MAT1_AF100_12_5<- data.frame(read.fasta(file = "MAT1-1/fastas/AF100-12_5.MAT1.consensus.fa", 
                              seqtype = "DNA",as.string = FALSE, set.attributes = FALSE, forceDNAtolower = TRUE))
MAT1_AF293<- data.frame(read.fasta(file = "MAT1-1/fastas/AF293.MAT1.consensus.fa", 
                             seqtype = "DNA",as.string = FALSE, set.attributes = FALSE, forceDNAtolower = TRUE))
MAT1_AF90<- data.frame(read.fasta(file = "MAT1-1/fastas/AF90.MAT1.consensus.fa", 
                             seqtype = "DNA",as.string = FALSE, set.attributes = FALSE, forceDNAtolower = TRUE))
MAT1_Afu_343_P_11<- data.frame(read.fasta(file = "MAT1-1/fastas/Afu_343-P-11.MAT1.consensus.fa", 
                             seqtype = "DNA",as.string = FALSE, set.attributes = FALSE, forceDNAtolower = TRUE))
MAT1_B7586_CDC_30<- data.frame(read.fasta(file = "MAT1-1/fastas/B7586_CDC-30.MAT1.consensus.fa", 
                             seqtype = "DNA",as.string = FALSE, set.attributes = FALSE, forceDNAtolower = TRUE))
MAT1_DMC2_AF100_1_18<- data.frame(read.fasta(file = "MAT1-1/fastas/DMC2_AF100-1_18.MAT1.consensus.fa", 
                             seqtype = "DNA",as.string = FALSE, set.attributes = FALSE, forceDNAtolower = TRUE))
MAT1_DMC2_AF100_1_3<- data.frame(read.fasta(file = "MAT1-1/fastas/DMC2_AF100-1_3.MAT1.consensus.fa", 
                             seqtype = "DNA",as.string = FALSE, set.attributes = FALSE, forceDNAtolower = TRUE))
MAT1_IFM_59359<- data.frame(read.fasta(file = "MAT1-1/fastas/IFM_59359.MAT1.consensus.fa", 
                             seqtype = "DNA",as.string = FALSE, set.attributes = FALSE, forceDNAtolower = TRUE))
MAT1_IFM_61407<- data.frame(read.fasta(file = "MAT1-1/fastas/IFM_61407.MAT1.consensus.fa", 
                             seqtype = "DNA",as.string = FALSE, set.attributes = FALSE, forceDNAtolower = TRUE))
MAT1_NCPF_7816<- data.frame(read.fasta(file = "MAT1-1/fastas/NCPF-7816.MAT1.consensus.fa", 
                             seqtype = "DNA",as.string = FALSE, set.attributes = FALSE, forceDNAtolower = TRUE))
MAT1_AF100_12_7G<- data.frame(read.fasta(file = "MAT1-1/fastas/AF100-12_7G.MAT1.consensus.fa", 
                             seqtype = "DNA",as.string = FALSE, set.attributes = FALSE, forceDNAtolower = TRUE))
MAT1_SF2S9<- data.frame(read.fasta(file = "MAT1-1/fastas/SF2S9.MAT1.consensus.fa", 
                             seqtype = "DNA",as.string = FALSE, set.attributes = FALSE, forceDNAtolower = TRUE))



#read in MAT2 fastas
MAT2_08_36_03_25<- data.frame(read.fasta(file = "MAT1-2/fastas/08-36-03-25.MAT2.consensus.fa", 
                              seqtype = "DNA",as.string = FALSE, set.attributes = FALSE, forceDNAtolower = TRUE))
MAT2_AF100_12_5<- data.frame(read.fasta(file = "MAT1-2/fastas/AF100-12_5.MAT2.consensus.fa", 
                             seqtype = "DNA",as.string = FALSE, set.attributes = FALSE, forceDNAtolower = TRUE))
MAT2_AF293<- data.frame(read.fasta(file = "MAT1-2/fastas/AF293.MAT2.consensus.fa", 
                        seqtype = "DNA",as.string = FALSE, set.attributes = FALSE, forceDNAtolower = TRUE))
MAT2_AF90<- data.frame(read.fasta(file = "MAT1-2/fastas/AF90.MAT2.consensus.fa", 
                       seqtype = "DNA",as.string = FALSE, set.attributes = FALSE, forceDNAtolower = TRUE))
MAT2_Afu_343_P_11<- data.frame(read.fasta(file = "MAT1-2/fastas/Afu_343-P-11.MAT2.consensus.fa", 
                               seqtype = "DNA",as.string = FALSE, set.attributes = FALSE, forceDNAtolower = TRUE))
MAT2_B7586_CDC_30<- data.frame(read.fasta(file = "MAT1-2/fastas/B7586_CDC-30.MAT2.consensus.fa", 
                               seqtype = "DNA",as.string = FALSE, set.attributes = FALSE, forceDNAtolower = TRUE))
MAT2_DMC2_AF100_1_18<- data.frame(read.fasta(file = "MAT1-2/fastas/DMC2_AF100-1_18.MAT2.consensus.fa", 
                                 seqtype = "DNA",as.string = FALSE, set.attributes = FALSE, forceDNAtolower = TRUE))
MAT2_DMC2_AF100_1_3<- data.frame(read.fasta(file = "MAT1-2/fastas/DMC2_AF100-1_3.MAT2.consensus.fa", 
                                 seqtype = "DNA",as.string = FALSE, set.attributes = FALSE, forceDNAtolower = TRUE))
MAT2_IFM_59359<- data.frame(read.fasta(file = "MAT1-2/fastas/IFM_59359.MAT2.consensus.fa", 
                            seqtype = "DNA",as.string = FALSE, set.attributes = FALSE, forceDNAtolower = TRUE))
MAT2_IFM_61407<- data.frame(read.fasta(file = "MAT1-2/fastas/IFM_61407.MAT2.consensus.fa", 
                            seqtype = "DNA",as.string = FALSE, set.attributes = FALSE, forceDNAtolower = TRUE))
MAT2_NCPF_7816<- data.frame(read.fasta(file = "MAT1-2/fastas/NCPF-7816.MAT2.consensus.fa", 
                            seqtype = "DNA",as.string = FALSE, set.attributes = FALSE, forceDNAtolower = TRUE))
MAT2_AF100_12_7G<- data.frame(read.fasta(file = "MAT1-2/fastas/AF100-12_7G.MAT2.consensus.fa", 
                            seqtype = "DNA",as.string = FALSE, set.attributes = FALSE, forceDNAtolower = TRUE))
MAT2_SF2S9<- data.frame(read.fasta(file = "MAT1-2/fastas/SF2S9.MAT2.consensus.fa", 
                            seqtype = "DNA",as.string = FALSE, set.attributes = FALSE, forceDNAtolower = TRUE))



#read in all depth files
dp <- list.files(path="MAT1-1/depth",pattern="*",full.names=T)
dp
MAT1_08_36_03_25_depth<- read.table("MAT1-1/depth/08-36-03-25.MAT1.coverage", header = FALSE, sep = "", colClasses=c("NULL", "NULL", NA))
MAT1_AF100_12_5_depth<- read.table("MAT1-1/depth/AF100-12_5.MAT1.coverage", header = FALSE, sep = "", colClasses=c("NULL", "NULL", NA))
MAT1_AF293_depth<- read.table("MAT1-1/depth/AF293.MAT1.coverage", header = FALSE, sep = "", colClasses=c("NULL", "NULL", NA))
MAT1_AF90_depth<- read.table("MAT1-1/depth/AF90.MAT1.coverage", header = FALSE, sep = "", colClasses=c("NULL", "NULL", NA))
MAT1_Afu_343_P_11_depth<- read.table("MAT1-1/depth/Afu_343-P-11.MAT1.coverage", header = FALSE, sep = "", colClasses=c("NULL", "NULL", NA))
MAT1_B7586_CDC_30_depth<- read.table("MAT1-1/depth/B7586_CDC-30.MAT1.coverage", header = FALSE, sep = "", colClasses=c("NULL", "NULL", NA))
MAT1_DMC2_AF100_1_18_depth<- read.table("MAT1-1/depth/DMC2_AF100-1_18.MAT1.coverage", header = FALSE, sep = "", colClasses=c("NULL", "NULL", NA))
MAT1_DMC2_AF100_1_3_depth<- read.table("MAT1-1/depth/DMC2_AF100-1_3.MAT1.coverage", header = FALSE, sep = "", colClasses=c("NULL", "NULL", NA))
MAT1_IFM_59359_depth<- read.table("MAT1-1/depth/IFM_59359.MAT1.coverage", header = FALSE, sep = "", colClasses=c("NULL", "NULL", NA))
MAT1_IFM_61407_depth<- read.table("MAT1-1/depth/IFM_61407.MAT1.coverage", header = FALSE, sep = "", colClasses=c("NULL", "NULL", NA))
MAT1_NCPF_7816_depth<- read.table("MAT1-1/depth/NCPF-7816.MAT1.coverage", header = FALSE, sep = "", colClasses=c("NULL", "NULL", NA))
MAT1_SF2S9_depth<- read.table("MAT1-1/depth/SF2S9.MAT1.coverage", header = FALSE, sep = "", colClasses=c("NULL", "NULL", NA))
MAT1_AF100_12_7G_depth<- read.table("MAT1-1/depth/AF100-12_7G.MAT1.coverage", header = FALSE, sep = "", colClasses=c("NULL", "NULL", NA))


MAT2_08_36_03_25_depth<- read.table("MAT1-2/depth/08-36-03-25.MAT2.coverage", header = FALSE, sep = "", colClasses=c("NULL", "NULL", NA))
MAT2_AF100_12_5_depth<- read.table("MAT1-2/depth/AF100-12_5.MAT2.coverage", header = FALSE, sep = "", colClasses=c("NULL", "NULL", NA))
MAT2_AF293_depth<- read.table("MAT1-2/depth/AF293.MAT2.coverage", header = FALSE, sep = "", colClasses=c("NULL", "NULL", NA))
MAT2_AF90_depth<- read.table("MAT1-2/depth/AF90.MAT2.coverage", header = FALSE, sep = "", colClasses=c("NULL", "NULL", NA))
MAT2_Afu_343_P_11_depth<- read.table("MAT1-2/depth/Afu_343-P-11.MAT2.coverage", header = FALSE, sep = "", colClasses=c("NULL", "NULL", NA))
MAT2_B7586_CDC_30_depth<- read.table("MAT1-2/depth/B7586_CDC-30.MAT2.coverage", header = FALSE, sep = "", colClasses=c("NULL", "NULL", NA))
MAT2_DMC2_AF100_1_18_depth<- read.table("MAT1-2/depth/DMC2_AF100-1_18.MAT2.coverage", header = FALSE, sep = "", colClasses=c("NULL", "NULL", NA))
MAT2_DMC2_AF100_1_3_depth<- read.table("MAT1-2/depth/DMC2_AF100-1_3.MAT2.coverage", header = FALSE, sep = "", colClasses=c("NULL", "NULL", NA))
MAT2_IFM_59359_depth<- read.table("MAT1-2/depth/IFM_59359.MAT2.coverage", header = FALSE, sep = "", colClasses=c("NULL", "NULL", NA))
MAT2_IFM_61407_depth<- read.table("MAT1-2/depth/IFM_61407.MAT2.coverage", header = FALSE, sep = "", colClasses=c("NULL", "NULL", NA))
MAT2_NCPF_7816_depth<- read.table("MAT1-2/depth/NCPF-7816.MAT2.coverage", header = FALSE, sep = "", colClasses=c("NULL", "NULL", NA))
MAT2_SF2S9_depth<- read.table("MAT1-2/depth/SF2S9.MAT2.coverage", header = FALSE, sep = "", colClasses=c("NULL", "NULL", NA))
MAT2_AF100_12_7G_depth<- read.table("MAT1-2/depth/AF100-12_7G.MAT2.coverage", header = FALSE, sep = "", colClasses=c("NULL", "NULL", NA))


#bind
MAT1_08_36_03_25_all<- data.frame(cbind("seq" = MAT1_08_36_03_25[,1], "depth" = as.numeric(MAT1_08_36_03_25_depth$V3)))
MAT1_AF100_12_5_all<- data.frame(cbind("seq" = MAT1_AF100_12_5[,1], "depth" = as.numeric(MAT1_AF100_12_5_depth$V3)))
MAT1_AF293_all<- data.frame(cbind("seq" = MAT1_AF293[,1], "depth" = as.numeric(MAT1_AF293_depth$V3)))
MAT1_AF90_all<- data.frame(cbind("seq" = MAT1_AF90[,1], "depth" = as.numeric(MAT1_AF90_depth$V3)))
MAT1_Afu_343_P_11_all<- data.frame(cbind("seq" = MAT1_Afu_343_P_11[,1], "depth" = as.numeric(MAT1_Afu_343_P_11_depth$V3)))
MAT1_B7586_CDC_30_all<- data.frame(cbind("seq" = MAT1_B7586_CDC_30[,1], "depth" = as.numeric(MAT1_B7586_CDC_30_depth$V3)))
MAT1_DMC2_AF100_1_18_all<- data.frame(cbind("seq" = MAT1_DMC2_AF100_1_18[,1], "depth" = as.numeric(MAT1_DMC2_AF100_1_18_depth$V3)))
MAT1_DMC2_AF100_1_3_all<- data.frame(cbind("seq" = MAT1_DMC2_AF100_1_3[,1], "depth" = as.numeric(MAT1_DMC2_AF100_1_3_depth$V3)))
MAT1_IFM_59359_all<- data.frame(cbind("seq" = MAT1_IFM_59359[,1], "depth" = as.numeric(MAT1_IFM_59359_depth$V3)))
MAT1_IFM_61407_all<- data.frame(cbind("seq" = MAT1_IFM_61407[,1], "depth" = as.numeric(MAT1_IFM_61407_depth$V3)))
MAT1_NCPF_7816_all<- data.frame(cbind("seq" = MAT1_NCPF_7816[,1], "depth" = as.numeric(MAT1_NCPF_7816_depth$V3)))
MAT1_SF2S9_all<- data.frame(cbind("seq" = MAT1_SF2S9[,1], "depth" = as.numeric(MAT1_SF2S9_depth$V3)))
MAT1_AF100_12_7G_all<- data.frame(cbind("seq" = MAT1_AF100_12_7G[,1], "depth" = as.numeric(MAT1_AF100_12_7G_depth$V3)))


MAT2_08_36_03_25_all<- data.frame(cbind("seq" = MAT2_08_36_03_25[,1], "depth" = as.numeric(MAT2_08_36_03_25_depth$V3)))
MAT2_AF100_12_5_all<- data.frame(cbind("seq" = MAT2_AF100_12_5[,1], "depth" = as.numeric(MAT2_AF100_12_5_depth$V3)))
MAT2_AF293_all<- data.frame(cbind("seq" = MAT2_AF293[,1], "depth" = as.numeric(MAT2_AF293_depth$V3)))
MAT2_AF90_all<- data.frame(cbind("seq" = MAT2_AF90[,1], "depth" = as.numeric(MAT2_AF90_depth$V3)))
MAT2_Afu_343_P_11_all<- data.frame(cbind("seq" = MAT2_Afu_343_P_11[,1], "depth" = as.numeric(MAT2_Afu_343_P_11_depth$V3)))
MAT2_B7586_CDC_30_all<- data.frame(cbind("seq" = MAT2_B7586_CDC_30[,1], "depth" = as.numeric(MAT2_B7586_CDC_30_depth$V3)))
MAT2_DMC2_AF100_1_18_all<- data.frame(cbind("seq" = MAT2_DMC2_AF100_1_18[,1], "depth" = as.numeric(MAT2_DMC2_AF100_1_18_depth$V3)))
MAT2_DMC2_AF100_1_3_all<- data.frame(cbind("seq" = MAT2_DMC2_AF100_1_3[,1], "depth" = as.numeric(MAT2_DMC2_AF100_1_3_depth$V3)))
MAT2_IFM_59359_all<- data.frame(cbind("seq" = MAT2_IFM_59359[,1], "depth" = as.numeric(MAT2_IFM_59359_depth$V3)))
MAT2_IFM_61407_all<- data.frame(cbind("seq" = MAT2_IFM_61407[,1], "depth" = as.numeric(MAT2_IFM_61407_depth$V3)))
MAT2_NCPF_7816_all<- data.frame(cbind("seq" = MAT2_NCPF_7816[,1], "depth" = as.numeric(MAT2_NCPF_7816_depth$V3)))
MAT2_SF2S9_all<- data.frame(cbind("seq" = MAT2_SF2S9[,1], "depth" = as.numeric(MAT2_SF2S9_depth$V3)))
MAT2_AF100_12_7G_all<- data.frame(cbind("seq" = MAT2_AF100_12_7G[,1], "depth" = as.numeric(MAT2_AF100_12_7G_depth$V3)))


#graph
MAT1_08_36_03_25_all$pos = as.numeric(rownames(MAT1_08_36_03_25_all))
MAT1_08_36_03_25_all$depth<- as.numeric(MAT1_08_36_03_25_all$depth)
MAT2_08_36_03_25_all$pos = as.numeric(rownames(MAT2_08_36_03_25_all))
MAT2_08_36_03_25_all$depth<- as.numeric(MAT2_08_36_03_25_all$depth)

plot08_36_03_25<- ggplot(MAT1_08_36_03_25_all, aes(x=pos, y=depth)) + 
  geom_point(shape=20, size=1, aes(color = "grey")) + theme_classic() +
  #add BPs
  #scale_x_discrete(limits = as.factor(MAT1_08_36_03_25_all$pos), 
   #                labels = MAT1_08_36_03_25_all$seq) +
  geom_point(data = MAT2_08_36_03_25_all, shape=20, size=1, aes(color = "red"))+
  scale_colour_manual(values=c("red", "grey"),
                      labels = c('MAT1', 'MAT2'))
  
plot08_36_03_25_final<- plot08_36_03_25+ theme(aspect.ratio=1/8, axis.text=element_text(size=6), legend.title = element_blank()) + ggtitle("Afum_08_36_03_25")
plot08_36_03_25_final

#AF293
MAT1_AF293_all$pos = as.numeric(rownames(MAT1_AF293_all))
MAT1_AF293_all$depth<- as.numeric(MAT1_AF293_all$depth)
MAT2_AF293_all$pos = as.numeric(rownames(MAT2_AF293_all))
MAT2_AF293_all$depth<- as.numeric(MAT2_AF293_all$depth)

plotAF293<- ggplot(MAT1_AF293_all, aes(x=pos, y=depth)) + 
  geom_point(shape=20, size=1, aes(color = "grey")) + theme_classic() +
  #add BPs
  #scale_x_discrete(limits = as.factor(MAT1_AF293_all$pos), 
  #                labels = MAT1_AF293_all$seq) +
  geom_point(data = MAT2_AF293_all, shape=20, size=1, aes(color = "red"))+
  scale_colour_manual(values=c("red", "grey"),
                      labels = c('MAT1', 'MAT2'))

plotAF293_final<- plotAF293+ theme(aspect.ratio=1/8, axis.text=element_text(size=6), legend.title = element_blank()) + ggtitle("Afum_AF293")
plotAF293_final

#AF90
MAT1_AF90_all$pos = as.numeric(rownames(MAT1_AF90_all))
MAT1_AF90_all$depth<- as.numeric(MAT1_AF90_all$depth)
MAT2_AF90_all$pos = as.numeric(rownames(MAT2_AF90_all))
MAT2_AF90_all$depth<- as.numeric(MAT2_AF90_all$depth)

plotAF90<- ggplot(MAT1_AF90_all, aes(x=pos, y=depth)) + 
  geom_point(shape=20, size=1, aes(color = "grey")) + theme_classic() +
  #add BPs
  #scale_x_discrete(limits = as.factor(MAT1_AF90_all$pos), 
  #                labels = MAT1_AF90_all$seq) +
  geom_point(data = MAT2_AF90_all, shape=20, size=1, aes(color = "red"))+
  scale_colour_manual(values=c("red", "grey"),
                      labels = c('MAT1', 'MAT2'))

plotAF90_final<- plotAF90+ theme(aspect.ratio=1/8, axis.text=element_text(size=6), legend.title = element_blank()) + ggtitle("Afum_AF90")
plotAF90_final

#AF100_12_5
MAT1_AF100_12_5_all$pos = as.numeric(rownames(MAT1_AF100_12_5_all))
MAT1_AF100_12_5_all$depth<- as.numeric(MAT1_AF100_12_5_all$depth)
MAT2_AF100_12_5_all$pos = as.numeric(rownames(MAT2_AF100_12_5_all))
MAT2_AF100_12_5_all$depth<- as.numeric(MAT2_AF100_12_5_all$depth)

plotAF100_12_5<- ggplot(MAT1_AF100_12_5_all, aes(x=pos, y=depth)) + 
  geom_point(shape=20, size=1, aes(color = "grey")) + theme_classic() +
  #add BPs
  #scale_x_discrete(limits = as.factor(MAT1_AF100_12_5_all$pos), 
  #                labels = MAT1_AF100_12_5_all$seq) +
  geom_point(data = MAT2_AF100_12_5_all, shape=20, size=1, aes(color = "red"))+
  scale_colour_manual(values=c("red", "grey"),
                      labels = c('MAT1', 'MAT2'))

plotAF100_12_5_final<- plotAF100_12_5+ theme(aspect.ratio=1/8, axis.text=element_text(size=6), legend.title = element_blank()) + ggtitle("Afum_AF100_12_5")
plotAF100_12_5_final

#Afu_343_P_11
MAT1_Afu_343_P_11_all$pos = as.numeric(rownames(MAT1_Afu_343_P_11_all))
MAT1_Afu_343_P_11_all$depth<- as.numeric(MAT1_Afu_343_P_11_all$depth)
MAT2_Afu_343_P_11_all$pos = as.numeric(rownames(MAT2_Afu_343_P_11_all))
MAT2_Afu_343_P_11_all$depth<- as.numeric(MAT2_Afu_343_P_11_all$depth)

plotAfu_343_P_11<- ggplot(MAT1_Afu_343_P_11_all, aes(x=pos, y=depth)) + 
  geom_point(shape=20, size=1, aes(color = "grey")) + theme_classic() +
  #add BPs
  #scale_x_discrete(limits = as.factor(MAT1_Afu_343_P_11_all$pos), 
  #                labels = MAT1_Afu_343_P_11_all$seq) +
  geom_point(data = MAT2_Afu_343_P_11_all, shape=20, size=1, aes(color = "red"))+
  scale_colour_manual(values=c("red", "grey"),
                      labels = c('MAT1', 'MAT2'))

plotAfu_343_P_11_final<- plotAfu_343_P_11+ theme(aspect.ratio=1/8, axis.text=element_text(size=6), legend.title = element_blank()) + ggtitle("Afum_Afu_343_P_11")
plotAfu_343_P_11_final

#B7586_CDC
MAT1_B7586_CDC_30_all$pos = as.numeric(rownames(MAT1_B7586_CDC_30_all))
MAT1_B7586_CDC_30_all$depth<- as.numeric(MAT1_B7586_CDC_30_all$depth)
MAT2_B7586_CDC_30_all$pos = as.numeric(rownames(MAT2_B7586_CDC_30_all))
MAT2_B7586_CDC_30_all$depth<- as.numeric(MAT2_B7586_CDC_30_all$depth)

plotB7586_CDC_30<- ggplot(MAT1_B7586_CDC_30_all, aes(x=pos, y=depth)) + 
  geom_point(shape=20, size=1, aes(color = "grey")) + theme_classic() +
  #add BPs
  #scale_x_discrete(limits = as.factor(MAT1_B7586_CDC_30_all$pos), 
  #                labels = MAT1_B7586_CDC_30_all$seq) +
  geom_point(data = MAT2_B7586_CDC_30_all, shape=20, size=1, aes(color = "red"))+
  scale_colour_manual(values=c("red", "grey"),
                      labels = c('MAT1', 'MAT2'))

plotB7586_CDC_30_final<- plotB7586_CDC_30+ theme(aspect.ratio=1/8, axis.text=element_text(size=6), legend.title = element_blank()) + ggtitle("Afum_B7586_CDC_30")
plotB7586_CDC_30_final

#DMC2_AF100_1_18
MAT1_DMC2_AF100_1_18_all$pos = as.numeric(rownames(MAT1_DMC2_AF100_1_18_all))
MAT1_DMC2_AF100_1_18_all$depth<- as.numeric(MAT1_DMC2_AF100_1_18_all$depth)
MAT2_DMC2_AF100_1_18_all$pos = as.numeric(rownames(MAT2_DMC2_AF100_1_18_all))
MAT2_DMC2_AF100_1_18_all$depth<- as.numeric(MAT2_DMC2_AF100_1_18_all$depth)

plotDMC2_AF100_1_18<- ggplot(MAT1_DMC2_AF100_1_18_all, aes(x=pos, y=depth)) + 
  geom_point(shape=20, size=1, aes(color = "grey")) + theme_classic() +
  #add BPs
  #scale_x_discrete(limits = as.factor(MAT1_DMC2_AF100_1_18_all$pos), 
  #                labels = MAT1_DMC2_AF100_1_18_all$seq) +
  geom_point(data = MAT2_DMC2_AF100_1_18_all, shape=20, size=1, aes(color = "red"))+
  scale_colour_manual(values=c("red", "grey"),
                      labels = c('MAT1', 'MAT2'))

plotDMC2_AF100_1_18_final<- plotDMC2_AF100_1_18+ theme(aspect.ratio=1/8, axis.text=element_text(size=6), legend.title = element_blank()) + ggtitle("Afum_AF100_1_18")
plotDMC2_AF100_1_18_final

#AF100_1_3
MAT1_DMC2_AF100_1_3_all$pos = as.numeric(rownames(MAT1_DMC2_AF100_1_3_all))
MAT1_DMC2_AF100_1_3_all$depth<- as.numeric(MAT1_DMC2_AF100_1_3_all$depth)
MAT2_DMC2_AF100_1_3_all$pos = as.numeric(rownames(MAT2_DMC2_AF100_1_3_all))
MAT2_DMC2_AF100_1_3_all$depth<- as.numeric(MAT2_DMC2_AF100_1_3_all$depth)

plotDMC2_AF100_1_3<- ggplot(MAT1_DMC2_AF100_1_3_all, aes(x=pos, y=depth)) + 
  geom_point(shape=20, size=1, aes(color = "grey")) + theme_classic() +
  #add BPs
  #scale_x_discrete(limits = as.factor(MAT1_DMC2_AF100_1_3_all$pos), 
  #                labels = MAT1_DMC2_AF100_1_3_all$seq) +
  geom_point(data = MAT2_DMC2_AF100_1_3_all, shape=20, size=1, aes(color = "red"))+
  scale_colour_manual(values=c("red", "grey"),
                      labels = c('MAT1', 'MAT2'))

plotDMC2_AF100_1_3_final<- plotDMC2_AF100_1_3+ theme(aspect.ratio=1/8, axis.text=element_text(size=6), legend.title = element_blank()) + ggtitle("Afum_AF100_1_3")
plotDMC2_AF100_1_3_final

#IFM_59359
MAT1_IFM_59359_all$pos = as.numeric(rownames(MAT1_IFM_59359_all))
MAT1_IFM_59359_all$depth<- as.numeric(MAT1_IFM_59359_all$depth)
MAT2_IFM_59359_all$pos = as.numeric(rownames(MAT2_IFM_59359_all))
MAT2_IFM_59359_all$depth<- as.numeric(MAT2_IFM_59359_all$depth)

plotIFM_59359<- ggplot(MAT1_IFM_59359_all, aes(x=pos, y=depth)) + 
  geom_point(shape=20, size=1, aes(color = "grey")) + theme_classic() +
  #add BPs
  #scale_x_discrete(limits = as.factor(MAT1_IFM_59359_all$pos), 
  #                labels = MAT1_IFM_59359_all$seq) +
  geom_point(data = MAT2_IFM_59359_all, shape=20, size=1, aes(color = "red"))+
  scale_colour_manual(values=c("red", "grey"),
                      labels = c('MAT1', 'MAT2'))

plotIFM_59359_final<- plotIFM_59359+ theme(aspect.ratio=1/8, axis.text=element_text(size=6), legend.title = element_blank()) + ggtitle("Afum_IFM_59359")
plotIFM_59359_final

#IFM_61407
MAT1_IFM_61407_all$pos = as.numeric(rownames(MAT1_IFM_61407_all))
MAT1_IFM_61407_all$depth<- as.numeric(MAT1_IFM_61407_all$depth)
MAT2_IFM_61407_all$pos = as.numeric(rownames(MAT2_IFM_61407_all))
MAT2_IFM_61407_all$depth<- as.numeric(MAT2_IFM_61407_all$depth)

plotIFM_61407<- ggplot(MAT1_IFM_61407_all, aes(x=pos, y=depth)) + 
  geom_point(shape=20, size=1, aes(color = "grey")) + theme_classic() +
  #add BPs
  #scale_x_discrete(limits = as.factor(MAT1_IFM_61407_all$pos), 
  #                labels = MAT1_IFM_61407_all$seq) +
  geom_point(data = MAT2_IFM_61407_all, shape=20, size=1, aes(color = "red"))+
  scale_colour_manual(values=c("red", "grey"),
                      labels = c('MAT1', 'MAT2'))

plotIFM_61407_final<- plotIFM_61407+ theme(aspect.ratio=1/8, axis.text=element_text(size=6), legend.title = element_blank()) + ggtitle("Afum_IFM_61407")
plotIFM_61407_final

#NCPF_7816
MAT1_NCPF_7816_all$pos = as.numeric(rownames(MAT1_NCPF_7816_all))
MAT1_NCPF_7816_all$depth<- as.numeric(MAT1_NCPF_7816_all$depth)
MAT2_NCPF_7816_all$pos = as.numeric(rownames(MAT2_NCPF_7816_all))
MAT2_NCPF_7816_all$depth<- as.numeric(MAT2_NCPF_7816_all$depth)

plotNCPF_7816<- ggplot(MAT1_NCPF_7816_all, aes(x=pos, y=depth)) + 
  geom_point(shape=20, size=1, aes(color = "grey")) + theme_classic() +
  #add BPs
  #scale_x_discrete(limits = as.factor(MAT1_NCPF_7816_all$pos), 
  #                labels = MAT1_NCPF_7816_all$seq) +
  geom_point(data = MAT2_NCPF_7816_all, shape=20, size=1, aes(color = "red"))+
  scale_colour_manual(values=c("red", "grey"),
                      labels = c('MAT1', 'MAT2'))

plotNCPF_7816_final<- plotNCPF_7816+ theme(aspect.ratio=1/8, axis.text=element_text(size=6), legend.title = element_blank()) + ggtitle("Afum_NCPF_7816")
plotNCPF_7816_final

#SF2S9
MAT1_SF2S9_all$pos = as.numeric(rownames(MAT1_SF2S9_all))
MAT1_SF2S9_all$depth<- as.numeric(MAT1_SF2S9_all$depth)
MAT2_SF2S9_all$pos = as.numeric(rownames(MAT2_SF2S9_all))
MAT2_SF2S9_all$depth<- as.numeric(MAT2_SF2S9_all$depth)

plotSF2S9<- ggplot(MAT1_SF2S9_all, aes(x=pos, y=depth)) + 
  geom_point(shape=20, size=1, aes(color = "grey")) + theme_classic() +
  #add BPs
  #scale_x_discrete(limits = as.factor(MAT1_SF2S9_all$pos), 
  #                labels = MAT1_SF2S9_all$seq) +
  geom_point(data = MAT2_SF2S9_all, shape=20, size=1, aes(color = "red"))+
  scale_colour_manual(values=c("red", "grey"),
                      labels = c('MAT1', 'MAT2'))

plotSF2S9_final<- plotSF2S9+ theme(aspect.ratio=1/8, axis.text=element_text(size=6), legend.title = element_blank()) + ggtitle("Afum_SF2S9")
plotSF2S9_final

#AF100_12_7G
MAT1_AF100_12_7G_all$pos = as.numeric(rownames(MAT1_AF100_12_7G_all))
MAT1_AF100_12_7G_all$depth<- as.numeric(MAT1_AF100_12_7G_all$depth)
MAT2_AF100_12_7G_all$pos = as.numeric(rownames(MAT2_AF100_12_7G_all))
MAT2_AF100_12_7G_all$depth<- as.numeric(MAT2_AF100_12_7G_all$depth)

plotAF100_12_7G<- ggplot(MAT1_AF100_12_7G_all, aes(x=pos, y=depth)) + 
  geom_point(shape=20, size=1, aes(color = "grey")) + theme_classic() +
  #add BPs
  #scale_x_discrete(limits = as.factor(MAT1_AF100_12_7G_all$pos), 
  #                labels = MAT1_AF100_12_7G_all$seq) +
  geom_point(data = MAT2_AF100_12_7G_all, shape=20, size=1, aes(color = "red"))+
  scale_colour_manual(values=c("red", "grey"),
                      labels = c('MAT1', 'MAT2'))

plotAF100_12_7G_final<- plotAF100_12_7G+ theme(aspect.ratio=1/8, axis.text=element_text(size=6), legend.title = element_blank()) + ggtitle("Afum_AF100_12_7G")
plotAF100_12_7G_final




#install.packages("ggpubr")
library(ggpubr)
#plot all 
p<-ggarrange(plotAF90_final, plotAF293_final, 
          plot08_36_03_25_final,
          plotNCPF_7816_final,
          plotAF100_12_5_final,
          plotAfu_343_P_11_final,
          plotB7586_CDC_30_final,
          plotDMC2_AF100_1_18_final,
          plotDMC2_AF100_1_3_final,
          plotIFM_59359_final,
          plotIFM_61407_final,
          plotSF2S9_final,
          plotAF100_12_7G_final,
labels = c("MAT-1 control", "MAT-2 control", "", "", "", "", ""),
          ncol = 2, nrow = 7)
#p
ggsave("MAT_alignments_11.Jan.2022.pdf",p, width=10, height=9, units="in")

