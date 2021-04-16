##this script runs allele variation analysis to determine support for diploidy for nine A.fumigatus strains with full length alignments to both MAT1 and MAT2
##Last updated 16.Apr.2021

##set wd
setwd("~/Desktop/Project_Afum_pangenome/")

##load libraries
library('vcfR')
library('pinfsc50')
library('ggpubr')
library('gridGraphics')
library('grid')


#sessionInfo() 

##set seed for reproducibility
set.seed(666)

##read in data
vcf<- read.vcfR("262_strains.selected.SNP.NO_TEs.vcf.gz") #W/o TE's and w 262 strains
Afum_grp<-read.delim("grp_temp.txt", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE)

knitr::kable(vcf@gt[c(1:2,11,30),1:4])

#depth each allele was sequence at 
ad <- extract.gt(vcf, element = 'AD')

#use function mansplit to extract first and second allele
allele1 <- masplit(ad, record = 1)
allele2 <- masplit(ad, record = 2)
ad1 <- allele1 / (allele1 + allele2)
ad2 <- allele2 / (allele1 + allele2)

#plot allele frequencies in hist.

#first AF293 control:
pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"AF293"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "AF293", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"AF293"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
Afum_AF293_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
Afum_AF293_all
#very clean, either matches ref or does not- haploid. 


pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"08-36-03-25"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "08_36_03_25", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"08-36-03-25"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
Afum_08_36_03_25_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
Afum_08_36_03_25_all
#little shoulders - looks like noise?

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"AF100-12_5"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "Afum_AF100_12_5", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"AF100-12_5"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
Afum_AF100_12_5_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
Afum_AF100_12_5_all
#even smaller shoulders- looks like noise.

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"Afu_343-P/11"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "Afu_343_P_11", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"Afu_343-P/11"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
Afum_Afu_343_P_11_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
Afum_Afu_343_P_11_all
#shoulders that look like noise, but slightly further in - no near 1/2 though. 

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"B7586_CDC-30"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "Afu_B7586_CDC_30", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"B7586_CDC-30"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
Afum_Afu_B7586_CDC_30_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
Afum_Afu_B7586_CDC_30_all
#very clean, almost no shoulders - haploid

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"DMC2_AF100-1_18"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "Afum_AF100_1_18", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"DMC2_AF100-1_18"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
Afum_DMC2_AF100_1_18_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
Afum_DMC2_AF100_1_18_all
#very clean, only very slight shoulders, no where hear 1/2 - haploid

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"DMC2_AF100-1_3"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "Afum_AF100_1_3", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"DMC2_AF100-1_3"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
Afum_DMC2_AF100_1_3_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
Afum_DMC2_AF100_1_3_all
#slight bump at 1/2 
#trim and investigate!

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"IFM_59359"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "IFM_59359", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"IFM_59359"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
Afum_IFM_59359_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
Afum_IFM_59359_all
#slight bump at 1/2 
#trim and investigate!

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"IFM_61407"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "IFM_61407", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"IFM_61407"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
Afum_IFM_61407_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
Afum_IFM_61407_all
#slight bump at 1/2 
#trim and investigate!

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2[,"NCPF-7816"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", main = "NCPF_7816", xlab ="allele frequency", ylab = "depth frequency")
hist(ad1[,"NCPF-7816"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
Afum_NCPF_7816_all <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
Afum_NCPF_7816_all
#slight shoulders



#remove homozygotes that overwhelm the plot, to focus on heterozygotes.
ad2_df<- data.frame(ad2)
ad2_df[ad2_df ==0] <- NA
ad1_df<- data.frame(ad1)
ad1_df[ad1_df ==1] <- NA

pdf(NULL)
dev.control(displaylist="enable")
#first AF293 control:
hist(ad2_df[,"AF293"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", xlab =NA, ylab = NA)
hist(ad1_df[,"AF293"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
Afum_AF293_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
Afum_AF293_zoom 
#very low frequency outside of 0/1

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2_df[,"X08.36.03.25"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n",  xlab =NA, ylab = NA)
hist(ad1_df[,"X08.36.03.25"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
Afum_08_36_03_25_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
Afum_08_36_03_25_zoom 
#shoulders - still look like noise.

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2_df[,"AF100.12_5"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", xlab =NA, ylab = NA)
hist(ad1_df[,"AF100.12_5"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
Afum_AF100_12_5_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
Afum_AF100_12_5_zoom
#looks like noise.

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2_df[,"Afu_343.P.11"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", xlab =NA, ylab = NA)
hist(ad1_df[,"Afu_343.P.11"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
Afum_Afu_343_P_11_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
Afum_Afu_343_P_11_zoom
#shoulders that look like noise, but slightly further in - no near 1/2 though. 

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2_df[,"B7586_CDC.30"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", xlab =NA, ylab = NA)
hist(ad1_df[,"B7586_CDC.30"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
Afum_B7586_CDC_30_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
Afum_B7586_CDC_30_zoom
#very clean, almost no shoulders - haploid

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2_df[,"DMC2_AF100.1_18"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", xlab =NA, ylab = NA)
hist(ad1_df[,"DMC2_AF100.1_18"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
Afum_DMC2_AF100_1_18_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
Afum_DMC2_AF100_1_18_zoom
#very clean, only very slight shoulders

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2_df[,"DMC2_AF100.1_3"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", xlab =NA, ylab = NA)
hist(ad1_df[,"DMC2_AF100.1_3"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
Afum_DMC2_AF100_1_3_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
Afum_DMC2_AF100_1_3_zoom
#slight bump at 1/2 

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2_df[,"IFM_59359"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", xlab =NA, ylab = NA)
hist(ad1_df[,"IFM_59359"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
Afum_IFM_59359_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
Afum_IFM_59359_zoom
#slight bump at 1/2 

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2_df[,"IFM_61407"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", xlab =NA, ylab = NA)
hist(ad1_df[,"IFM_61407"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
Afum_IFM_61407_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
Afum_IFM_61407_zoom
#slight bump at 1/2 

pdf(NULL)
dev.control(displaylist="enable")
hist(ad2_df[,"NCPF.7816"], breaks = seq(0,1,by=0.02), col = "#A9A578", xaxt="n", xlab =NA, ylab = NA)
hist(ad1_df[,"NCPF.7816"], breaks = seq(0,1,by=0.02), col = "#7D7A70", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
Afum_NCPF_7816_zoom <- recordPlot()
invisible(dev.off())
grid::grid.newpage()
Afum_NCPF_7816_zoom
#shoulders


#plot all 
dev.off()
p<-ggarrange(Afum_08_36_03_25_all, 
             Afum_IFM_61407_all,
             Afum_IFM_59359_all,
             Afum_AF100_12_5_all,
             Afum_DMC2_AF100_1_18_all,
             Afum_DMC2_AF100_1_3_all,
             Afum_Afu_B7586_CDC_30_all,
             Afum_Afu_343_P_11_all,
             Afum_NCPF_7816_all,
             ncol = 3, nrow = 3)
ggsave("ploidy_by_allele_diversity_all.pdf",p, width=28, height=32, units="in")

#plot zoom
q<-ggarrange(Afum_08_36_03_25_zoom, 
             Afum_IFM_61407_zoom,
             Afum_IFM_59359_zoom,
             Afum_AF100_12_5_zoom,
             Afum_DMC2_AF100_1_18_zoom,
             Afum_DMC2_AF100_1_3_zoom,
             Afum_B7586_CDC_30_zoom,
             Afum_Afu_343_P_11_zoom,
             Afum_NCPF_7816_zoom,
             ncol = 3, nrow = 3)
ggsave("ploidy_by_allele_diversity_zoom.pdf",q, width=18, height=22, units="in")

