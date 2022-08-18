#This script generates figures and runs statistical analysis to look for enrichment of core/accessory/singleton genes relative to genomic location (distance to telo. ends)

#load libraries
library(tidyverse)
library(fs)
library(dplyr)
library(cowplot)
library(scales)

#set wd
setwd("~/Desktop/Afumigatus_pangenome-master/Pangenome_analysis/Subtelo_analysis/")

#set pallet
singleton_col<- "#316A6E"
#accessory_col<- "#BA9141"
#w/ increased contrast
accessory_col<- "#E6B450"
#core_col<- "#806633"
#w/ increased contrast
core_col<- "#665229"
#check
show_col(accessory_col)
show_col("#BA9141")

data_dir = "reports"
obs_tsv_files <- fs::dir_ls(data_dir, regexp = "\\.observed_50k.tsv$")
obsdata <- obs_tsv_files %>% map_dfr(read_tsv,show_col_types = FALSE)

obsAdd <- obsdata %>% mutate(SOURCE = "SubTelomere")

rand_tsv_files <- fs::dir_ls(data_dir, regexp = "\\.random_50k.tsv$")
rand_data <- rand_tsv_files %>% map_dfr(read_tsv,show_col_types = FALSE)

randAdd <- rand_data %>% mutate(SOURCE="Random")

teloGeneType <- bind_rows(obsAdd,randAdd) %>% filter (TOTAL >0) %>% mutate(ratio = CORE / ACCESSORY,
                                                                           COREP = CORE / TOTAL,
                                                                           ACCESSORYP = ACCESSORY / TOTAL,
                                                                           SINGLEP = SINGLETON / TOTAL,
                                                                           DISPENSABLE=ACCESSORY+SINGLETON,
                                                                           DISPENSABLEP=(ACCESSORY+SINGLETON)/TOTAL)



###PLOTS 

###SINGLETONS###

##significant?
#check normalcy
teloGeneType %>% filter(SINGLETON!= 0) %>% filter(SOURCE== "SubTelomere") %>% shapiro_test(SINGLETON)
#not normal
teloGeneType %>% filter(SINGLETON!= 0) %>% filter(SOURCE== "Random") %>% shapiro_test(SINGLETON)
#not normal

#test variance
var.test(SINGLETON ~ SOURCE, data = teloGeneType %>% filter(SINGLETON!= 0), 
         alternative = "two.sided") #var not equal

#data is not normal, use wilcox rank sum test
wilcox_singleton<- wilcox.test(SINGLETON ~ SOURCE,
                               data   = teloGeneType %>% filter(SINGLETON!= 0), 
                               alternative = "two.sided", 
                               paired = FALSE)

wilcox_singleton$p.value
#P_val<- as.character(paste("",format(wilcox_singleton$p.value, scientific = TRUE, digits = 2)))


#plot
singleton_densi<- ggplot(teloGeneType %>% filter(SINGLETON!= 0), aes(x=SINGLETON, fill=SOURCE)) + 
    geom_density(alpha=.8) + 
  xlab("n singleton gene fams in 50kb windows")+
  scale_fill_manual(values=c("#4D4D4D", singleton_col))+ #need to change scale - very difficult to read. 
theme(text=element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
      legend.position = "right", 
      legend.title=element_blank())+
  labs(y = "Density")
  #add pvalue
temp<- teloGeneType %>% filter(SINGLETON!= 0)
x_pos<- median(unique(temp$SINGLETON))
singleton_densi<- singleton_densi+ annotate(geom="text", label = paste0("italic(p)", "==", format(wilcox_singleton$p.value, scientific = TRUE, digits = 2)), parse = TRUE,x=3,y=7.5,size=3)

singleton_densi



singleton_box <- ggplot(teloGeneType %>% filter(SINGLEP!= 0),aes(x=SOURCE,y=SINGLEP,fill=SOURCE)) + geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE, alpha=.8) +
  scale_fill_manual(values=c("#4D4D4D", singleton_col))+ 
  ylab("% singleton gene fams in 50kb windows") + xlab(NULL) +
  theme(text=element_text(size=12), panel.grid.minor = element_blank(),
        legend.position = "right", 
        legend.key = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black",
                                        size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                        colour = "light grey"),
        legend.title=element_blank())
  singleton_box




###Accessory###
  ##significant?
  #check normalcy
  teloGeneType %>% filter(ACCESSORY!= 0) %>% filter(SOURCE== "SubTelomere") %>% shapiro_test(ACCESSORY)
  #not normal
  teloGeneType %>% filter(ACCESSORY!= 0) %>% filter(SOURCE== "Random") %>% shapiro_test(ACCESSORY)
  #dataframe too large (see table) but SubTelomere not normal anyway, use wilcox
  test<- table(teloGeneType$ACCESSORY, teloGeneType$SOURCE)
  
  #test variance
  var.test(ACCESSORY ~ SOURCE, data = teloGeneType %>% filter(ACCESSORY!= 0), 
           alternative = "two.sided") #var not equal
  
  #data is not normal, use wilcox rank sum test
  wilcox_accessory<- wilcox.test(ACCESSORY ~ SOURCE,
                                 data   = teloGeneType %>% filter(ACCESSORY!= 0), 
                                 alternative = "two.sided", 
                                 paired = FALSE)
  wilcox_accessory
  wilcox_accessory$p.value

  
#generate plot  
accessory_densi<- ggplot(teloGeneType %>% filter(ACCESSORY!= 0), aes(x=ACCESSORY, fill=SOURCE))+ 
geom_density(alpha=.8) + 
  xlab("n accessory gene fams in 50kb windows")+
  scale_fill_manual(values=c("#4D4D4D", accessory_col))+ #need to change scale - very difficult to read. 
  theme(text=element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        legend.position = "right",
        legend.title=element_blank()) +
  labs(y = "Density")
temp<- teloGeneType %>% filter(ACCESSORY!= 0)
x_pos<- median(unique(temp$ACCESSORY))
accessory_densi<- accessory_densi+ annotate(geom="text", label = paste0("italic(p)", "==", signif(wilcox_accessory$p.value, digits = 3)), parse = TRUE,x=10,y=0.9,size=3)

accessory_densi 


#box plots
accessory_box <- ggplot(teloGeneType %>% filter(ACCESSORYP!= 0),aes(x=SOURCE,y=ACCESSORYP,fill=SOURCE)) + geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE, alpha=.8) + 
ylab("% accessory gene fams in 50kb windows") + xlab(NULL)+   scale_fill_manual(values=c("#4D4D4D", accessory_col))+ 
  theme(text=element_text(size=12), panel.grid.minor = element_blank(),
        legend.position = "right", 
        legend.key = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black",
                                        size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                        colour = "light grey"),
        legend.title=element_blank())

accessory_box

#Accessory genes enriched in subtelomeres 



###CORE###
##significant?
#check normalcy
teloGeneType %>% filter(CORE!= 0) %>% filter(SOURCE== "SubTelomere") %>% shapiro_test(CORE)
#not normal
teloGeneType %>% filter(CORE!= 0) %>% filter(SOURCE== "Random") %>% shapiro_test(CORE)
#dataframe too large (see table) but SubTelomere not normal anyway, use wilcox
test<- table(teloGeneType$CORE, teloGeneType$SOURCE)

#test variance
var.test(CORE ~ SOURCE, data = teloGeneType %>% filter(CORE!= 0), 
         alternative = "two.sided") #var not equal

#data is not normal, use wilcox rank sum test
wilcox_core<- wilcox.test(CORE ~ SOURCE,
                          data   = teloGeneType %>% filter(CORE!= 0), 
                          alternative = "two.sided", 
                          paired = FALSE)
#difference is significant
wilcox_core
wilcox_core$p.value


#render plot
core_densi <- ggplot(teloGeneType %>% filter(CORE!= 0), aes(x=COREP, fill=SOURCE)) + 
  geom_density(alpha=.8) + 
  xlab("n core gene fams in 50kb windows")+
  scale_fill_manual(values=c("#4D4D4D", core_col))+ #need to change scale - very difficult to read. 
  theme(text=element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        legend.position = "right",
        legend.title=element_blank())+
  labs(y = "Density")
temp<- teloGeneType %>% filter(CORE!= 0)
#x_pos<- median(unique(temp$CORE))
core_densi<- core_densi+ annotate(geom="text", label = paste0("italic(p)", "==", signif(wilcox_core$p.value, digits = 3)), parse = TRUE,x=0.4,y=10,size=3)
core_densi


core_box <- ggplot(teloGeneType %>% filter(COREP!= 0),aes(x=SOURCE,y=COREP,fill=SOURCE)) + geom_boxplot(outlier.colour="black", outlier.shape=16,outlier.size=2, notch=FALSE, alpha=.8)+
ylab("% core gene fams in 50kb windows") + xlab(NULL)+   scale_fill_manual(values=c("#4D4D4D", core_col))+ 
  theme(text=element_text(size=12), panel.grid.minor = element_blank(),
        legend.position = "right", 
        legend.key = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black",
                                        size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                        colour = "light grey"),
        legend.title=element_blank())
core_box
#core genes are depleted in the subtelomeres



###make compound figure###
p<- cowplot::plot_grid(core_densi, core_box,
                       accessory_densi, accessory_box,
                       singleton_densi, singleton_box,
                       align = 'h',
                       nrow = 3, 
                       ncol = 2,
                       rel_widths = c(.8,.8,.8), 
                       labels = c('A', 'B', 'C', 'D', 'E', 'F'))
p
ggsave(file="../sub_teleo_plot_50.pdf",device="pdf", p, width=10, height=10, units="in")



#get total counts in sub telo. vs overall. 
core_random<- data.frame(teloGeneType %>% filter(CORE!= 0) %>% filter(SOURCE== "Random"))
n_CORE_random<- sum(core_random$CORE)
n_CORE_random #2,095,777

core_SubTelomere<- data.frame(teloGeneType %>% filter(CORE!= 0) %>% filter(SOURCE== "SubTelomere"))
n_CORE_SubTelomere<- sum(core_SubTelomere$CORE)
n_CORE_SubTelomere #6,472
#out of 2,102,249 core genes detected on 50kb + scaffolds 6,472 were subtelomeric, which is significantly depleted over random. 
#assuming there are 8866 core genes what kind of ratio do we get for 50kb scaffolds?
2102249/8866 #237.114


#accessory
accessory_random<- data.frame(teloGeneType %>% filter(ACCESSORY!= 0) %>% filter(SOURCE== "Random"))
n_ACCESSORY_random<- sum(accessory_random$ACCESSORY)
n_ACCESSORY_random #260,332

accessory_SubTelomere<- data.frame(teloGeneType %>% filter(ACCESSORY!= 0) %>% filter(SOURCE== "SubTelomere"))
n_ACCESSORY_SubTelomere<- sum(accessory_SubTelomere$ACCESSORY)
n_ACCESSORY_SubTelomere #1,867
#out of 262,199 accessory genes detected on 50kb + scaffolds 1,867 were subtelomeric, which is significantly enriched over random. 

#assuming there are 4334 accessory genes what kind of ratio do we get for 50kb scaffolds?
262199/4334 #60.498


#singleton
singleton_random<- data.frame(teloGeneType %>% filter(SINGLETON!= 0) %>% filter(SOURCE== "Random"))
n_SINGLETON_random<- sum(singleton_random$SINGLETON)
n_SINGLETON_random #1,277

singleton_SubTelomere<- data.frame(teloGeneType %>% filter(SINGLETON!= 0) %>% filter(SOURCE== "SubTelomere"))
n_SINGLETON_SubTelomere<- sum(singleton_SubTelomere$SINGLETON)
n_SINGLETON_SubTelomere #14
#out of 1241 singletons detected on 50kb + scaffolds 14 were subtelomeric, which is significantly enriched over random. 

#assuming there are 2109 singleton genes what kind of ratio do we get for 50kb scaffolds?
(1241/2109) *260 #0.588
