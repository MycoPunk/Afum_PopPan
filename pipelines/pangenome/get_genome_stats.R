
#get genome statistics 

#run other genome stats quick
genome_stats<-read.delim("data/genome_stats.csv", header = TRUE, sep = ",", fill = TRUE, strip.white = TRUE)
names(genome_stats)

#busco stats
mean(genome_stats$BUSCO.score)
max(genome_stats$BUSCO.score)
min(genome_stats$BUSCO.score)

#depth stats
mean(genome_stats$sequence.depth..bp.)
max(genome_stats$sequence.depth..bp.)
min(genome_stats$sequence.depth..bp.)

#L50 stats
mean(genome_stats$L50)
max(genome_stats$L50)
min(genome_stats$L50)

#N50 stats
mean(genome_stats$N50)
max(genome_stats$N50)
min(genome_stats$N50)

#assembily size stats
mean(genome_stats$Assembly.size..bp.)
max(genome_stats$Assembly.size..bp.)
min(genome_stats$Assembly.size..bp.)

#is there a relationship between genome size and the number of singleton / accessory genes?
n_singletons_by_strain<-read.delim("data/n_singletons_by_strain_OF.txt", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE)
n_accessory_by_strain<-read.delim("data/n_accessory_by_strain_OF.txt", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE)

#fix names so that they match
genome_stats<- data.frame(sapply(genome_stats, gsub, pattern = "-", replacement = "_"))
n_singletons_by_strain<- data.frame(sapply(n_singletons_by_strain, gsub, pattern = "-", replacement = "_"))
n_accessory_by_strain<- data.frame(sapply(n_accessory_by_strain, gsub, pattern = "-", replacement = "_"))
genome_stats<- data.frame(sapply(genome_stats, gsub, pattern = "A_fum_", replacement = ""))

#get stragglers 
genome_stats<- data.frame(sapply(genome_stats, gsub, pattern = "DMC2_", replacement = ""))
n_singletons_by_strain<- data.frame(sapply(n_singletons_by_strain, gsub, pattern = "DMC2_", replacement = ""))
n_accessory_by_strain<- data.frame(sapply(n_accessory_by_strain, gsub, pattern = "DMC2_", replacement = ""))

genome_stats<- data.frame(sapply(genome_stats, gsub, pattern = "AFIS1435CDC_6", replacement = "AFIS1435_CDC_6"))
genome_stats<- data.frame(sapply(genome_stats, gsub, pattern = "EXP", replacement = "EP"))
genome_stats<- data.frame(sapply(genome_stats, gsub, pattern = "Afu_343_P/11", replacement = "Afu_343_P_11"))
genome_stats<- data.frame(sapply(genome_stats, gsub, pattern = "117535A_1.1", replacement = "117535A_1_1"))
genome_stats<- data.frame(sapply(genome_stats, gsub, pattern = "B7594CDC_27", replacement = "B7594_CDC_27"))
genome_stats<- data.frame(sapply(genome_stats, gsub, pattern = "B5258CDC_22", replacement = "B5258_CDC_22"))
genome_stats<- data.frame(sapply(genome_stats, gsub, pattern = "B6894CDC_5", replacement = "B6894_CDC_5"))
genome_stats<- data.frame(sapply(genome_stats, gsub, pattern = "B5960CDC_29", replacement = "B5960_CDC_29"))

#check that names all match
setdiff(n_singletons_by_strain$strain, genome_stats$strain.name)

colnames(n_singletons_by_strain)<- c("n_singletons", "strain.name")
colnames(n_accessory_by_strain)<- c("n_accessory", "strain.name")

#combine count data by strain name 
genome_stats<- merge(genome_stats, n_singletons_by_strain[, c("strain.name", "n_singletons")], by="strain.name")
genome_stats<- merge(genome_stats, n_accessory_by_strain[, c("strain.name", "n_accessory")], by="strain.name")

#stats
lm_singletons<-lm(as.numeric(n_singletons) ~ as.numeric(Assembly.size..bp.), data = genome_stats)
summary(lm_singletons) #significant

lm_accessory<-lm(as.numeric(n_accessory) ~ as.numeric(Assembly.size..bp.), data = genome_stats)
summary(lm_accessory) #significant

#format long
genome_stats_subset<- as.data.frame(cbind(strain.name = genome_stats$strain.name, genome_size = genome_stats$Assembly.size..bp., n_singletons = genome_stats$n_singletons, n_accessory = genome_stats$n_accessory))

genome_stats_subset_wide<- genome_stats_subset %>%
  pivot_longer(cols = c(n_singletons, n_accessory), names_to = "group", values_to = "count")


#singleton = "#316A6E",
#accessory = "#BA9141"
# plot
plot<- ggplot(genome_stats_subset_wide, aes(x=as.numeric(genome_size), y=as.numeric(count),color=group)) + 
  geom_point(size=2, alpha = 0.5) +
  theme_minimal()
to_print<- plot + facet_wrap(~group, scales = "free") +
  scale_color_manual(values = c("#BA9141", "#316A6E"))+
  scale_fill_manual(values = c("#BA9141", "#316A6E"))+
  scale_color_manual(values=c("#BA9141", "#316A6E")) + geom_smooth(method=lm, fullrange=TRUE, se=FALSE)
  
to_print

ggsave("genome_size.pdf",to_print, width=8, height=4, units="in")


