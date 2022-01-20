library(dplyr)
library(stringr)
library(ggplot2)
library(purrr)
library(data.table)

setwd("~/bigdata/pop_genomics/LD_decay_4/LD_out/")

#Note- first we define a variable in the BASH submission script called name, representing the file names we want to process
#for N in LD_out/test/*; do \
#name=$(basename "$N")
#Rscript pipeline/05_average_summaries_ea.R $name --save
#done

#then, run the R script with just these two calls to look at the arguments list
args <- commandArgs()
print(args)
#look at the log file for the .sh run, it will tell you which argument number the variable name is (this will be 6, unless you're making doubble calls or have defined something else already)

#set argument number to pull in variable (file to loop over) form the bash script
name <- args[6]
name

#read in the file
input_df<- fread(name)

#set column names
colnames(input_df) <- c("dist","rsq")

#get mean distances R2 values for each distance across the replicates
input_df_rsq_means<- aggregate(rsq~dist, data=input_df, FUN=function(x) c(mean=mean(x), count=length(x)))

#clean up memory
rm(input_df)

#define output file
output_file_name <- paste(name,"_ea.tab",sep="")

#write output file
write.table(input_df_rsq_means, output_file_name, sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
