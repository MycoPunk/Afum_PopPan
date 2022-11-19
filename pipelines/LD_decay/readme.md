This pipeline calculates Linkage Disequilibrium Decay. 

In the main directory, we have a set of .txt files (Clade_1_members.txt, Clade_2_members.txt, Clade_3_members.txt, and all_members.txt) that consist of lists of species names for each clade (clade membership determined previously using DAPCA) or all strains regardless of clade

```
>cat Clade_3_members.txt
CM2495
CM2730
CM2733
CM3249
CM3249b
CM3262
CM3720
...
```

We use these txt files to create 20 random subsets of 12 strains per clade from the inclusive vcf using step 00, then step 01 changes each vcf subset into plink format, step 02 runs LD decay, step 03 creates a summay file of pairwise distances, step 04 averages the summary files for position for individual runs, and step 05 averages the summaries for each position across all 20 reps.

**00_rarify_vcf.sh**<br>
This script randomly sorts the .txt files of strain names, and takes the first 12 lines from each random assortment (repeated 20 times for each clade). Then, using these lists, the sampled strains are subset from the inclusive VCF file, invariant sites removed, and each file named by group (clade 1-3 or all) and sample number (1-20). 

**01_format_as_plink.sh**<br>
This script recodes each subset vcf file into plink format. 

**02_run_LD_decay.sh**<br>
This script runs LD decay in plink

**03_create_summary.sh**<br>
This script creates a summary of pairwise distances w/ distance in the first col, and R2 in the second. <br>

**04_average_summaries_ea.sh** runs th script **04_average_summaries_ea.R** which generates the mean R2 value for each distance, generating 20 files per group (clades 1-3 or all). The BASH script itterates over each file produced in step 03 by creating a unique file name variable that is iteratively piped into the R script to significantly speed up the processing. Note, on smaller LD files, this step is not necesary, and everything can be averaged across all runs for all positions in one step (just skip to step 05 and remove the "ea" match from the input search (pattern = glob2rx("<name>*ea*")) to read in the input from step 03. However, at a sample size of n=12 X 20 reps on isolates with significant population structure, memory requirements are unreasonable and exceed 800GB for aggregation, so step 04 must be run first. Note that running or skipping step 04 will give you slightly different estimates for LD50 though, with increased differences depending on the amount of variance in the datset <br>

**05_average_summaries_all.sh** runs th script **05_average_summaries_all.R** which averages the mean R2 value for each distance for each itteration, generating 1 files per group (clades 1-3 or all).

**06_plot_LD_decay.sh** runs the script **06_plot_LD_decay.R** to calculate LD decay in BP and plot the figures presented in the paper. 
