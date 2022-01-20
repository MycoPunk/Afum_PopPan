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

**04_average_summaries.sh** runs th script **04_average_summaries.R** which generates the mean R2 value for each distance for each itteration, generating 20 files per group (clades 1-3 or all)

**06_plot_LD.sh** runs two plotting scripts: <br>
**05_plot_LD_decay_zoomed_in.R** plots a zoomed in view of LD decay, over 10bp windows, and <br>
**05_plot_LD_decay_LD50.R** plots LD decay over the entier length of the dataset, including arrows at LD50 for each clade
