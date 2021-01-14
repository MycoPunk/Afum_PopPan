This pipeline calculates Linkage Disequilibrium Decay. 

In our main file, we have a set of .txt files (Clade_1_members.txt, #Clade_2_members.txt, Clade_3_members.txt) that consist of lists of species names for each clade (clad emembership determined previously using DAPCA)

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

We use these txt files to create 50 random subsets of 8 strains per clade from the inclusive vcf using step 00, then step 01 changes each vcf subset into plink format, step 02 runs LD decay, step 03 creates a summay file of pairdise distances. 

**00_rarify_vcf.sh**
This script randomly sorts the .txt files of strain names, and takes the first 8 lines from each random assortment (repeated 50 times for each clade). Then, using these lists, the sampled strains are subset from the inclusive VCF file, invariant sites removed, and each file named by clade (1-3) and sample number (1-50). 

**01_format_as_plink.sh** 
This script recodes each subset vcf file into plink format. 

**02_run_LD_decay.sh** 
This script runs LD decay in plink

**03_create_summary.sh**
This script creates a summary of pairwise distances w/ distance in the first col, and R2 in the second. 
