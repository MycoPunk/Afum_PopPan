This pipeline calculates Linkage Disequilibrium Decay. 

In our main file, we have a set of .txt files (Clade_1_members.txt, #Clade_2_members.txt, Clade_3_members.txt) that consist of lists of species names for each clade (clad emembership determined previously using DAPCA)
>cat Clade_3_members.txt
CM2495
CM2730
CM2733
CM3249
CM3249b
CM3262
CM3720
CM4602
CM4946
CM7560
F15927
LMB-35Aa
MO79587EXP
TP32
W72310

#we will use these txt files to create random subsets of the inclusive vcf using the script:

*00_rarify_vcf.sh*

This script randomly sorts the .txt files of strain names, and takes the first 8 lines from each random assortment (repeated 50 times for each clade). Then, using these lists, the sampled strains are subset from the inclusive VCF file, invariant sites removed, and each file named by clade (1-3) and sample number (1-50). 
