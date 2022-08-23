This pipeline runs Treemix and optM to look at historical introgression events and the direction of gene flow 

This pipeline requires a VCF file with an appropriate outgroup, and a population file in the form: 

```
>cat strains.clust3
Afischeri	Afischeri	outgroup
08-19-02-30	08-19-02-30	Clade_1.4
09-7500806	09-7500806	Clade_1.4
117535A-1.1	117535A-1.1	Clade_1.2
...
```

**01_remove_missing_sites.sh**<br>
This script removes missing sites from the main VCF (TreeMix does not like missing data).<br>

**02_ldPruning.sh**<br>
This script removes high LD sites at a 0.2 threshold value.<br>

**03_downsample.vcf**<br>
This script downsamples the VCF to 80% of the SNP sites replicated ten times - these VCFs are ultimately used to calculate the optimum m value in OptM, which will error if any SD values between subsamples are 0 (which can happen, particularly if your data is really clean). Taking random draws ensures there will be some variation.<br>

**04_vcf_to_TreeMix_downsamples.sh**<br>
This script converts the vcf to the treemix format required by treemix. It uses 04_plink2treemix.py. Note, this is the loop version to run on the downsampled reps.<br>

**04_plink2treemix.py**<br>
This is a.py (v2) script modified from one written by Joana Meier for format conversion.<br>

**05_runTreeMix_downsample.sh**<br> 
This script runs TreeMix with no migration edges as an input for OptM. Note, setting k has little effect here as we've already pruned out high LD sites, so we just exclude the blocking call.<br>

**06_run_optM.R**<br> 
This script runs optM and produces graphs used to determine the optimum number of migration edges.<br>

**07_vcf_to_TreeMix_optimum.sh**<br>
This script now runs the conversion for the full corrected VCF.<br>

**08_runTreeMix_on_optimum_m.sh**<br>
This script runs TreeMix on the full corrected VCF - NOTE you need to update the script to input your the m value predicted by OptM.<br>

**09_graphTreeMix.R**<br> 
This script graphs the treemix tree. Tequires that you have the Treemix plotting functions installed. (plotting_funcs.R). 
