**00_make_symlinks.sh**<br>The set up file to link the strains into the working folder.


**01_BLASTn_and_make_fastas.sh**<br>Runs BLASTn on the MAT idiomorphs and formats the significant hits as fasta files.<br>
  NOTE: for MAT1 alignments this will be pretty binary (hit or no hit), but for MAT2 alignments, there is a trailing ~270 bp region that is conserved between MAT1 and MAT2, and will always align, MAT2 idomorph strains will have a a full alignment over the (~1078bp) MAT2 reference, whereas MAT1 idomorph strains will only align over the ~270bp conserved region of the MAT2 reference.
  
  
**02_align_multi_fastas.sh**<br>This script runs a MAFT alignment of the above sequences, that you can view in your preferred viewer.
  
  
**03_further_eval_strains_with_hits_to_bothMAT.sh**<br>If all alignments are confidently MAT1 or MAT2, you don't need to run this. In our set, there are 9 strains with significant alignments to both MAT1 and MAT2. To make sure that this isn't some artifact in the scaffolds, we can look closer at the idiomorph regions by aligning the raw reads onto the MAT reference sequences. This script does that and pulls out new fasta files with any basepair changes included, adding "N"s for any position that fails to align.<br>
NOTE: symlink the raw reads (R1 and R2) for the strains of interest into the working directory first.
  
  
**04_fixnames_and_align_multifasta.sh**<br>This script will attach strain names to each of the output fasta files so you can tell them apart in the viewer, and then preform a multiple sequence alignment in MAFT as in step 2.
