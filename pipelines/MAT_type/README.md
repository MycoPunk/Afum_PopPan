**00_make_symlinks.sh** = the set up file to link the strains into the working folder<br>
**01_BLASTn_and_make_fastas.sh** = runs BLASTn on the MAT idiomorphs and formats the significant hits as fasta files.<br>
  NOTE: for MAT1 alignments this will be pretty binary (hit or no hit), but for MAT2 alignments, there is a trailing ~270 bp region that is conserved between MAT1 and MAT2, and will always align, MAT2 idomorph strains will have a a full alignment over the (~1078bp) MAT2 reference, whereas MAT1 idomorph strains will only align over the ~270bp conserved region of the MAT2 reference.<br>
  **02_align_multi_fastas.sh** = this script runs a MAFT alignment of the above sequences, that you can view in your prefered viewer. 
