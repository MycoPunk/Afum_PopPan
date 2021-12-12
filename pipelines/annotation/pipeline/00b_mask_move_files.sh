#!/bin/bash
for d in *.mask.*; do 
	b=$(echo $d | perl -p -e 's/\.mask\.\d+//'); 
	mv $d/funannotate-mask.log logs/mask.$b.log; 
	mv $d/repeatmodeler-library.*.fasta repeat_library/$b.repeatmodeler-library.fasta; 
	rmdir $d
done

