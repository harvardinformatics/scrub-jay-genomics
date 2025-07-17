grep PROJECTION loss_summ_data.tsv | awk '{if ($3 == "I" || $3 == "PI" || $3 == "UL") print $2}' | sort > I_PI_UL.txt

sort -k4,4 filt.query_annotation.bed -o filt.query_annotation.sorted.bed

join -1 1 -2 4 I_PI_UL.txt filt.query_annotation.sorted.bed -t '\t' > filt.query_annotation.I_PI_UL.bed