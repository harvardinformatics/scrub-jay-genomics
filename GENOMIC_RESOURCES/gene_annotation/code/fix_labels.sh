while read line
do
    sampID=$(echo "$line" | awk 'BEGIN{FS="\t"} {print $1}')
   # echo $sampID
    replace=$(echo "$line" | awk 'BEGIN{FS="\t"} {print $2}')
#    echo "$sampID     $replace"
    sed -i "s/$sampID/$replace/g" temp.r

done < RNAseq_sampleID_key.tsv
