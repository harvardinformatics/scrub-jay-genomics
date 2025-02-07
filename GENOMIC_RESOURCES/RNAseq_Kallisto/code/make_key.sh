
for abun in ../results/*/abundance.tsv
do
   # echo $abun
    sampID=$(echo $abun | awk 'BEGIN{FS="/"} {split($3,a,"_"); print a[4]}')
   # echo $sampID
    rnaID=$(echo $abun | awk 'BEGIN{FS="/"} {print $3}' | sed 's/_kallisto//g')
    tissue=$(echo $abun | awk 'BEGIN{FS="/"} {split($3,a,"_"); print a[5]}')
   # echo $tissue
    fullID=$(grep $sampID ../../assemblies/sample_sheet.txt | cut -f1)
  #  echo $fullID
    spp=$(echo $fullID | awk 'BEGIN{FS="_"} {print $1}')
    #echo "$rnaID $fullID $tissue"
    printf "$rnaID\t$fullID.$tissue\n"

done
