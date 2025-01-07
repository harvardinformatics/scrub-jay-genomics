#!/bin/bash
#SBATCH -n 2                # Number of cores
#SBATCH -N 1
#SBATCH -t 3-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared,holy-smokes   # Partition to submit to
#SBATCH --mem=32G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o collate.out # File to which STDOUT will be written
#SBATCH -e collate.err # File to which STDERR will be written

for abun in ../results/*_redo/abundance.tsv
do
   # echo $abun
    sampID=$(echo $abun | awk 'BEGIN{FS="/"} {split($3,a,"_"); print a[4]}')
   # echo $sampID
    tissue=$(echo $abun | awk 'BEGIN{FS="/"} {split($3,a,"_"); print a[5]}')
   # echo $tissue
    fullID=$(grep $sampID ../../assemblies/sample_sheet.txt | cut -f1)
  #  echo $fullID
    spp=$(echo $fullID | awk 'BEGIN{FS="_"} {print $1}')
   # echo "$sampID $fullID $tissue $spp"

    count=1
    while read line
    do
        transID=$(echo $line | cut -d ' ' -f1)
        #geneID=$(grep -w $transID ../../gene_annotation/results/merged_TOGA/gene_ids.tsv | cut -f2)
	geneID=$(grep -w $transID ../../gene_annotation/results/fixed_merged_TOGA/merged_TOGA_tgutRef_noDup.tsv | cut -f8)
        if [[ $count -ge 2 ]]
        then
            printf "$line\t$fullID\t$tissue\t$spp\t$geneID\n"
        fi
        ((count+=1))
    done < $abun

done
