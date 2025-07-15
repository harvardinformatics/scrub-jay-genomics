[ -e results/assembly_qc/ASSEMBLY_STATS.tsv ] && rm results/assembly_qc/ASSEMBLY_STATS.tsv
echo -e "Sample\tpri.N50\tpri.Num_contigs\tpri.Tot_len\thap1.N50\thap1.Num_contigs\thap1.Tot_len\thap2.N50\thap2.Num_contigs\thap2.Tot_len" >> results/assembly_qc/ASSEMBLY_STATS.tsv

while read sample; do

    pN50=$(sed 's/ \+ /\t/g' results/assembly_qc/${sample}-p_ctg-quast/report.txt | grep 'N50' | cut -f2)
    pnumCon=$(sed 's/ \+ /\t/g' results/assembly_qc/${sample}-p_ctg-quast/report.txt | grep '# contigs'$'\t' | cut -f2)
    ptotLen=$(sed 's/ \+ /\t/g' results/assembly_qc/${sample}-p_ctg-quast/report.txt | grep 'Total length'$'\t' | cut -f2)


    h1N50=$(sed 's/ \+ /\t/g' results/assembly_qc/${sample}-hap1.p_ctg-quast/report.txt | grep 'N50' | cut -f2)
    h1numCon=$(sed 's/ \+ /\t/g' results/assembly_qc/${sample}-hap1.p_ctg-quast/report.txt | grep '# contigs'$'\t' | cut -f2)
    h1totLen=$(sed 's/ \+ /\t/g' results/assembly_qc/${sample}-hap1.p_ctg-quast/report.txt | grep 'Total length'$'\t' | cut -f2)


    h2N50=$(sed 's/ \+ /\t/g' results/assembly_qc/${sample}-hap2.p_ctg-quast/report.txt | grep 'N50' | cut -f2)
    h2numCon=$(sed 's/ \+ /\t/g' results/assembly_qc/${sample}-hap2.p_ctg-quast/report.txt | grep '# contigs'$'\t' | cut -f2)
    h2totLen=$(sed 's/ \+ /\t/g' results/assembly_qc/${sample}-hap2.p_ctg-quast/report.txt | grep 'Total length'$'\t' | cut -f2)
    
#    echo -e "${sample}${p2N50}\t${pnumCon}\t$ptotLen\t${h1N50}\t${h1numCon}\t$h1totLen\t${h2N50}\t${h2numCon}\t$h2totLen" >> results/assembly_qc/ASSEMBLY_STATS.tsv
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$sample" "$pN50" "$pnumCon" "$ptotLen" "$h1N50" "$h1numCon" "$h1totLen" "$h2N50" "$h2numCon" "$h2totLen" >> results/assembly_qc/ASSEMBLY_STATS.tsv

done < sample_IDs.txt
