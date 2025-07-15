#############
## Panacus ##
#############
cd /n/holyscratch01/edwards_lab/bfang/Scrub_jay/Panacus_SJ
module load python
source activate fasrc2

### merge communities (chromosomes)
## prep GFA files directories in a TSV file to merge
# for a in `cat Tag_communities.txt`; do
# ls /n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/results/PGGB/allbird_${a}/*final_nameFix.gfa >> GFA_31communities.txt
# done

odgi squeeze -f GFA_29communities.txt -o - -t 10 | odgi view -i - -g > GFA_29communities.gfa & odgi squeeze -f GFA_31communities.txt -o - -t 10 | odgi view -i - -g > GFA_31communities.gfa

### 31 communities
grep '^P' GFA_31communities.gfa | cut -f2 | grep 'AW' > 31communities_community_AW.txt
grep '^P' GFA_31communities.gfa | cut -f2 | grep 'AC' > 31communities_community_AC.txt
grep '^P' GFA_31communities.gfa | cut -f2 | grep 'AI' > 31communities_community_AI.txt
grep '^P' GFA_31communities.gfa | cut -f2 | grep -E 'AW|AC|AI' | sort > 31communities_community_All.txt

RUST_LOG=info panacus histgrowth -t30 -l 1,1,1 -q 0,0.1,0.90 -S -c bp -a -s 31communities_community_AW.txt GFA_31communities.gfa > 31communities_AW_histgrowth_node.tsv
RUST_LOG=info panacus histgrowth -t30 -l 1,1,1 -q 0,0.1,0.90 -S -c bp -a -s 31communities_community_AC.txt GFA_31communities.gfa > 31communities_AC_histgrowth_node.tsv
RUST_LOG=info panacus histgrowth -t30 -l 1,1,1 -q 0,0.1,0.90 -S -c bp -a -s 31communities_community_AI.txt GFA_31communities.gfa > 31communities_AI_histgrowth_node.tsv
RUST_LOG=info panacus histgrowth -t30 -l 1,1,1 -q 0,0.1,0.90 -S -c bp -a -s 31communities_community_All.txt GFA_31communities.gfa > 31communities_All_histgrowth_node.tsv

panacus-visualize -e 31communities_AW_histgrowth_node.tsv > PDF_out/31communities_AW_histgrowth_node.pdf
panacus-visualize -e 31communities_AC_histgrowth_node.tsv > PDF_out/31communities_AC_histgrowth_node.pdf
panacus-visualize -e 31communities_AI_histgrowth_node.tsv > PDF_out/31communities_AI_histgrowth_node.pdf
panacus-visualize -e 31communities_All_histgrowth_node.tsv > PDF_out/31communities_All_histgrowth_node.pdf
