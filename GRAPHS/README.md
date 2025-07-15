To generate a PGGB graph, we first partioned the combined fasta files of all assembled individuals (pggb_partition_v2.slurm) into communities. 
Then, produced a graph for each community (pggb_communities_array_v2.slurm)
The final graph was produced by combining each community graph.

Additional statistics and analysis of individual communities was done using `get_communities*` scripts. 
Comparison of core and accessory sequence was done using Panacus, see: `panacus.sh`

The `depth` subdirectory contains an analysis of graph depth.