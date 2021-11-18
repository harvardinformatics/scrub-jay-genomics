#!/bin/bash

#SBATCH -N 1
#SBATCH -t 1-00:00:00
#SBATCH --mem 400G
#SBATCH -n 24
#SBATCH --partition holy-cow
#SBATCH -J hifiasm

./run_hifiasm.sh $1
