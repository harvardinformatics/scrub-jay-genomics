#!/bin/bash

## code to run assembly of scrub jay genome ##

#sample id is first, threads is second
SAMPLE=$1
THREADS=24 #fixed to ntasks from slurm script
SEQFILES="/n/holylfs05/LABS/informatics/Everyone/scrubjay/RAW_DATA/$SAMPLE/*/*.fastq.gz"

../../hifiasm/hifiasm -t $THREADS -o $SAMPLE $SEQFILES
