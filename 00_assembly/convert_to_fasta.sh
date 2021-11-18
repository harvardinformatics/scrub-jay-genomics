#!/bin/bash

# code to convert hifiasm output to fasta

SAMPLE=$1

awk '/^S/{print ">"$2;print $3}' ${SAMPLE}.bp.hap1.p_ctg.gfa | gzip > ${SAMPLE}.hap1.fa.gz
awk '/^S/{print ">"$2;print $3}' ${SAMPLE}.bp.hap2.p_ctg.gfa | gzip > ${SAMPLE}.hap2.fa.gz
awk '/^S/{print ">"$2;print $3}' ${SAMPLE}.bp.pri.p_ctg.gfa | gzip > ${SAMPLE}.pri.fa.gz
