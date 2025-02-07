#!/usr/bin/env python3

'''
This script takes a gtf file where both geneIDs and transcriptIDs are the same;
It replaces gene_id filed with IDs provided in a separate gene table (transcID \t geneID)
'''

import argparse
from collections import defaultdict


def read_isoforms(isoforms):
	## Reads file transcriptID \t geneID into a dict

	iso_dict = {}
	with open(isoforms, 'r') as inf:
		for line in inf.readlines():
			transc = line.split()[0]
			gene = line.split()[1]
			iso_dict[transc] = gene
	return iso_dict


def assign_geneid_to_gtf(annogtf, iso_dict):
	## Reads a 9-field gtf file;
	## replaces gene_id field with geneIDs from the isoform dict

	with open(annogtf, 'r') as inf:
		for line in inf.readlines():
			parts_to_keep = '\t'.join(line.split('\t')[:-1])
			transc_info = line.split('\t')[-1].rstrip('\n')

			if transc_info.startswith('gene_id'):
				transc_info_to_keep = ';'.join(transc_info.split(';')[1:])
				gene_id_field = transc_info.split(';')[0]
				gene_id = gene_id_field.replace(' gene_id ', '').strip('"')
				transc_id_field = transc_info.split(';')[1]
				transc_id = transc_id_field.replace(' transcript_id ', '').strip('"')

				if transc_id in iso_dict:
					new_gene_id = iso_dict[transc_id]
				else:
					new_gene_id = gene_id	
				new_gtf_line = '{}\tgene_id "{}";{}'.format(parts_to_keep, new_gene_id, transc_info_to_keep)
				print(new_gtf_line)


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-a', '--annogtf', type=str, help='annotation file in gtf format')
	parser.add_argument('-i', '--isoforms', type=str, help='isoforms file: transcriptID \t geneID')
	args = parser.parse_args()

	## Parse arguments
	annogtf = args.annogtf
	isoforms = args.isoforms

	## Read isoforms
	iso_dict = read_isoforms(isoforms)

	## Output new gtf file with geneIDs
	assign_geneid_to_gtf(annogtf, iso_dict)


if __name__ == "__main__":
	main()