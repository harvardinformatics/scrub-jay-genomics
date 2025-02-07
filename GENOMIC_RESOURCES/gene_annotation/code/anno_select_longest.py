#!/usr/bin/env python3

import sys
import argparse
from collections import defaultdict
from operator import itemgetter


__author__ = "Ekaterina Osipova, 2022."


def read_into_dict(file):
    ## Read isoformes into dict

    iso_dict = defaultdict(list)
    with open(file, 'r') as inf:
        for line in inf.readlines():
           gene = line.rstrip().split()[0]
           trans = line.rstrip().split()[1]
           iso_dict[gene].append(trans)
    return iso_dict


def read_anno_into_dict(anno, field):
    ## Reads bed12 annotation file into a dictionary ID : (LEN, trancs_info)

    anno_dict = {}
    with open(anno, 'r') as inf:
        for line in inf.readlines():
            if len(line.split()) < field - 1:
                print('There is no field {} in the annotation! Abort'.format(field))
                sys.exit(1)
            else:
                id = line.split()[field - 1]
                exon_lens = line.split()[10].rstrip(',').split(',')
                transc_len = sum([int(i) for i in exon_lens])
                transc_info = line.rstrip()
                anno_dict[id] = (transc_len, transc_info)
    return anno_dict


def main():
    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--iso', type=str, help='isoformes file : gene \t transcript')
    parser.add_argument('-a', '--anno', type=str, help='annotationin bed12 format') 
    args = parser.parse_args()

    ## Read isoforms into a dict
    iso_dict = read_into_dict(args.iso)
    
    ## Read annotation into a dict
    field = 4
    anno_dict = read_anno_into_dict(args.anno, field)
    
    ## Output the longest isoform
    for gene in iso_dict:
        try:
            iso_info_list = [anno_dict[iso] for iso in iso_dict[gene]]
            longest_iso = max(iso_info_list, key=itemgetter(0))[1]
            print(longest_iso)
        except KeyError:
            print(gene + " not found")

if __name__ == "__main__":
    main()
