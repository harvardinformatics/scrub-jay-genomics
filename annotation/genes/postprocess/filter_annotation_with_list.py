#!/usr/bin/env python3

import sys
import argparse
from collections import defaultdict


__author__ = "Ekaterina Osipova, 2022."


def main():
    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--list', type=str, help='list of gene/transcript IDs to filter for')
    parser.add_argument('-a', '--annotation', type=str, help='annotation of transcripts: define the column with IDs!')
    parser.add_argument('-c', '--column', type=int, default=4, help='column with transcript IDs; default: 4 (bed12)')
    parser.add_argument('-s', '--suffix', action='store_true', help='specify if anno IDs have .version youd like to remove')
    parser.add_argument('-b', '--but', action='store_true', help='ALL-BUT-LIST: specify if you want to filter OUT IDs from the list')
    args = parser.parse_args()

    ## Read IDs into a list
    id_list = [] 
    with open(args.list, 'r') as inf:
        for line in inf.readlines():
           id_list.append(line.rstrip())
    
    ## Read annotation, output only IDs present in the list 
    column = args.column
    with open(args.annotation, 'r') as inf:
        for line in inf.readlines():
            anno_id = line.split()[column - 1]
            if args.suffix:
                anno_id = anno_id.split('.')[0]
            if args.but:
                if anno_id not in id_list:
                    print(line.rstrip())
            else:
                if anno_id in id_list:
                    print(line.rstrip())



if __name__ == "__main__":
    main()
