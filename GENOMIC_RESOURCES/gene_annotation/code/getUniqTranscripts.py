#!/usr/bin/env python3


# To get unique transcripts of a file run e.g:\
# e.g usage: getUniqTrascripts.py -f HLparMaj1.ncbi.bed > uniq.HLparMaj1.ncbi.bed


import argparse
import sys

__author__ = "Ekaterina Osipova, 2019."


def get_uniq_transcripts(bed_file, a_score, a_color):
    ## Makes a dictionary of unique transcripts considering which fields are important

    overlap_transcripts = {}
    with open(bed_file, 'r') as inf:
        for line in inf.readlines():
            transc_name, score, color = line.split('\t')[3], line.split('\t')[4], line.split('\t')[8]

            # get transcript info considering which fields are important
            if a_color and a_score:
                # both score and color are considered
                transc_info = '\t'.join(line.split('\t')[:3] + line.split('\t')[4:])
            elif a_color and not a_score:
                # only color is considered
                transc_info = '\t'.join(line.split('\t')[:3] + line.split('\t')[5:])
            elif a_score and not a_color:
                # only score is considered
                transc_info = '\t'.join(line.split('\t')[:3] + line.split('\t')[4:8] + line.split('\t')[9:])
            else:
                # both color and score are ignored
                transc_info = '\t'.join(line.split('\t')[:3] + line.split('\t')[5:8] + line.split('\t')[9:])


            if (transc_info in overlap_transcripts):
                sys.stderr.write('DUPLICATION: {}\t{}'.format(transc_name, '\t'.join(transc_info.split('\t'))+'\n'))
            else:
                overlap_transcripts[transc_info] = (transc_name, score, color)
    return  overlap_transcripts


def output_uniq_transcripts(overlap_transcripts, a_color, a_score):
    ## Outputs unique elements of the overlap_transcripts dictionary

    for transc_info in overlap_transcripts:
        transc_var_parts = overlap_transcripts[transc_info]

        # assemble a new annotation line considering fields excluded before
        if a_color and a_score:
            # both score and color are considered
            newBedLine = '\t'.join(transc_info.split('\t')[:3]) + '\t' + transc_var_parts[0] + '\t' + \
                         transc_var_parts[1] + '\t' + '\t'.join(transc_info.split('\t')[4:7]) + '\t' + \
                         transc_var_parts[2] + '\t' + '\t'.join(transc_info.split('\t')[8:])
        elif a_color and not a_score:
            # only color is considered
            newBedLine = '\t'.join(transc_info.split('\t')[:3]) + '\t' + transc_var_parts[0] + '\t' + \
                         transc_var_parts[1] + '\t' + '\t'.join(transc_info.split('\t')[3:6]) + '\t' + \
                         transc_var_parts[2] + '\t' + '\t'.join(transc_info.split('\t')[7:])
        elif a_score and not a_color:
            # only score is considered
            newBedLine = '\t'.join(transc_info.split('\t')[:3]) + '\t' + transc_var_parts[0] + '\t' + \
                         transc_var_parts[1] + '\t' + '\t'.join(transc_info.split('\t')[4:7]) + '\t' + \
                         transc_var_parts[2] + '\t' + '\t'.join(transc_info.split('\t')[7:])
        else:
            # both color and score are ignored
            newBedLine = '\t'.join(transc_info.split('\t')[:3]) + '\t' + transc_var_parts[0] + '\t' + \
                         transc_var_parts[1] + '\t' + '\t'.join(transc_info.split('\t')[3:6]) + '\t' + \
                         transc_var_parts[2] + '\t' + '\t'.join(transc_info.split('\t')[6:])
        print(newBedLine.rstrip('\n'))


def main():
    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--filebed', type=str, help='bed12 file')
    parser.add_argument('-c', '--color', action='store_true',
                        help='if specified, RGB column must also be identical in duplicated transcripts')
    parser.add_argument('-s', '--score', action='store_true',
                        help='if specified, score column must also be identical in duplicated transcripts')
    args = parser.parse_args()

    ## Make a dictionary of unique transcripts
    overlap_transcripts = get_uniq_transcripts(args.filebed, args.score, args.color)
    
    ## Print unique transcripts
    output_uniq_transcripts(overlap_transcripts, args.score, args.color)


if __name__ == "__main__":
    main()
