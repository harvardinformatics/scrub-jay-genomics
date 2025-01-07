#!/usr/bin/env python3
#
# This script takes a file with numeric values in some columns; calculates sum per row;
# outputs sums or outputs rows with sum above threshold if parameter -c is provided

import argparse
import math


__author__ = "Ekaterina Osipova, 2022."


def sum_one_row(line):
    ## Returns a sum of all numberic values in a row

    numbers = [float(i) for i in line.split() if i.isdigit()]
    return sum(numbers)


def main():
    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, help='input file with numeric columns')
    parser.add_argument('-c', '--cutoff', type=float, default=-math.inf, help='cutoff for the sum of a row to stay')
    args = parser.parse_args()
    
    ## Read input file, count sum per row -> output
    cutoff = args.cutoff
    with open(args.input, 'r') as inf:
        for line in inf.readlines():
            row_sum = sum_one_row(line)
            if cutoff == -math.inf:
                print(row_sum)
            else:
                if row_sum >= cutoff:
                    print(line.rstrip())


if __name__ == '__main__':
    main()
