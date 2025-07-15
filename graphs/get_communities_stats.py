#!/usr/bin/env python3
import sys

#Give FAI for community as argument 1
comm = {}

with open(sys.argv[1],"r") as f:
    for line in f:
        entry = line.split("\t")
        length = entry[1]

        #get sample#haplo
        contig = entry[0].split("#")
        haplo = contig[0] + "#" + contig[1]
#        print (haplo, "\t", length)

        if haplo in comm:
            val = comm[haplo]
            comm[haplo] += int(length)
        else:
            comm[haplo] = int(length)
    f.close()

for key in comm:
    val = comm[key]
    print (sys.argv[1], "\t", key, "\t", val)
