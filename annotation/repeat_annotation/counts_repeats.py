import argparse
import json

argParser = argparse.ArgumentParser()
argParser.add_argument("-i", "--input", dest='fofn', help="FoFN of GFF files")

args = argParser.parse_args()


def add_keys_nested_dict(d, file, repName, size):
    if file not in d:
        d[file] = {}
    if repName not in d[file]:
        d[file][repName] = size
    else:
        value = d[file][repName]
        d[file][repName] = value + size


repCounts = {}
fofn = open(args.fofn, 'r')
for file in fofn:
    file = file.rstrip()
    gff = open(file, 'r')
    for entry in gff:
        if not entry.startswith('#'):
            repeat = entry.split("\t")
            size = int(repeat[4])-int(repeat[3])
            id = repeat[8].split(" ")
            id[1] = id[1].replace("Motif:","")
            repName = id[1].replace("\"","")

            add_keys_nested_dict(repCounts, file, repName, size)
            value = repCounts[file][repName]

for key in repCounts:
    for key2 in repCounts[key]:
        value = repCounts[key][key2]
        print(key + "\t" + key2 + "\t" + str(value))
