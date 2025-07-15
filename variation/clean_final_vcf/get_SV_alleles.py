import argparse
import re

argParser = argparse.ArgumentParser()
argParser.add_argument("-i", "--input", dest='tab', help=".tab summary file of SVs")
argParser.add_argument("-m", "--min", dest='minSize', help="minimum size SV to keep")

args = argParser.parse_args()

tab = open(args.tab, 'r')

for line in tab:
    line = line.rstrip()
    sv = line.split("\t")
    scaff = sv[0]
    scaffShort = re.sub('aphWoo1#REF#','',scaff)
    scaffShort = re.sub('.*ch', 'ch', scaffShort)
 #   print(scaff+"\t"+scaffShort, file=open('scaffID_shortened_key.txt','a'))
    start = sv[1]
    end = sv[2]
    svType = sv[3]
    refAl = sv[8]
    altAl = sv[9]
    #make arrays to hold multiple alleles
    refMulti = refAl.split(",")
    altMulti = altAl.split(",")

    refLongest = max(refMulti, key=len)
    altLongest = max(altMulti, key=len)
  #  print(refLongest+"\t"+altLongest)
  #  print(len(refLongest)+len(altLongest))

    #If reference allele at least threshold size and longer than alt:
    if len(refLongest) >= int(args.minSize) and len(refLongest) > len(altLongest):
        print(">"+scaffShort+":"+start+"-"+end+"\n"+refLongest)
    #If alt allele > threshold size and longer than ref
    elif len(altLongest) >= int(args.minSize) and len(altLongest) > len(refLongest):
        print(">"+scaffShort+":"+start+"-"+end+"\n"+altLongest)
