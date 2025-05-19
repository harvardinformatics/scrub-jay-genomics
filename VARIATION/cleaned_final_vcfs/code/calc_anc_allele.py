def calc_ancestral_allele(line):
    ## Calculate the ancestral allele

    # Split genotype (line[4]) on slash or pipe
    haps = line[4].split('/') if '/' in line[4] else line[4].split('|')

    # If the genotype is homozygous and not missing, the ancestral allele is the same as the genotype
    if haps[0] == haps[1] and haps[0] != '.':
        aa = int(haps[0])
    # If both missing, the ancestral allele is missing
    elif haps[0] == '.' and haps[1] == '.':
        return '.'
    # If one is missing, the ancestral allele is the other one
    elif haps[0] == '.':
        aa = int(haps[1])
    elif haps[1] == '.':
        aa = int(haps[0])
    # If the genotype is heterozygous, the ancestral allele is missing
    else:
        return '.'

    # now get the genotypes
    alts = line[3].split(',')
    ref = [line[2]]
    genotypes = ref + alts
    return genotypes[aa]


def process_file(input_file_path, output_file_path):
    with open(input_file_path, 'r') as infile, open(output_file_path, 'w') as outfile:
        for line in infile:
            elements = line.strip().split('\t')
            first_four_columns = elements[:4]  # Extract the first four columns
            result = calc_ancestral_allele(elements)
            outfile.write('\t'.join(first_four_columns) + '\t' + result + '\n')

# Example usage
input_file = 'aa.tab'  # Replace with your input file path
output_file = 'aa.processed.tab'  # Replace with your desired output file path
process_file(input_file, output_file)
