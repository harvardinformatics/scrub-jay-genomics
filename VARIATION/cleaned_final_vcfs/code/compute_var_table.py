def calc_indel_length(ref,alts,aa,use_longest=True):
    alleles = [ref] + alts
    if aa != ".":
        base = aa
        var = [x for x in alleles if x != aa]
    else:
        base = ref
        var = alts

    biggest_alt = max([len(x) for x in var])
    smallest_alt = min([len(x) for x in var])

    if use_longest:
       return biggest_alt
    else:
       return smallest_alt

def calc_pop_numbers(genotypes, samples, alleles, aa):
    ## assume samples are listed as POP_SAMPLE, so split on underscore
    pop_dict = {}
    for sample in samples:
        pop = sample.split("_")[0]
        if pop not in pop_dict:
            pop_dict[pop] = {"derived": 0, "called": 0, "missing": 0}
        genotype = genotypes[samples.index(sample)].split("/") if "/" in genotypes[samples.index(sample)] else genotypes[samples.index(sample)].split("|")
        for hap in genotype:
            if hap == ".":
                pop_dict[pop]["missing"] += 1
            elif alleles[int(hap)] != aa and aa != ".":
                pop_dict[pop]["derived"] += 1
                pop_dict[pop]["called"] += 1
            elif int(hap) != 0 and aa == ".":
                pop_dict[pop]["derived"] += 1
                pop_dict[pop]["called"] += 1
            else:
                pop_dict[pop]["called"] += 1
    return pop_dict

def calc_subtype(type, polarized, base_allele_length, alt_len_min, alt_len_max, allele_count):
    max_bp = max(alt_len_min, alt_len_max, base_allele_length)
    min_bp = min(alt_len_min, alt_len_max, base_allele_length)
    if type == "SNP" and allele_count == 2:
        return "SNP"
    elif (type == "SNP" and allele_count > 2) or (type == "MNP"):
        return "SNP_Complex"
    elif allele_count == 2:
        if polarized and max_bp >= 50:
            if base_allele_length == 1 and alt_len_max >= 50:
                return "SVINS"
            elif base_allele_length >= 50 and alt_len_max == 1:
                return "SVDEL"
            elif base_allele_length > alt_len_max:
                return "SVDEL_Complex"
            elif base_allele_length < alt_len_max:
                return "SVINS_Complex"
            else:
                return "SV_Complex"
        if polarized and max_bp < 50:
            if base_allele_length == 1 and alt_len_max > 1:
                return "INS"
            elif base_allele_length > 1 and alt_len_max == 1:
                return "DEL"
            elif base_allele_length > alt_len_max:
                return "DEL_Complex"
            elif base_allele_length < alt_len_max:
                return "INS_Complex"
            else:
                return "INDEL_Complex"
        if not polarized and max_bp >= 50:
            return "SV_Complex"
        if not polarized and max_bp < 50:
            return "INDEL_Complex"
    elif allele_count > 2:
        if polarized and max_bp >= 50:
            if base_allele_length == 1 and alt_len_min >= 50:
                return "SVINS"
            elif base_allele_length >= 50 and alt_len_min == 1:
                return "SVDEL"
            elif base_allele_length > alt_len_min:
                return "SVDEL_Complex"
            elif base_allele_length < alt_len_max:
                return "SVINS_Complex"
            else:
                return "SV_Complex"
        if polarized and max_bp < 50:
            if base_allele_length == 1 and alt_len_min > 1:
                return "INS"
            elif base_allele_length > 1 and alt_len_min == 1:
                return "DEL"
            elif base_allele_length > alt_len_min:
                return "DEL_Complex"
            elif base_allele_length < alt_len_max:
                return "INS_Complex"
            else:
                return "INDEL_Complex"
        if not polarized and max_bp >= 50:
            return "SV_Complex"
        if not polarized and max_bp < 50:
            return "INDEL_Complex"
    else:
            return "Complex"

def process_file(input_file_path, output_file_path, regions, repeats, samples):
    with open(input_file_path, 'r') as infile, open(output_file_path, 'w') as outfile:
        outfile.write("chrom\tbedStart\tbedEnd\ttype\toverlap\trepeat\trepeat_family\tsubtype\tref\talt\taa\tinv\tpolarized\tbase_allele_len\talt_len_max\talt_len_min\tallele_count\tAW_DAC\tAW_AN\tAW_MISS\tAC_DAC\tAC_AN\tAC_MISS\tAI_DAC\tAI_AN\tAI_MISS\t" + "\t".join(samples) + "\n")
        for line in infile:
            elements = line.strip().split('\t')
            region_key = "-".join(elements[:3])
            if regions.get(region_key):
                elements[4] = regions[region_key]
            if repeats.get(region_key):
                repeats_info = repeats[region_key]
            else:
                repeats_info = "none\tnone"
            pos_info = elements[:5] + [repeats_info]
            ref = elements[5]
            alt = elements[6]
            aa = elements[7]
            inv = "NO" if elements[8] == "." else elements[8].upper()
            genotypes = elements[9:]
           
            ## compute allele info ##

            alleles = [ref] + alt.split(",")
            alts = alt.split(",")
            
            # set polarized flag
            polarized = False if aa == "." else True

            # compute indel length
            alt_len_max = calc_indel_length(ref,alts,aa,True)
            alt_len_min = calc_indel_length(ref,alts,aa,False)
            base_allele_length = len(aa) if aa != "." else len(ref)

            # compute allele count
            allele_count = len(alleles)

            subtype = calc_subtype(elements[3], polarized, base_allele_length, alt_len_min, alt_len_max, allele_count)

            # allele info
            allele_info = [str(subtype), str(ref), str(alt), str(inv), str(polarized), str(aa), str(base_allele_length), str(alt_len_max), str(alt_len_min), str(allele_count)]

            # process samples
            pop_info = calc_pop_numbers(genotypes, samples, alleles, aa)
            pop_info_output = []

            # for each pop calculate derived count, called count, missing count
            # POP = AW, AI, AC
            for pop in ["AW", "AC", "AI"]:
                pop_derived = pop_info[pop]["derived"]
                pop_called = pop_info[pop]["called"]
                pop_missing = pop_info[pop]["missing"]
                pop_info_output.append([str(pop_derived), str(pop_called), str(pop_missing)])
            
            # flatten pop_info_output
            pop_info_output = [item for sublist in pop_info_output for item in sublist]

            outfile.write('\t'.join(str(v) for v in pos_info) + '\t' + '\t'.join(str(a) for a in allele_info) + '\t' + '\t'.join(pop_info_output) + '\t' + '\t'.join(genotypes) + '\n')

def calc_genomic_region(region_file, split_csv=False):
    with open(region_file, 'r') as infile:
        region = {}
        for line in infile:
            elements = line.strip().split('\t')
            key = "-".join(elements[:3])
            elements[3] = "none,none" if elements[3] == "." else elements[3]
            if split_csv:
                if key not in region:
                    elements[3] = "\t".join(elements[3].split(",")) if elements[3] != "none,none" else "none\tnone"
                    region[key] = elements[3]
                else:
                    old_parts = region[key].split("\t")
                    new_parts = elements[3].split(",")
                    merged_first = ",".join(sorted(set([old_parts[0],new_parts[0]])))
                    merged_second = ",".join(sorted(set([old_parts[1],new_parts[1]])))
                    region[key] = "\t".join([merged_first, merged_second])
            else:
                if key not in region:
                    region[key] = ",".join(sorted(set(elements[3].split(","))))
                else:
                    parts = region[key].split(",")
                    parts = parts + elements[3].split(",")
                    region[key] = ",".join(sorted(set(parts)))
    return region

def get_samples(sample_file):
    with open(sample_file, 'r') as infile:
        samples = []
        for line in infile:
            samples.append(line.strip().split('\t'))
    # flatten samples
    samples = [item for sublist in samples for item in sublist]
    return samples

input_file = 'pggb_variation.tab'  # Replace with your input file path
output_file = 'variation_final_v2.tab'  # Replace with your desired output file path
region_file = 'pggb_variation_genomic_overlaps.tab' # genomic regions
repeats_file = 'pggb_variation_repeat_overlaps.tab' # repeats
sample_file = 'samples.txt' # sample file
regions = calc_genomic_region(region_file)
repeats = calc_genomic_region(repeats_file, split_csv=True)
samples = get_samples(sample_file)

process_file(input_file, output_file, regions, repeats, samples)
