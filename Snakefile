
rule hifi:
    input:
        "../RAW_DATA/{reads}.fastq.gz"
    output:
        "results/hifiasm_output_RAW/{sample}.ctg.gfa"
    shell:
        "../../hifiasm/hifiasm -t $THREADS -o {output} {input}" #how do threads??

rule gfa_to_fasta:
    input:
        "results/hifiasm_output_RAW/{sample}.ctg.gfa"
    output:
        "results/assemblies/{sample}.fa.gz"
    shell:
        "awk '/^S/{print \">\"$2;print $3}' {input} | gzip > {output}"

rule quast:
    input:

    output:

    shell:
