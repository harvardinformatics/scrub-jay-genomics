import gffutils
import argparse
import os
from collections import defaultdict

# parser.add_argument('--ref', type=str, help='Reference GFF/GTF file')
# parser.add_argument('--query', type=str, help='Query GFF/GTF file')

# args = parser.parse_args()

# refFile = open(args.ref, 'r')
# queryFile = open(args.query, 'r')

refFile = 'both_ref.gtf'
queryFile = 'merged_TOGA_noRef_noDup.combined.labelFix.gtf'

orthoFile = open('combined_orthology_classification.noMissing.tsv', 'r')
orthoDict = defaultdict(dict)

#Make dict from ortho file where transID = gene name (from ref spp)
for trans in orthoFile:
    #Populate orthology hash
    #Key is the transcript ID in the TOGA output
    transID = trans.split('\t')[3]

    orthoDict[transID] = trans.split('\t')[0]

# #Create the gffutils DBs if they do not exist already:
# if os.path.exists('refAnnotation.db'):
#     refDB = gffutils.FeatureDB('refAnnotation.db')
# else:
#     refDB = gffutils.create_db(refFile, dbfn="refAnnotation.db", force=True, keep_order=True, merge_strategy="merge", sort_attribute_values=True, disable_infer_transcripts=True, disable_infer_genes=True)

# if os.path.exists('queryAnnotation.db'):
#     queryDB = gffutils.FeatureDB('queryAnnotation.db')
# else:
#     queryDB = gffutils.create_db(queryFile, dbfn="queryAnnotation.db", force=True, keep_order=True, merge_strategy="merge", sort_attribute_values=True, disable_infer_transcripts=True)


refDB = gffutils.FeatureDB('refAnnotation.db')
queryDB = gffutils.FeatureDB('queryAnnotation.db')

#MANUALLY making gtf dialect dictionary, for printing
gtfDia = {'field separator': '; ',
 'fmt': 'gtf',
 'keyval separator': ' ',
 'leading semicolon': False,
 'multival separator': ',',
 'quoted GFF2 values': True,
 'repeated keys': False,
 'trailing semicolon': True}

#gffutils.constants.always_return_list = False
#Make a check string to avoid annoying printing gene multiple times
geneNameCheck = "bingus"

for querTrans in queryDB.features_of_type('transcript'):
    #print its parent gene annotation (INFERRED by gffutils during db construction)
    for gene in queryDB.parents(querTrans, featuretype='gene', order_by='start'):
        #Get name of current gene
        currentGeneID = gene.id

        #As long an new gene (have not printed yet...)
        if currentGeneID != geneNameCheck:
            attrs = dict(gene.attributes)
            #make empty list to hold different source gene names (in case multiple)
            sourceGeneList = []

            #Goign to modify by adding gene_name to attr; need to get all childen to gene first
            for childTrans in queryDB.children(gene, featuretype='transcript', order_by='start'):
                childID = childTrans.id
                if childID in orthoDict:
                    sourceGeneID = orthoDict[childID]
                    if sourceGeneID not in sourceGeneList:
                        sourceGeneList.append(sourceGeneID)
            #Once gone thru all the transcripts for that gene, add the list of gene names to attr dict
            attrs['gene_name'] = sourceGeneList

            #build modified gene feature manually and print it:
            gf = gffutils.Feature(
                seqid=gene.chrom,
                source='gffutils_merge',
                featuretype=gene.featuretype,
                start=gene.start,
                end=gene.end,
                strand=gene.strand,
                attributes=attrs,
                dialect=gtfDia
            )
            print(gf)

        #rename for next cycle
        geneNameCheck = currentGeneID

    #output transcript
    print(querTrans)

    #output the babies
    for child in queryDB.children(querTrans, featuretype='exon', order_by='start'):
        print(child)

    #Get matching transcript from ref...
    refTrans = refDB[querTrans.id]
    #Get all the CDS annotations for it
    for refFeat in refDB.children(refTrans, featuretype=('CDS', 'stop_codon', 'start_codon'), order_by='start'):
        #change the attributes of the CDS from ref; first pull out as dict
        attrs = dict(refFeat.attributes)
        #Modify the gene ID to be XLOC from query
        attrs['gene_id'] = querTrans.attributes['gene_id']

        #build feature to manually add:
        f = gffutils.Feature(
            seqid=refFeat.chrom,
            source='gffutils_merge',
            featuretype=refFeat.featuretype,
            start=refFeat.start,
            end=refFeat.end,
            strand=refFeat.strand,
            attributes=attrs,
            dialect=gtfDia
        )
        #output CDS
        print(f)

