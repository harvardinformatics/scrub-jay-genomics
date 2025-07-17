import gffutils
import argparse
import os

# parser = argparse.ArgumentParser(description="Replace transcript_id with oId from ref")
# parser.add_argument('--gtf', type=str, help='Reference GFF/GTF file')

# args = parser.parse_args()

#refFile = open(args.gtf, 'r')
refFile = 'merged_TOGA_noRef_noDup.combined.gtf'

#Create the gffutils DBs if they do not exist already:
# if os.path.exists('refAnnotation.db'):
#     refDB = gffutils.FeatureDB('refAnnotation.db')
# else:
#     refDB = gffutils.create_db(refFile, dbfn="refAnnotation.db", force=True, keep_order=True, merge_strategy="merge", sort_attribute_values=True, disable_infer_transcripts=True)

refDB = gffutils.create_db(refFile, ':memory:', force=True, keep_order=True, merge_strategy="merge", sort_attribute_values=True, disable_infer_transcripts=True)

#MANUALLY making gtf dialect dictionary, for printing
gtfDia = {'field separator': '; ',
 'fmt': 'gtf',
 'keyval separator': ' ',
 'leading semicolon': False,
 'multival separator': ',',
 'quoted GFF2 values': True,
 'repeated keys': False,
 'trailing semicolon': True}


for transcript in refDB.features_of_type('transcript'):
    attrs = dict(transcript.attributes)
    
    newLabel = attrs['oId']
    #save the old label just in case
    attrs['gffcompare_id'] = attrs['transcript_id']
    #Fix the transcript label
    attrs['transcript_id'] = newLabel

     #build feature to manually add:
    f = gffutils.Feature(
        seqid=transcript.chrom,
        source='gffutils_fix',
        featuretype=transcript.featuretype,
        start=transcript.start,
        end=transcript.end,
        strand=transcript.strand,
        attributes=attrs,
        dialect=gtfDia
    )
    print(f)

    #fix the baby...
    for exon in refDB.children(transcript, featuretype='exon', order_by='start'):
        exonAttrs = dict(exon.attributes)
        exonAttrs['gffcompare_id'] = exonAttrs['transcript_id']
        #Fix the transcript label
        exonAttrs['transcript_id'] = newLabel

        e = gffutils.Feature(
            seqid=exon.chrom,
            source='gffutils_fix',
            featuretype=exon.featuretype,
            start=exon.start,
            end=exon.end,
            strand=exon.strand,
            attributes=exonAttrs,
            dialect=gtfDia
        )
        print(e)