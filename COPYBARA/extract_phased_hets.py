import bisect
import csv
import string
import pysam
import numpy as np
import vcf
import sys

def loadcsv(filename, criterion1, criterion2, repeat_units):
    with open(filename, "rb") as csvfile:
        datareader = csv.reader(csvfile, delimiter="\t")
        for row in datareader:
            if int(row[5]) >= criterion1 and int(row[5]) <= criterion2 and \
                            int(row[4]) in repeat_units:
                yield row


#sample=sys.argv[1]
phased_vcf=sys.argv[1]
#cohort=sys.argv[2]


import vcf
vcf_reader = vcf.Reader(open(phased_vcf, 'r'))

sample=vcf_reader.samples[0]

for record in vcf_reader:
    now = record.get_hets()
    if record.is_snp == True and len(record.get_hets()) > 0:# and record.genotype(sample).phased is True:
        for sample in record.samples:
            try:
                try:
                    PS=sample['PS']
                except:
                    PS=""
                    pass
                GT=sample['GT']
                DP=sample['DP']
                AD=sample['AD']
                print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(record.CHROM, record.POS, record.POS, record.REF, record.ALT[0], GT, PS))#, record.INFO["PS"]))
            except:
                pass 


