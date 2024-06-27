# import bisect
# import csv
# import string
# import pysam
# import numpy as np
# import vcf
import cyvcf2
import sys

phased_vcf=sys.argv[1]

vcf_reader = cyvcf2.VCF(phased_vcf)
hets = []
for variant in vcf_reader:
    # if variant.is_snp == True and len(record.get_hets()) > 0:# and record.genotype(sample).phased is True:
    # check if variant is a SNP
    if variant.is_snp:
        # iterate through each genotype
        for s_idx, gt in enumerate(variant.genotypes):
            # only get heterozygous snps
            if gt[0] != gt[1]:
                # print(s_idx, gt)
                sample = vcf_reader.samples[s_idx]
                ps = int(variant.format('PS')[s_idx]) if 'PS' in variant.FORMAT else None
                gt_str = f"{gt[0]}|{gt[1]}" if gt[2] == 1 else f"{gt[0]}/{gt[1]}"
                # dp = int(variant.format('DP')[s_idx]) if 'DP' in variant.FORMAT else None
                # ad = int(variant.format('AD')[s_idx]) if 'AD' in variant.FORMAT else None
                var_out = [variant.CHROM,str(variant.POS),str(variant.POS),variant.REF,variant.ALT[0],gt_str,str(ps)]
                hets.append(var_out)
                # print(f"{variant.CHROM}\t{variant.POS}\t{variant.POS}\t{variant.REF}\t{variant.ALT[0]}\t{gt_str}\t{ps}")



