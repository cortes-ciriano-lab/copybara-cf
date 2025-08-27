"""
Script to detect/analyse focal amplifications/ecDNA
Created: 18/08/2025
Python 3.9.7
Carolin Sauer
"""

import pysam
import random
from intervaltree import IntervalTree
import re

#----
# Define functions
#----

def sample_region(size, chrom_sizes, blacklists, seed_val):
    '''Function to randomly select random regions of equal size across the genome.'''
    chosen = {chrom: IntervalTree() for chrom in chrom_sizes}
    random.seed(seed_val)
    while True:
        chrom = random.choices(list(chrom_sizes.keys()), weights=chrom_sizes.values())[0]
        max_start = chrom_sizes[chrom] - size
        start = random.randint(0, max_start)
        end = start + size
        # Check blacklist + already chosen
        if not blacklists[chrom].overlaps(start, end) and not chosen[chrom].overlaps(start, end):
            chosen[chrom].addi(start, end)
            # calculate gc content 
            return chrom, start, end
        
def chr_key(chr_name):
    '''function to allow sorting of chromosomes'''
    # If it's already an integer, return it directly
    if isinstance(chr_name, int):
        return chr_name
    # If it's a string, normalize
    m = re.match(r'chr(\d+|X|Y|M)', chr_name)
    if m:
        val = m.group(1)
        if val == 'X':
            return 23
        elif val == 'Y':
            return 24
        elif val == 'M':
            return 25
        else:
            return int(val)
    return 999  # in case of unusual names


def generate_bins(outdir, sample, ref, roi, roi_buffer, n_regions, chromosomes, blacklist): #, threads):
    '''
    Main function to process fasta file and generate bins for focal analysis based on regions of interest
    '''
    # # check and define threads
    # new_threads = min(threads, cpu_count())
    # print(f"... Bin generator will use threads = {new_threads}. (threads = {threads} defined; threads = {cpu_count()} available) ...")
    # threads = new_threads
    # Extract chormosome info and pass fasta file into chr_in list
    fasta = pysam.FastaFile(ref)
    # contigs = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']
    contigs = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22']
    ref_contigs = fasta.references
    if 'chr1' in ref_contigs:
        contigs = [f'chr{x}' for x in contigs]
    if chromosomes != 'all':
        chr_names = [contigs[(int(x)-1)] for x in chromosomes]
    else:
        chr_names = contigs
    #   define chrom sizes
    chrom_sizes = {}
    for chrom in chr_names:
        chrom_sizes[chrom] = fasta.get_reference_length(chrom)
    # chr_in = [[chrom,fasta.get_reference_length(chrom),ref,bin_size,blacklist,blacklist_buffer] for chrom in chr_names]
    # fasta.close()
    roi_name = roi[3]
    bin_size = roi[4]
    # prep blacklist
    blacklisted_regions = {chrom: IntervalTree() for chrom in chrom_sizes}
    # add in blacklist to only include regions without blacklist overlap if provided
    if (blacklist is not None):
        with open(blacklist) as f:
            for line in f:
                chrom, start, end = line.strip().split()[:3]
                if chrom in blacklisted_regions.keys():
                    blacklisted_regions[chrom].addi(int(start), int(end))
    else:
        pass
    blacklisted_regions[roi[0]].addi(int(roi[1])-roi_buffer, int(roi[2])+roi_buffer)
    # randomly choose regions
    regions = [sample_region(bin_size, chrom_sizes, blacklisted_regions, _) for _ in range(n_regions)]
    # add in region of interest
    regions.append((roi[0], roi[1], roi[2]))
    # sort regions
    regions = sorted(regions, key=lambda r: (chr_key(r[0]), r[1], r[2]))
    # prepare output
    List_of_bins = []
    for r in regions:
        chrom,start,end=r
        overlap_pct = float(0)
        chr_bin_seq = fasta.fetch(chrom, start, end)
        acgt_len = 0
        gc_len = 0
        for base in chr_bin_seq:
            if base.upper() in 'ACTG':
                acgt_len += 1
            if base.upper() in 'CG':
                gc_len += 1
        chr_bases = acgt_len / bin_size * 100
        gc_content = gc_len / bin_size
        #out
        if blacklist is not None:
            List_of_bins.append([chrom, str(start), str(end), str(gc_content), str(chr_bases), str(overlap_pct)])
        else:
            List_of_bins.append([chrom, str(start), str(end), str(gc_content), str(chr_bases)])
    # Concat results into a single bed file
    outfile_name = f"{outdir}/focal_regions_{sample}.bed"
    outfile = open(outfile_name, "w")
    for r in List_of_bins:
        Line = '\t'.join(r) + '\n'
        outfile.write(Line)
    outfile.close()
    #
    return outfile_name
