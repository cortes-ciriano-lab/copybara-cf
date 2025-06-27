"""
Script to generate panel of normal for CN analysis normalisation
Created: 05/02/2025
Python 3.9.7
Carolin Sauer
"""

from multiprocessing import Pool
from multiprocessing import cpu_count
import sys
import pysam
import pybedtools
# import copy
import statistics
import math
# from scipy.stats import pearsonr
# import numpy as np
# from statsmodels.nonparametric.smoothers_lowess import lowess
# from scipy.interpolate import interp1d

def count_reads_in_curr_bin(bam, chrom, start, end, readcount_mapq):
    chunk_read_count = 0
    for read in bam.fetch(chrom, start, end, multiple_iterators = True):
        if read.is_secondary or read.is_supplementary or read.mapping_quality < readcount_mapq:
            continue
        else: 
            read_start = read.reference_start
            read_end = read.reference_end
            if read_start >= start and read_start <= end:
                    chunk_read_count += 1
            elif read_end >= start and read_end <= end:
                    chunk_read_count += 1
    return chunk_read_count

def binned_read_counting(curr_chunk, bed_chunk, aln_files, reader, readcount_mapq):
    print(f"    Read counting {curr_chunk} ...")
    # prepare pon list 
    # aln_files = []
    # with open(pon_list, "r") as file:
    #     for line in file:
    #         aln = line.strip()
    #         aln_files.append(aln)
    # print(f'PoN is being generated from a total of {len(aln_files)} alignment files')
    # # define alignment file type (bam or cram)
    # if aln_files[0].endswith('bam'):
    #     reader = "rb"
    # elif aln_files[0].endswith('cram'):
    #     reader = "rc"
    # else:
    #     sys.exit('Unrecognized file extension. Input files must be BAM/CRAM.')
    # initiate read counter
    chr_read_counts = []
    for bin in bed_chunk:
        chrom, start, end, gc, bases = bin[0], int(bin[1]), int(bin[2]), float(bin[3]), float(bin[4])
        bin_name = f"{chrom}:{start}_{end}"
        # # Filtering of bins by assigning each bin to boolean 'use' variable
        # if blacklisting == False:
        #     if bases_filter == False:
        #         use = True
        #     elif bases_filter == True:
        #         if bases >= bases_threshold:
        #             use = True
        #         else:
        #             use = False

        # elif blacklisting == True:
        #     blacklist = float(bin[5])
        #     if bases_filter == False:
        #         if blacklist <= bl_threshold:
        #             use = True
        #         else:
        #             use = False
        #     elif bases_filter == True:
        #         if blacklist <= bl_threshold and bases >= bases_threshold:
        #             use = True
        #         else:
        #             use = False

        # counting reads in all panel of normal bam files provided
        bin_read_count_collector = []
        for file in aln_files:
        # Open the BAM file using pysam
            bam = pysam.AlignmentFile(file, reader)
            cur_bam_bin = count_reads_in_curr_bin(bam, chrom, start, end, readcount_mapq)
            bin_read_count_collector.append(cur_bam_bin)
            # Close current bam file
            bam.close()
        read_count_median = statistics.median(bin_read_count_collector)
        # return output
        # chr_read_counts.append([bin_name, str(gc), str(use), str(read_count_mean)])
        chr_read_counts.append([bin_name, str(read_count_median)])
    return(chr_read_counts)

# def gc_loess(read_gc, read_outlier, gc_outlier, sample_size):
#     read_counts = [float(x[0]) for x in read_gc]
#     gc_content = [float(x[1]) for x in read_gc]
#     # define read and gc range to remove outliers for loess fitting 
#     read_range = np.quantile(read_counts, [0, 1 - read_outlier])
#     gc_range = np.quantile(gc_content, [gc_outlier, 1 - gc_outlier])
#     # subset data based on ranges
#     counts_subset = [x for x in read_gc if float(x[0]) > read_range[0] and float(x[0]) <= read_range[1] and
#                      float(x[1]) >= gc_range[0] and float(x[1]) <= gc_range[1]] 
#     read_subset = [float(x[0]) for x in counts_subset]
#     gc_subset = [float(x[1]) for x in counts_subset]
#     # Estimate correlation between GC content and read count to determine whether or not to perform GC correction
#     pear_res = pearsonr(read_subset, gc_subset)
#     pear_R, pear_pval = pear_res[0],pear_res[1]
#     print(f"gc bias: pearsons's R={round(pear_R,3)}; p_val={round(pear_pval,6)}")
#     if (abs(pear_R)>=0.1 and pear_pval < 0.05):
#         print("         GC bias detected. LOESS correction being performed.")
#         # Rough LOESS fit
#         if len(read_subset) >= sample_size:
#             ind = np.where(read_subset)[0]
#             random_ind = np.random.choice(ind, min(len(ind), sample_size), replace=False)
#             rough_fit = lowess(np.asarray(read_subset)[random_ind], np.asarray(gc_subset)[random_ind], frac=0.03)
#         else:
#             rough_fit = lowess(read_subset, gc_subset, frac=0.03)
#         # Refine the fit
#         gc_grid = np.linspace(0, 1, 1000)  # Interpolation points for refined LOESS
#         rough_interp = interp1d(rough_fit[:, 0], rough_fit[:, 1], bounds_error=False, fill_value=np.nan)
#         rough_pred = rough_interp(gc_grid)
#         # Perform the final LOESS fit
#         final_fit = lowess(rough_pred, gc_grid, frac=0.3)
#         # Interpolate the final fit
#         # final_interp = interp1d(final_fit[:, 0], final_fit[:, 1], bounds_error=False, fill_value=np.nan)
#         final_interp = interp1d(final_fit[:, 0], final_fit[:, 1], bounds_error=False, fill_value="extrapolate")
#         final_pred = final_interp(gc_content)  # Predicted values for all bins
#         # Correct read counts
#         cor_gc = read_counts / final_pred # GC correct read counts
#         ind = np.where(np.isfinite(cor_gc))[0]
#         pear_res = pearsonr(cor_gc[ind], np.asarray(gc_content)[ind])
#         pear_R, pear_pval = pear_res[0],pear_res[1]
#         print(f"gc bias after correction: pearsons's R={round(pear_R,3)}; p_val={round(pear_pval,6)}")
#         cor_gc = list(cor_gc)
#     else:
#         cor_gc = read_counts
#     return cor_gc

# def gc_correct_counts(filtered_counts, read_outlier = 0.01, gc_outlier = 0.001, sample_size = 50000):
#     # extract read counts and gc content from filtered counts
#     read_gc = [[float(x[-1]),float(x[1])] for x in filtered_counts]
#     cor_gc = gc_loess(read_gc, read_outlier, gc_outlier, sample_size)
#     # Prepare output
#     gc_cor_counts = copy.deepcopy(filtered_counts)
#     for id,r in enumerate(gc_cor_counts):
#         r[-1] = float(cor_gc[id])
#     return gc_cor_counts

# def filter_correct(countData):
#     """ perform filtering, normalisation and transformation """
#     ## filtering: Remove bins with 0 reads (no sequencing data in this region). - need to be excluded for segmentation and log2 transformation
#     # filtered_counts = [x for x in countData if x[-2] == 'True' and float(x[-1]) != 0]
#     # gc_cor_counts = gc_correct_counts(filtered_counts) # gc correct raw read counts
#     # return out
#     # counts_out = [[str(x[0]), str(x[-1])] for x in gc_cor_counts]
#     counts_out = [[str(x[0]), str(x[-1])] for x in countData]
#     return counts_out

def chunkify_bed(bed_file, chunk_size):
    """ Divides bed files into chunks based on threads available for multiprocessing."""
    chunks = []
    current_chunk = []
    for i, feature in enumerate(bed_file):
        current_chunk.append(feature)
        if (i + 1) % chunk_size == 0:
            chunks.append(pybedtools.BedTool(current_chunk))
            current_chunk = []
    if current_chunk:
        chunks.append(pybedtools.BedTool(current_chunk))
    return chunks

def count_reads(outdir, pon_list, pon_name, bin_annotations_path, readcount_mapq, threads):
    """ Perform binned read counting on bam file/files per chromosome """

    # check and define threads
    new_threads = min(threads, cpu_count())
    print(f"... Bin read counter will use threads = {new_threads}. (threads = {threads} defined; threads = {cpu_count()} available) ...")
    threads = new_threads

    aln_files = []
    with open(pon_list, "r") as file:
        for line in file:
            aln = line.strip()
            aln_files.append(aln)
    print(f'PoN is generated from a total of {len(aln_files)} alignment files')
    # define alignment file type (bam or cram)
    if aln_files[0].endswith('bam'):
        reader = "rb"
    elif aln_files[0].endswith('cram'):
        reader = "rc"
    else:
        sys.exit('Unrecognized file extension. Input files must be BAM/CRAM.')

    # Define bed_file chunks
    bed_file =  pybedtools.BedTool(bin_annotations_path)
    bed_length = sum(1 for _ in bed_file)
    chunk_size = bed_length // threads + (bed_length % threads > 0)
    print(f"    splitting bed file into n = {math.ceil(bed_length / chunk_size)} chunks ...")
    # Split the BED file into chunks
    chunks = chunkify_bed(bed_file, chunk_size)

    # only use multiprocessing if more than 1 thread available/being used.
    if threads == 1:
        # loop through chromosomes
        print("multithreading skipped.")
        countData = []
        for idx,bed_chunk in enumerate(chunks):
            curr_chunk = f"chunk {idx+1}"
            counts_binned_chr = binned_read_counting(curr_chunk, bed_chunk, aln_files, reader, readcount_mapq)
            countData.append(counts_binned_chr)
        countData = [x for xs in countData for x in xs]

    else:
        print(f"multithreading using {threads} threads.")
        args_in = [[str(f"chunk {idx+1}"),bed_chunk, aln_files, reader, readcount_mapq] for idx,bed_chunk in enumerate(chunks)]
        # print(args_in)
        with Pool(processes=threads) as pool:
            countData = [x for xs in list(pool.starmap(binned_read_counting, args_in)) for x in xs]

    # counts_out = filter_correct(countData=countData)

    #----
    # 4. Get results and write out
    #----
    pon_read_counts = f"{outdir}/{pon_name}_read_counts.tsv"
    outfile = open(pon_read_counts, "w")
    header=['bin', 'PoN_read_count']
    outfile.write('\t'.join(header)+'\n')
    for r in countData:
        Line = '\t'.join(r) + '\n'
        outfile.write(Line)
    outfile.close()

    return pon_read_counts

if __name__ == "__main__":
    print("PoN generator")