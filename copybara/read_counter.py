"""
Script to count reads copy number estimation
Created: 12/09/2023
Python 3.9.7
Carolin Sauer
"""

from multiprocessing import Pool
from multiprocessing import cpu_count
import pysam
import pybedtools
import copy
import statistics
import math
import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess

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

def binned_read_counting(curr_chunk, bed_chunk, alns, nmode, blacklisting, bl_threshold, bases_filter, bases_threshold, readcount_mapq):
    print(f"    Read counting {curr_chunk} ...")

    # Open the BAM file using pysam
    bam_T = pysam.AlignmentFile(alns['tumour'], "rb")
    bam_N = pysam.AlignmentFile(alns['normal'], "rb") if nmode == "mnorm" and len(alns) == 2 else None
    # PoN = pysam.AlignmentFile(alns['pon'], "r") if nmode == "pon" and len(alns) == 2 else None
    if nmode == "pon" and len(alns) == 2:
        PoN = []
        with open(alns['pon'], "r") as file:
            next(file)
            for line in file:
                fields = line.strip().split("\t")
                PoN.append(fields)
   
    # initiate read counter
    chr_read_counts = []
    for bin in bed_chunk:
        chrom, start, end, gc, bases = bin[0], int(bin[1]), int(bin[2]), float(bin[3]), float(bin[4])
        bin_name = f"{chrom}:{start}_{end}"
        # Filtering of bins by assigning each bin to boolean 'use' variable
        if blacklisting == False:
            if bases_filter == False:
                use = True
            elif bases_filter == True:
                if bases >= bases_threshold:
                    use = True
                else:
                    use = False

        elif blacklisting == True:
            blacklist = float(bin[5])
            if bases_filter == False:
                if blacklist <= bl_threshold:
                    use = True
                else:
                    use = False
            elif bases_filter == True:
                if blacklist <= bl_threshold and bases >= bases_threshold:
                    use = True
                else:
                    use = False
        # counting reads in tumour (and if provided) normal bam files in each bin
        chunk_read_count_T = count_reads_in_curr_bin(bam_T, chrom, start, end, readcount_mapq)
        chunk_read_count_N = count_reads_in_curr_bin(bam_N, chrom, start, end, readcount_mapq) if bam_N else None

        if blacklisting == True:
            chr_read_counts.append([bin_name, chrom, str(start), str(end), str(gc), str(bases), str(blacklist), str(use), str(chunk_read_count_T), str(chunk_read_count_N)])
        else:
            chr_read_counts.append([bin_name, chrom, str(start), str(end), str(gc), str(bases), str(use), str(chunk_read_count_T), str(chunk_read_count_N)])

        # replace None values (lack of matched normal) with PoN count values if provided. PoN will need to be generated seperately, using the same bin size.
        if nmode == "pon": 
            print('         PoN read counts being added for normalisation')
            if len(PoN) != len(chr_read_counts):
                print('ERROR: PoN will have to be generated using the same bin size as used for read counting. See documentation on how to generate PoN.')
            if len(PoN) == len(chr_read_counts):
                for id,r in enumerate(chr_read_counts):
                    r[-1] = PoN[id][-1]

    # Close BAM file and return out
    bam_T.close()
    if bam_N is not None:
        bam_N.close()
    return(chr_read_counts)

# def correct_counts(gc_content,counts,smoothness=0.1):
#     # Smoothing spline
#     spline = UnivariateSpline(gc_content, counts, s=smoothness)
#     # Calculate expected counts based on the spline
#     # expected_counts_interim = spline(gc_content)
#     expected_counts = spline(gc_content)
#     print(expected_counts)
#     # # Avoid division by zero
#     # expected_counts = []
#     # for value in expected_counts_interim:
#     #     if value == 0:
#     #         value = float('NaN')
#     #     expected_counts.append(value)
#     # expected_counts[expected_counts == 0] = float('NaN')
#     # Correct read counts
#     corrected_counts = []
#     for id,val in enumerate(counts):
#         corrected = val / expected_counts[id]
#         corrected_counts.append(corrected)
#     return corrected_counts


# def correct_gc_bias(nmode,countData,smoothness=0.1):
#     '''
#     Model GC bias to estimate GC correction for each bin and correct raw read counts
#     '''
#     # add 'row number' to start of count data for later sorting
#     for id,r in enumerate(countData):
#         r.insert(0,id)
#     # sort by gc content
#     countData.sort(key = lambda x: float(x[5]))
#     # countData.sort(key = lambda x: (float(x[5]),float(x[-2])))
#     #
#     # Prepare input for smoothing spline to model GC bias
#     gc_content = [float(x[5]) for x in countData]
#     t_counts = [int(x[-2]) for x in countData]
#     #
#     corrected_counts_t = correct_counts(gc_content,t_counts,smoothness)
#     #
#     # replace count values with corrected values
#     for id,r in enumerate(countData):
#         r[-2] = str(corrected_counts_t[id])
#     # 
#     # repeat same for normal counts (if available)
#     if nmode == "mnorm" or nmode == "pon":
#         n_counts = [int(x[-1]) for x in countData]
#         corrected_counts_n = correct_counts(gc_content,n_counts,smoothness)
#         # replace count values with corrected values
#         for id,r in enumerate(countData):
#             r[-1] = str(corrected_counts_n[id])
#     #
#     # undo sorting by gc content
#     countData.sort(key = lambda x: x[0])
#     countData = [x[1:] for x in countData]
#     return countData

# def filter_and_normalise(nmode, countData):
#     """ perform filtering, normalisation and transformation """
#     ## filtering: Remove bins with 0 reads (no sequencing data in this region). - need to be excluded for segmentation and log2 transformation
#     if nmode == "self":
#         filtered_counts = [x for x in countData if x[-3] == 'True' and float(x[-2]) != 0 and x[-2] != 'nan']
#         med_self = statistics.median([int(x[-2]) for x in filtered_counts]) #estimate genome wide median for selfnormalisation
#         # Normalise and log2 transform
#         normalised_counts = copy.deepcopy(filtered_counts)
#         for r in normalised_counts:
#             r[-2] = str(math.log2(int(r[-2])/med_self)) # median normalise readcounts
#             del r[-1]
#     elif nmode == "mnorm" or nmode == "pon":
#         filtered_counts = [x for x in countData if x[-3] == 'True' and int(x[-2]) != 0 and int(x[-1]) != 0]
#         # cov_scaler = statistics.median([math.log2(int(x[-2])/int(x[-1])) for x in filtered_counts])
#         cov_scaler = math.log2(statistics.median([(int(x[-2])/int(x[-1])) for x in filtered_counts]))
#         normalised_counts = copy.deepcopy(filtered_counts)
#         for r in normalised_counts:
#             n = str(math.log2(int(r[-2])/int(r[-1])) - cov_scaler)
#             del r[-2:]
#             r.append(n)
#     # return out
#     return filtered_counts, normalised_counts
def gc_loess(read_gc, read_outlier, gc_outlier):
    read_counts = [int(x[0]) for x in read_gc]
    gc_content = [float(x[1]) for x in read_gc]
    # define read and gc range to remove outliers for loess fitting 
    read_range = np.quantile(read_counts, [0, 1 - read_outlier])
    gc_range = np.quantile(gc_content, [gc_outlier, 1 - gc_outlier])
    # subset data based on ranges
    counts_subset = [x for x in read_gc if int(x[0]) >= read_range[0] and int(x[0]) <= read_range[1] and
                     float(x[1]) >= gc_range[0] and float(x[1]) <= gc_range[1]] 
    read_subset = [int(x[0]) for x in counts_subset]
    gc_subset = [float(x[1]) for x in counts_subset]
    # Rough LOESS fit
    rough_loess = lowess(read_subset, gc_subset, frac=0.03, return_sorted=False)
    # Refine the fit
    gc_grid = np.linspace(0, 1, 1000)  # Interpolation points for refined LOESS
    rough_pred = np.interp(gc_grid, np.sort(gc_subset), np.sort(rough_loess))
    final_loess = lowess(rough_pred, gc_grid, frac=0.3, return_sorted=False)
    # Correct read counts
    smooth_gc = np.interp(gc_content, gc_grid, final_loess)  # Predicted values for all bins
    cor_gc = read_counts / smooth_gc # GC correct read counts
    return cor_gc

def gc_correct_counts(filtered_counts, nmode, read_outlier = 0.01, gc_outlier = 0.001):
    # extract read counts and gc content from filtered counts
    # tumour reads
    read_gc = [[int(x[-2]),float(x[4])] for x in filtered_counts]
    cor_gc = gc_loess(read_gc, read_outlier, gc_outlier)
    # Prepare output
    gc_cor_counts = copy.deepcopy(filtered_counts)
    for id,r in enumerate(gc_cor_counts):
        r[-2] = float(cor_gc[id])
    # repeat same for normal read counts if nmode = "mnorm" or "pon"
    if nmode == "mnorm" or nmode == "pon":
        read_gc = [[int(x[-1]),float(x[4])] for x in filtered_counts]
        cor_gc = gc_loess(read_gc, read_outlier, gc_outlier)
        for id,r in enumerate(gc_cor_counts):
            r[-1] = float(cor_gc[id])
    return gc_cor_counts

def filter_correct_normalise(nmode, countData):
    """ perform filtering, normalisation and transformation """
    ## filtering: Remove bins with 0 reads (no sequencing data in this region). - need to be excluded for segmentation and log2 transformation
    if nmode == "self":
        filtered_counts = [x for x in countData if x[-3] == 'True' and float(x[-2]) != 0]
        gc_cor_counts = gc_correct_counts(filtered_counts, nmode) # gc correct raw read counts
        med_self = statistics.median([float(x[-2]) for x in gc_cor_counts]) #estimate genome wide median for selfnormalisation
        # Normalise and log2 transform
        normalised_counts = copy.deepcopy(gc_cor_counts)
        for r in normalised_counts:
            r[-2] = str(math.log2(float(r[-2])/med_self)) # median normalise readcounts
            del r[-1]
    elif nmode == "mnorm" or nmode == "pon":
        filtered_counts = [x for x in countData if x[-3] == 'True' and int(x[-2]) != 0 and int(x[-1]) != 0]
        # cov_scaler = statistics.median([math.log2(int(x[-2])/int(x[-1])) for x in filtered_counts])
        cov_scaler = math.log2(statistics.median([(int(x[-2])/int(x[-1])) for x in filtered_counts]))
        normalised_counts = copy.deepcopy(filtered_counts)
        for r in normalised_counts:
            n = str(math.log2(int(r[-2])/int(r[-1])) - cov_scaler)
            del r[-2:]
            r.append(n)
    # return out
    return normalised_counts

def chunkify_bed(bed_file, chunk_size):
    '''
    Divides bed files into chunks based on threads available for multiprocessing.
    '''
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

def count_reads(outdir, tumour, normal, panel_of_normals, sample, bin_annotations_path, readcount_mapq, blacklisting, bl_threshold, bases_filter, bases_threshold, threads):
    """ Perform binned read counting on bam file/files per chromosome """

    # check and define threads
    new_threads = min(threads, cpu_count())
    print(f"... Bin read counter will use threads = {new_threads}. (threads = {threads} defined; threads = {cpu_count()} available) ...")
    threads = new_threads

    if normal is not None and panel_of_normals is None:
        nmode = "mnorm" #matched normal will be used for logR normalisation
        aln_files = {
                'tumour': tumour,
                'normal': normal
            }
    elif normal is None and panel_of_normals is not None:
        nmode = "pon" #panel of normals will be used for logR normalisation
        aln_files = {
                'tumour': tumour,
                'pon': panel_of_normals
            }
    elif normal is None and panel_of_normals is None:
        nmode = "self"
        aln_files = {
                'tumour': tumour
            }

    print(f"Normalisation mode: {nmode}")

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
            counts_binned_chr = binned_read_counting(curr_chunk, bed_chunk, aln_files, nmode, blacklisting, bl_threshold, bases_filter, bases_threshold, readcount_mapq)
            countData.append(counts_binned_chr)
        countData = [x for xs in countData for x in xs]

    else:
        print(f"multithreading using {threads} threads.")
        args_in = [[str(f"chunk {idx+1}"),bed_chunk,aln_files, nmode, blacklisting, bl_threshold, bases_filter, bases_threshold, readcount_mapq] for idx,bed_chunk in enumerate(chunks)]
        # print(args_in)
        with Pool(processes=threads) as pool:
            countData = [x for xs in list(pool.starmap(binned_read_counting, args_in)) for x in xs]

    # countData_corrected = correct_gc_bias(nmode,countData,smoothness=0.1)

    # ##################################
    # # test write out for DEV ###
    # outfile = open(f"{outdir}/{sample}_CORRECTED_read_counts.tsv", "w")
    # if blacklisting == True:
    #     header=['bin', 'chromosome','start','end', 'gc_content', 'known_bases', 'overlap_blacklist', 'use_bin', 'tumour_read_count_cor', 'normal_read_count_cor']
    # else: 
    #     header=['bin', 'chromosome','start','end', 'gc_content', 'known_bases', 'use_bin', 'tumour_read_count_cor', 'normal_read_count_cor']
    # outfile.write('\t'.join(header)+'\n')
    # for r in countData_corrected:
    #     Line = '\t'.join(r) + '\n'
    #     outfile.write(Line)
    # outfile.close()
    # ##################################


    normalised_counts = filter_correct_normalise(nmode=nmode, countData=countData)

    #----
    # 4. Get results and write out
    #----
    outfile = open(f"{outdir}/{sample}_raw_read_counts.tsv", "w")
    if blacklisting == True:
        header=['bin', 'chromosome','start','end', 'gc_content', 'known_bases', 'overlap_blacklist', 'use_bin', 'tumour_read_count', 'normal_read_count']
    else: 
        header=['bin', 'chromosome','start','end', 'gc_content', 'known_bases', 'use_bin', 'tumour_read_count', 'normal_read_count']
    outfile.write('\t'.join(header)+'\n')
    for r in countData:
        Line = '\t'.join(r) + '\n'
        outfile.write(Line)
    outfile.close()

    
    
    # outfile2 = open(f"{outdir}/{sample}_read_counts_filtered.tsv", "w")
    # for r in filtered_counts:
    #     Line = '\t'.join(r) + '\n'
    #     outfile2.write(Line)
    # outfile2.close()

    log2_ratio_readcounts_path = f"{outdir}/{sample}_read_counts_{nmode}_log2r.tsv"
    outfile3 = open(log2_ratio_readcounts_path, "w")
    if blacklisting == True:
        header=['bin', 'chromosome','start','end', 'gc_content', 'known_bases', 'overlap_blacklist', 'use_bin', 'log2r_copynumber']
    else: 
        header=['bin', 'chromosome','start','end', 'gc_content', 'known_bases', 'use_bin', 'log2r_copynumber']
    outfile3.write('\t'.join(header)+'\n')
    for r in normalised_counts:
        Line = '\t'.join(r) + '\n'
        outfile3.write(Line)
    outfile3.close()

    return log2_ratio_readcounts_path

if __name__ == "__main__":
    print("Read bin counter")