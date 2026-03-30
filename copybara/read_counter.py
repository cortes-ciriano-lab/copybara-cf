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
from scipy.stats import pearsonr
import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.interpolate import interp1d

def count_reads_in_curr_bin(bam, chrom, start, end, readcount_mapq, size_select, min_read_size, max_read_size):
    chunk_read_count = 0
    for read in bam.fetch(chrom, start, end, multiple_iterators = True):
        if read.is_secondary or read.is_supplementary or read.mapping_quality < readcount_mapq:
            continue
        else: 
            read_start = read.reference_start
            read_end = read.reference_end
            # read_len = read_end - read_start + 1
            read_len=len(read.query_sequence)
            if size_select == True:
                c=0
                for id,x in enumerate(min_read_size):
                    if read_len >= x and read_len <= max_read_size[id]:
                        c+=1
                if c > 0:
                    print(read_len)
                    if read_start >= start and read_start <= end:
                        chunk_read_count += 1
                    elif read_end >= start and read_end <= end:
                        chunk_read_count += 1
            else:
                if read_start >= start and read_start <= end:
                        chunk_read_count += 1
                elif read_end >= start and read_end <= end:
                        chunk_read_count += 1
    return chunk_read_count

def binned_read_counting(curr_chunk, bed_chunk, alns, nmode, blacklisting, bl_threshold, bases_filter, bases_threshold, readcount_mapq, size_select, min_read_size, max_read_size):
    print(f"    Read counting {curr_chunk} ...")
    # Open the BAM file using pysam
    bam_T = pysam.AlignmentFile(alns['tumour'], "rb")
    bam_N = pysam.AlignmentFile(alns['normal'], "rb") if nmode == "mnorm" and len(alns) == 2 else None
    # process PoN data into dictionary
    if nmode == "pon" and len(alns) == 2:
        pon_dict = {}
        with open(alns['pon'], "r") as file:
            next(file)
            for line in file:
                fields = line.strip().split("\t")
                pon_dict[fields[0]] = fields[1]
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
        # 
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
        # counting reads in tumour (and if provided) matched normal bams in each bin or adding PoN read counts 
        chunk_read_count_T = count_reads_in_curr_bin(bam_T, chrom, start, end, readcount_mapq, size_select, min_read_size, max_read_size)
        if bam_N:
            chunk_read_count_N = count_reads_in_curr_bin(bam_N, chrom, start, end, readcount_mapq, False, min_read_size, max_read_size)
        elif nmode == "pon":
            chunk_read_count_N = pon_dict[bin_name]
        else:
            chunk_read_count_N = None
        # chunk_read_count_N = count_reads_in_curr_bin(bam_N, chrom, start, end, readcount_mapq) if bam_N else None
        if blacklisting == True:
            chr_read_counts.append([bin_name, chrom, str(start), str(end), str(gc), str(bases), str(blacklist), str(use), str(chunk_read_count_T), str(chunk_read_count_N)])
        else:
            chr_read_counts.append([bin_name, chrom, str(start), str(end), str(gc), str(bases), str(use), str(chunk_read_count_T), str(chunk_read_count_N)])
        # # replace None values (lack of matched normal) with PoN count values if provided. PoN will need to be generated seperately, using the same bin size.
        # if nmode == "pon": 
        #     print('         PoN read counts being added for normalisation')
        #     if len(PoN) != len(chr_read_counts):
        #         print('ERROR: PoN will have to be generated using the same bin size as used for read counting. See documentation on how to generate PoN.')
        #     if len(PoN) == len(chr_read_counts):
        #         for id,r in enumerate(chr_read_counts):
        #             r[-1] = PoN[id][-1]
    # Close BAM file and return out
    bam_T.close()
    if bam_N is not None:
        bam_N.close()
    return(chr_read_counts)

def gc_loess(read_gc, read_outlier, gc_outlier, sample_size):
    read_counts = [float(x[0]) for x in read_gc]
    gc_content = [float(x[1]) for x in read_gc]
    # define read and gc range to remove outliers for loess fitting 
    read_range = np.quantile(read_counts, [0, 1 - read_outlier])
    gc_range = np.quantile(gc_content, [gc_outlier, 1 - gc_outlier])
    # subset data based on ranges
    counts_subset = [x for x in read_gc if float(x[0]) > read_range[0] and float(x[0]) <= read_range[1] and
                     float(x[1]) >= gc_range[0] and float(x[1]) <= gc_range[1]] 
    read_subset = [float(x[0]) for x in counts_subset]
    gc_subset = [float(x[1]) for x in counts_subset]
    # Estimate correlation between GC content and read count to determine whether or not to perform GC correction
    pear_res = pearsonr(read_subset, gc_subset)
    pear_R, pear_pval = pear_res[0],pear_res[1]
    print(f"gc bias: pearsons's R={round(pear_R,3)}; p_val={round(pear_pval,6)}")
    if (abs(pear_R)>=0.1 and pear_pval < 0.05):
        print("         GC bias detected. LOESS correction being performed.")
        # Rough LOESS fit
        if len(read_subset) >= sample_size:
            ind = np.where(read_subset)[0]
            random_ind = np.random.choice(ind, min(len(ind), sample_size), replace=False)
            rough_fit = lowess(np.asarray(read_subset)[random_ind], np.asarray(gc_subset)[random_ind], frac=0.03)
        else:
            rough_fit = lowess(read_subset, gc_subset, frac=0.03)
        # Refine the fit
        gc_grid = np.linspace(0, 1, 1000)  # Interpolation points for refined LOESS
        rough_interp = interp1d(rough_fit[:, 0], rough_fit[:, 1], bounds_error=False, fill_value=np.nan)
        rough_pred = rough_interp(gc_grid)
        # Perform the final LOESS fit
        final_fit = lowess(rough_pred, gc_grid, frac=0.3)
        # Interpolate the final fit
        # final_interp = interp1d(final_fit[:, 0], final_fit[:, 1], bounds_error=False, fill_value=np.nan)
        final_interp = interp1d(final_fit[:, 0], final_fit[:, 1], bounds_error=False, fill_value="extrapolate")
        final_pred = final_interp(gc_content)  # Predicted values for all bins
        # Correct read counts
        cor_gc = read_counts / final_pred # GC correct read counts
        ind = np.where(np.isfinite(cor_gc))[0]
        pear_res = pearsonr(cor_gc[ind], np.asarray(gc_content)[ind])
        pear_R, pear_pval = pear_res[0],pear_res[1]
        print(f"gc bias after correction: pearsons's R={round(pear_R,3)}; p_val={round(pear_pval,6)}")
        cor_gc = list(cor_gc)
    else:
        cor_gc = read_counts
    return cor_gc

def gc_correct_counts(filtered_counts, nmode, read_outlier = 0.01, gc_outlier = 0.001, sample_size = 50000):
    # extract read counts and gc content from filtered counts
    # tumour reads
    read_gc = [[float(x[-2]),float(x[4])] for x in filtered_counts]
    cor_gc = gc_loess(read_gc, read_outlier, gc_outlier, sample_size)
    # Prepare output
    gc_cor_counts = copy.deepcopy(filtered_counts)
    for id,r in enumerate(gc_cor_counts):
        r[-2] = float(cor_gc[id])
    # repeat same for normal read counts if nmode = "mnorm" or "pon"
    if nmode == "mnorm" or nmode == "pon":
        print(f'gc correction {nmode}...')
        read_gc = [[float(x[-1]),float(x[4])] for x in filtered_counts]
        cor_gc = gc_loess(read_gc, read_outlier, gc_outlier, sample_size)
        for id,r in enumerate(gc_cor_counts):
            r[-1] = float(cor_gc[id])
    return gc_cor_counts

def filter_correct_normalise(nmode, countData):
    """ perform filtering, normalisation and transformation """
    ## filtering: Remove bins with 0 reads (no sequencing data in this region). - need to be excluded for segmentation and log2 transformation
    if nmode == "self":
        filtered_counts = [x for x in countData if x[-3] == 'True' and float(x[-2]) != 0]
        gc_cor_counts = gc_correct_counts(filtered_counts, nmode) # gc correct raw read counts
        gc_cor_counts = [x for x in gc_cor_counts if float(x[-2]) > 0] # remove regions > 0 after gc correction
        med_self = statistics.median([float(x[-2]) for x in gc_cor_counts]) #estimate genome wide median for selfnormalisation
        # Normalise and log2 transform
        normalised_counts = copy.deepcopy(gc_cor_counts)
        for r in normalised_counts:
            r[-2] = str(math.log2(float(r[-2])/med_self)) # median normalise readcounts
            del r[-1]
    elif nmode == "mnorm" or nmode == "pon":
        filtered_counts = [x for x in countData if x[-3] == 'True' and float(x[-2]) != 0 and float(x[-1]) != 0]
        gc_cor_counts = gc_correct_counts(filtered_counts, nmode) # gc correct raw read counts
        gc_cor_counts = [x for x in gc_cor_counts if float(x[-2]) > 0] # remove regions < 0 after gc correction
        # cov_scaler = statistics.median([math.log2(int(x[-2])/int(x[-1])) for x in filtered_counts])
        cov_scaler = math.log2(statistics.median([(float(x[-2])/float(x[-1])) for x in gc_cor_counts]))
        normalised_counts = copy.deepcopy(gc_cor_counts)
        for r in normalised_counts:
            n = str(math.log2(float(r[-2])/float(r[-1])) - cov_scaler)
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

def estimate_coverage(tumour, ref, readcount_mapq, size_select, min_read_size, max_read_size):
    bam_T = pysam.AlignmentFile(tumour, "rb")
    fasta = pysam.FastaFile(ref)
    # Calculate Total Mapped Bases
    total_mapped_bases = 0
    for read in bam_T.fetch():
        if not read.is_unmapped and not read.is_secondary and not read.is_supplementary and read.mapping_quality >= readcount_mapq:
            if size_select == True:
                read_len=len(read.query_sequence)
                c=0
                for id,x in enumerate(min_read_size):
                    if read_len >= x and read_len <= max_read_size[id]:
                        c+=1
                if c > 0:
                    print(read_len)
                    total_mapped_bases += read.query_length 
            else:
                total_mapped_bases += read.query_length 
    # Calculate Genome Size
    genome_size = sum(fasta.get_reference_length(contig) for contig in fasta.references)
    # Compute Whole Genome Coverage
    coverage = total_mapped_bases / genome_size
    return coverage

def count_reads(outdir, tumour, normal, panel_of_normals, sample, bin_annotations_path, ref, readcount_mapq, blacklisting, bl_threshold, bases_filter, bases_threshold, size_select, min_read_size, max_read_size, threads):
    """ Perform binned read counting on bam file/files per chromosome """

    # estimate tumour bam coverage as needed later
    coverage = estimate_coverage(tumour, ref, readcount_mapq, size_select, min_read_size, max_read_size)
    print(f"Genome wide coverage of tumour bam: {coverage}")

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
            counts_binned_chr = binned_read_counting(curr_chunk, bed_chunk, aln_files, nmode, blacklisting, bl_threshold, bases_filter, bases_threshold, readcount_mapq, size_select, min_read_size, max_read_size)
            countData.append(counts_binned_chr)
        countData = [x for xs in countData for x in xs]

    else:
        print(f"multithreading using {threads} threads.")
        args_in = [[str(f"chunk {idx+1}"),bed_chunk,aln_files, nmode, blacklisting, bl_threshold, bases_filter, bases_threshold, readcount_mapq, size_select, min_read_size, max_read_size] for idx,bed_chunk in enumerate(chunks)]
        # print(args_in)
        with Pool(processes=threads) as pool:
            countData = [x for xs in list(pool.starmap(binned_read_counting, args_in)) for x in xs]

    normalised_counts = filter_correct_normalise(nmode=nmode, countData=countData)

    # Get results and write out
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

    return log2_ratio_readcounts_path,nmode,coverage

if __name__ == "__main__":
    print("Read bin counter")