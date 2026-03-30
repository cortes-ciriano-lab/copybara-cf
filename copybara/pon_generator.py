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
import statistics
import math


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


    # 4. Get results and write out
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