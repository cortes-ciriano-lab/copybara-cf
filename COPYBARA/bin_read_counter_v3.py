"""
Script to count reads copy number estimation
Created: 12/09/2023
Python 3.9.7
Carolin Sauer
"""

from multiprocessing import Pool
from multiprocessing import cpu_count
import pysam
import argparse
import os
import timeit
import pybedtools
import copy
import statistics
import math


#----
# 1. Parse command line arguments
#----

parser = argparse.ArgumentParser(description="Count reads in given bam files across bins (using bin annotation bed file).")
parser.add_argument('-b', '--bam', type=str, help='Path to bam file.', required=True)

normalisation_group = parser.add_mutually_exclusive_group()
normalisation_group.add_argument('-nb', '--normal_bam', type=str, help='Path to matched normal bam file.', required=False)
normalisation_group.add_argument('-pon', '--panel_of_normals', type=str, help='Path to panel-of-normals (PoN) file.', required=False)

parser.add_argument('-s', '--sample_prefix', default='bam', type=str, help='Sample name of prefix to use for out files.', required=False)
parser.add_argument('-a', '--bin_annotations', type=str, help='Path to bed file with bin annotations (need to be generated prior to this). ', required=True)
parser.add_argument('-q', '--mapping_quality', type=int,  default=5, help='Mapping quality threshold used for read counting', required=False)
parser.add_argument('--no_blacklist', dest='blacklisting', action='store_false')
parser.set_defaults(blacklisting=True)
parser.add_argument('-blt', '--bl_threshold', type=int,  default='5', help='Percentage overlap between bin and blacklist threshold to tolerate for read counting (default = 0, i.e. no overlap tolerated). Please specify percentage threshold as integer, e.g. "-t 5" ', required=False)
parser.add_argument('--no_basesfilter', dest='bases_filter', action='store_false')
parser.set_defaults(bases_filter=True)
parser.add_argument('-bt', '--bases_threshold', type=int,  default='75', help='Percentage of known bases per bin required for read counting (default = 0, i.e. no filtering). Please specify percentage threshold as integer, e.g. "-bt 95" ', required=False)
parser.add_argument('-t', '--threads', type=int,  default=24, help='number of threads to be used for multiprocessing of chromosomes. Use threads = 1 to avoid multiprocessing.', required=False)
parser.add_argument('-o', '--out_dir', type=str, default='.', help='Path to out directory. Default is working directory', required=False)
args = parser.parse_args()

bam_file_path = args.bam
nbam_file_path = args.normal_bam
pon_file_path = args.panel_of_normals
prefix = args.sample_prefix
bed_file_path = args.bin_annotations
mapq = args.mapping_quality
blacklisting = args.blacklisting
bl_threshold = args.bl_threshold
bases_filter = args.bases_filter
bases_threshold = args.bases_threshold
# threads = args.threads
outdir = args.out_dir

if outdir != '.':
    try:
        os.mkdir(outdir)
    except FileExistsError:
        pass
else:
    pass

# check and define threads
threads = min(args.threads, cpu_count())
print(f"... Bin read counter will use threads = {threads}. (threads = {args.threads} defined; threads = {cpu_count()} available) ...")

if args.normal_bam is not None:
    nmode = "mnorm" #matched normal will be used for logR normalisation
    aln_files = {
			'tumour': bam_file_path,
			'normal': nbam_file_path
		}
elif args.panel_of_normals is not None:
    nmode = "pon" #panel of normals will be used for logR normalisation
    aln_files = {
			'tumour': bam_file_path
		}
elif args.normal_bam is None and args.panel_of_normals is None:
    nmode = "self"
    aln_files = {
			'tumour': bam_file_path
		}

# print(f"Normalisation mode: {nmode}")

#----
# 2. Define functions
#----
def count_reads_in_curr_bin(bam, chrom, start, end):
    chunk_read_count = 0
    for read in bam.fetch(chrom, start, end, multiple_iterators = True):
        if read.is_secondary or read.is_supplementary or read.mapping_quality < mapq:
            continue
        else: 
            read_start = read.reference_start
            read_end = read.reference_end
            if read_start >= start and read_start <= end:
                    chunk_read_count += 1
            elif read_end >= start and read_end <= end:
                    chunk_read_count += 1
    return chunk_read_count


def binned_read_counting(curr_chunk,bed_chunk,alns):
    print(f"    Read counting {curr_chunk} ...")

    # bam_T = alns['tumour']
    # bam_N = alns['normal'] if nmode == "mnorm" and len(alns) == 2 else None

    # Open the BAM file using pysam
    bam_T = pysam.AlignmentFile(alns['tumour'], "rb")
    bam_N = pysam.AlignmentFile(alns['normal'], "rb") if nmode == "mnorm" and len(alns) == 2 else None
    
    # Read the BED file to get regions and start/end positions for each chromosome
    # bed_file =  pybedtools.BedTool(bed).filter(lambda b: b.chrom == chr)

    chr_read_counts = []
    for bin in bed_chunk:
        chrom, start, end, bases = bin[0], int(bin[1]), int(bin[2]), float(bin[3])
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
            blacklist = float(bin[4])
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
        chunk_read_count_T = count_reads_in_curr_bin(bam_T, chrom, start, end)
        chunk_read_count_N = count_reads_in_curr_bin(bam_N, chrom, start, end) if bam_N else None

        # if nmode == "mnorm" and len(aln_files) == 2:
        #     chr_read_counts.append([bin_name, chrom, str(start), str(end), str(bases), str(use), str(chunk_read_count_T), str(chunk_read_count_N)])
        # elif nmode != "mnorm" and len(aln_files) == 1:
        #     chr_read_counts.append([bin_name, chrom, str(start), str(end), str(bases), str(use), str(chunk_read_count_T)])
        chr_read_counts.append([bin_name, chrom, str(start), str(end), str(bases), str(use), str(chunk_read_count_T), str(chunk_read_count_N)])

    # Close BAM file and return out
    bam_T.close()
    if bam_N is not None:
        bam_N.close()    
    return(chr_read_counts)


def filter_and_normalise(nmode, countData):
    # perform filtering, normalisation and transformation
    ## filtering: Remove bins with 0 reads (no sequencing data in this region). - need to be excluded for segmentation and log2 transformation
    if nmode == "self":
        filtered_counts = [x for x in countData if x[-3] == 'True' and int(x[-2]) != 0]
        med_self = statistics.median([int(x[-2]) for x in filtered_counts]) #estimate genome wide median for selfnormalisation
        # Normalise and log2 transform
        normalised_counts = copy.deepcopy(filtered_counts)
        for r in normalised_counts:
            r[-2] = str(math.log2(int(r[-2])/med_self)) # median normalise readcounts
            del r[-1]
    elif nmode == "pon":
        print("This function is yet to be implemented...")
    elif nmode == "mnorm":
        filtered_counts = [x for x in countData if x[-3] == 'True' and int(x[-2]) != 0 and int(x[-1]) != 0]
        cov_scaler = statistics.median([math.log2(int(x[-2])/int(x[-1])) for x in filtered_counts])
        normalised_counts = copy.deepcopy(filtered_counts)
        for r in normalised_counts:
            n = str(math.log2(int(r[-2])/int(r[-1])) - cov_scaler)
            del r[-2:]
            r.append(n)
    # return out        
    return filtered_counts, normalised_counts


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


#----    
# 3. Perform binned read counting on bam file/files per chromosome
#----
def main():

    print(f"Normalisation mode: {nmode}")

    # Define contig names from annotation bed file
    # chr_names = list(dict.fromkeys([x[0] for x in pybedtools.BedTool(bed_file_path)]))

    # Define bed_file chunks
    bed_file =  pybedtools.BedTool(bed_file_path)
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
            counts_binned_chr = binned_read_counting(curr_chunk,bed_chunk, aln_files)
            countData.append(counts_binned_chr)
        countData = [x for xs in countData for x in xs]

    else:
        print(f"multithreading using {threads} threads.")
        args_in = [[str(f"chunk {idx+1}"),bed_chunk,aln_files] for idx,bed_chunk in enumerate(chunks)]
        # print(args_in)
        with Pool(processes=threads) as pool:
            countData = [x for xs in list(pool.starmap(binned_read_counting, args_in)) for x in xs]

    filtered_counts, normalised_counts = filter_and_normalise(nmode=nmode, countData=countData)

    #----
    # 4. Get results and write out
    #----
    outfile = open(f"{outdir}/{prefix}_read_counts.tsv", "w")
    for r in countData:
        Line = '\t'.join(r) + '\n'
        outfile.write(Line)
    outfile.close()           

    outfile2 = open(f"{outdir}/{prefix}_read_counts_filtered.tsv", "w")
    for r in filtered_counts:
        Line = '\t'.join(r) + '\n'
        outfile2.write(Line)
    outfile2.close()    

    outfile3 = open(f"{outdir}/{prefix}_read_counts_{nmode}_log2r.tsv", "w")
    for r in normalised_counts:
        Line = '\t'.join(r) + '\n'
        outfile3.write(Line)
    outfile3.close()    


if __name__ == "__main__":
    start_t = timeit.default_timer()
    main()
    stop = timeit.default_timer()
    Seconds = round(stop - start_t)
    print(f"Computation time: {Seconds} seconds\n") 

