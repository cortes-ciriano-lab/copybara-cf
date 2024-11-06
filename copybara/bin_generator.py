"""
Script to generate binned reference for copy number estimation
Created: 31/08/2023
Python 3.9.7
Carolin Sauer
"""

from multiprocessing import Pool
from multiprocessing import cpu_count
import pybedtools
import pysam
import gzip

#----
# Define functions
#----
def overlaps(a, b):
    '''
    Function to estimate overlaps between two regions
    '''
    return min(a[1], b[1]) - max(a[0], b[0]) + 1

def process_chromosome(chrom, chr_length, fasta_file_path, bin_size, blacklist):
    '''
    Function to process each chromosome of interest (either in a loop or in parallel) to generate bins based on bin_size with annotated blacklist and unknown bases values. 
    '''
    print(f"    Processing {chrom} ...")

    fasta_file = pysam.FastaFile(fasta_file_path)
    chr_starts = list(range(1, chr_length + 1, bin_size))

    # define bins
    chr_bins = []
    for i in range(len(chr_starts)-1):
        bin_start,bin_end = chr_starts[i], chr_starts[i+1]-1
        chr_bins.append([bin_start,bin_end])
    # add last bin
    chr_bins.append([chr_starts[-1],chr_length])

    # estimate overlap with blacklist if blacklist is provided
    if (blacklist is not None):
        blacklist_chr = pybedtools.BedTool(blacklist).filter(lambda b: b.chrom == chrom)
        blacklist_chr_list = [[int(row[1]), int(row[2])] for row in blacklist_chr]
    else:
        pass

    # Generate output
    List_of_bins = []
    for bin in chr_bins:
        start, end = bin
        cur_bin_size = float(end)-float(start)+1

        # Estimate overlap with blacklist for each bin
        # Using overlap function
        if blacklist is not None:
            overlap_list = [overlaps(bin,row) for row in blacklist_chr_list]
            overlap_list = [x for x in overlap_list if x > 0]
            overlap_pct = sum(overlap_list) / cur_bin_size * 100
            # print(overlap_pct)

        # Estimate overlap manually
        # bl_bin_bases = 0
        # for row in blacklist_chr_list:
        #     bed_start, bed_end = row
        #     if bed_start >= start and bed_start <= end:
        #         bl_length = min(bed_end, end) - bed_start + 1
        #         bl_bin_bases += bl_length
        #     elif bed_end >= start and bed_end <= end:
        #         bl_length = bed_end - min(bed_start,start) + 1
        #         bl_bin_bases += bl_length
        # overlap_pct = bl_bin_bases / float(cur_bin_size) * 100

        # Estimate number/fraction of known bases for each bin
        chr_bin_seq = fasta_file.fetch(chrom, start, end)
        acgt_len = 0
        for base in chr_bin_seq:
            if base.upper() in 'ACTG':
                acgt_len += 1
        chr_bases = acgt_len / cur_bin_size * 100

        if blacklist is not None:
            List_of_bins.append([chrom, str(start), str(end), str(chr_bases), str(overlap_pct)])
        else:
            List_of_bins.append([chrom, str(start), str(end), str(chr_bases)])

    fasta_file.close()
    return(List_of_bins)

def generate_bins(outdir, sample, ref, chromosomes, bin_size, blacklist, threads):
    '''
    Main function to process and bin fasta file
    '''
    # check and define threads
    new_threads = min(threads, cpu_count())
    print(f"... Bin generator will use threads = {new_threads}. (threads = {threads} defined; threads = {cpu_count()} available) ...")
    threads = new_threads
    bin_size = bin_size * int(1000) # in kb

    # a. Extract chormosome info and pass fasta file into chr_in list
    fasta = pysam.FastaFile(ref)
    contigs = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']
    ref_contigs = fasta.references
    if 'chr1' in ref_contigs:
        contigs = [f'chr{x}' for x in contigs]
    if chromosomes != 'all':
        chr_names = [contigs[(int(x)-1)] for x in chromosomes]
    else:
        chr_names = contigs
    #   define input
    chr_in = [[chrom,fasta.get_reference_length(chrom),ref,bin_size,blacklist] for chrom in chr_names]
    fasta.close()

    # c. Generate bins and count bases per bins for each chrom in parallel (run function using multithreader)
    ## only use multiprocessing if more than 1 thread available/being used.
    if threads == 1:
        # loop through chromosomes
        print("multithreading skipped.")
        chrData = []
        for contig in chr_in:
            chrom, chr_length, fasta_file_path = contig[0], contig[1], contig[2]
            binned_chr = process_chromosome(chrom, chr_length, fasta_file_path, bin_size, blacklist)
            chrData.append(binned_chr)
    else:
        print(f"multithreading using {threads} threads.")
        with Pool(processes=threads) as pool:
            chrData = list(pool.starmap(process_chromosome, chr_in))

    # c. Concat results into a single bed file
    if chromosomes != 'all':
        outfile_name = f"{outdir}/{int(bin_size/1000)}kbp_bin_ref_subset_{sample}.bed"
    else:
        outfile_name = f"{outdir}/{int(bin_size/1000)}kbp_bin_ref_all_{sample}.bed"

    outfile = open(outfile_name, "w")

    for obj in chrData:
        for r in obj:
            Line = '\t'.join(r) + '\n'
            outfile.write(Line)

    outfile.close()

    return outfile_name

if __name__ == "__main__":
    print("Bin generator")
