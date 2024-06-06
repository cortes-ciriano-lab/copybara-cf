from multiprocessing import Pool
import argparse
import os
import timeit
import pybedtools
import pysam

start_t = timeit.default_timer()

#----
# 1. Parse command line arguments
#----
parser = argparse.ArgumentParser(description="Generate bin annotations from fasta file with given bin size.")
#parser.add_argument('-w', '--bin_size', type=int, help='Bin window size in kbp', required=True)
#parser.add_argument('-f', '--fasta_file', type=str, help='Path to the reference FASTA file', required=True)
#parser.add_argument('-b', '--black_list', type=str, help='Path to the blacklist file', required=False)
parser.add_argument('-c', '--chromosomes', nargs='+', default='all', help='Chromosomes to analyse. To run on all chromosomes, leave unspecified (default). To run on a subset of chromosomes only, specify the chromosome numbers separated by spaces. For x and y chromosomes, use 23 and 24, respectively.  E.g. use "-c 1 4 23 24" to run chromosomes 1, 4, X and Y', required=False)
#parser.add_argument('-t', '--threads', type=int,  default=24, help='number of threads to be used for multiprocessing of chromosomes', required=False)
#parser.add_argument('-o', '--out_dir', type=str, default='.', help='Path to out directory. Default is working directory', required=False)
args = parser.parse_args()

#bin_size = args.bin_size * int(1000)
#fasta_file_path = args.fasta_file
#blacklist_file = args.black_list
contigs = args.chromosomes
#threads = args.threads
#outdir = args.out_dir

print(type(contigs))
print(contigs)

if len(contigs) >= 1:
    print("len contigs", len(contigs))

exit()


if outdir != '.':
    try:
        os.mkdir(outdir)
    except FileExistsError:
        pass
else:
    pass

#----
# 2. Define functions
#----
# Function to estimate overlaps 
def overlaps(a, b):
    return min(a[1], b[1]) - max(a[0], b[0]) + 1

# Function to process each chromosome of interest to run in parallel
def process_chromosome(chrom, chr_length, fasta_file_path):
# def process_chromosome(chr_in):
    # chrom, chr_length, fasta_file_path = chr_in
    print(f"    Processing {chrom} ...")

    fasta_file = pysam.FastaFile(fasta_file_path)
    chr_starts = list(range(1, chr_length + 1, bin_size))
    chr_bins = [[start,min((start + bin_size - 1),chr_length)] for start in chr_starts]
    
    if (blacklist_file != None):
        blacklist_chr = pybedtools.BedTool(blacklist_file).filter(lambda b: b.chrom == chrom)    
        blacklist_chr_list = [[int(row[1]), int(row[2])] for row in blacklist_chr]
    else: 
        pass

    # Define temp out file
    # temp_outfile_name = f"temp/{chrom}.{int(bin_size/1000)}kbp.tmp.tsv"
    # temp_outfile = open(temp_outfile_name, "w")
    List_of_bins = []
    for bin in chr_bins:
        # print(bin)
        start, end = bin 
        
        # Estimate overlap with blacklist for each bin
        # Using overlap function
        if blacklist_file != None:
            overlap_list = [overlaps(bin,row) for row in blacklist_chr_list]
            overlap_list = [x for x in overlap_list if x > 0]
            overlap_pct = sum(overlap_list) / float(bin_size) * 100
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
        # overlap_pct = bl_bin_bases / float(bin_size) * 100

        # Estimate number/fraction of known bases for each bin
        chr_bin_seq = fasta_file.fetch(chrom, start, end)
        acgt_len = 0
        for base in chr_bin_seq:
            if base.upper() in 'ACTG':
                acgt_len += 1
        chr_bases = acgt_len / float(bin_size) * 100

        # Line = f"{chrom}\t{start}\t{end}\t{chr_bases}" 
        # Line = Line + "\n" 
        # temp_outfile.write(Line)

        if blacklist_file != None:
            List_of_bins.append([chrom, str(start), str(end), str(chr_bases), str(overlap_pct)])
        else:
            List_of_bins.append([chrom, str(start), str(end), str(chr_bases)])

    # temp_outfile.close()
    fasta_file.close()
    # return(temp_outfile_name)
    return(List_of_bins)

#----    
# 3. Extract chormosome info and pass fasta file into chr_in list
#----
fasta = pysam.FastaFile(fasta_file_path)
if contigs != 'all':
    chr_names = [fasta.references[(int(x)-1)] for x in contigs]     
else:
    chr_names = fasta.references[0:24]
chr_in = [[chr,fasta.get_reference_length(chr),fasta_file_path] for chr in chr_names]
fasta.close()

#----
# 4. Generate bins and count bases per bins for each chrom in parallel (run function using multithreader)
#----

# only use multiprocessing if more than 1 threads available or being used. 

# if threads == 1:
    


# with ThreadPoolExecutor(max_workers=24) as executor:
#     chrData = list(executor.map(process_chromosome, chr_in))
with Pool(processes=threads) as pool:
    chrData = list(pool.starmap(process_chromosome, chr_in))

#----
# 5. Concat results into a single bed file
#----
if contigs != 'all':
    outfile_name = f"{outdir}/{int(bin_size/1000)}kbp_bin_ref_subset.bed"
else:
    outfile_name = f"{outdir}/{int(bin_size/1000)}kbp_bin_ref_all.bed"

outfile = open(outfile_name, "w")

for obj in chrData:
    for r in obj:
        Line = '\t'.join(r) + '\n'
        outfile.write(Line)
        
outfile.close()

#--------
stop = timeit.default_timer()
Seconds = round(stop - start_t)
print(f"Computation time (bin_ref_generator): {Seconds} seconds\n") 
