from multiprocessing import Pool
from multiprocessing import cpu_count
import argparse
import timeit
import cyvcf2
import pysam
import pybedtools
import math
import os

#----
# 1. Parse command line arguments
#----
parser = argparse.ArgumentParser(description="Count alleles/bases in tumour bam using heterozygous SNPs from phased normal.")
parser.add_argument('-b', '--bam', type=str, help='Path to tumour bam file.', required=True)
parser.add_argument('-v', '--phased_vcf', type=str, help='Path to phased vcf file to extract heterozygous SNPs for allele counting.', required=True)
# parser.add_argument('-hb', '--hets_bed', type=str, help='Path to bed file containing heterozygous SNPs extracted from normal phased.vcf.', required=True)
parser.add_argument('-s', '--sample_prefix', default='bam', type=str, help='Sample name of prefix to use for out files.', required=False)
# parser.add_argument('-c', '--chromosomes', nargs='+', default='all', help='Chromosomes to analyse. To run on all chromosomes, leave unspecified (default). To run on a subset of chromosomes only, specify the chromosome numbers separated by spaces. For x and y chromosomes, use 23 and 24, respectively.  E.g. use "-c 1 4 23 24" to run chromosomes 1, 4, X and Y', required=False)
parser.add_argument('-q', '--mapping_quality', type=int,  default=0, help='Mapping quality threshold used reads to be included in the allele counting (default = 0)', required=False)
parser.add_argument('-mr', '--min_reads', type=int,  default=10, help='Minimum number of reads required per het SNP site for allele counting (default = 10)', required=False)
parser.add_argument('-t', '--threads', type=int,  default=48, help='number of threads to be used for multiprocessing of chromosomes. Use threads = 1 to avoid multiprocessing.', required=False)
parser.add_argument('-o', '--out_dir', type=str, default='.', help='Path to out directory. Default is working directory', required=False)
args = parser.parse_args()

bam_file_path = args.bam
phased_vcf_path = args.phased_vcf
prefix = args.sample_prefix
# contigs = args.chromosomes
mapq = args.mapping_quality
min_reads = args.min_reads
outdir = args.out_dir

if outdir != '.':
    try:
        os.mkdir(outdir)
    except FileExistsError:
        pass
else:
    pass

hets_bed_path = f"{outdir}/{prefix}_phased_het_snps.bed"

# check and define threads
threads = min(args.threads, cpu_count())
print(f"... Allele counter will use threads = {threads}. (threads = {args.threads} defined; threads = {cpu_count()} available) ...")

#----
# 2. Define functions
#----
def extract_hets(phased_vcf):
    '''
    Extracts heterozygous SNPs from phased.vcf (e.g. from matched normal bam).
    '''
    print(f"    Extracting phased heterozygous SNPs from {phased_vcf} ...")
    vcf_reader = cyvcf2.VCF(phased_vcf)
    hets = []
    # iterate through variants
    for variant in vcf_reader:
        # check if variant is a SNP
        if variant.is_snp:
            # iterate through each genotype
            for s_idx, gt in enumerate(variant.genotypes):
                # only get heterozygous snps
                if gt[0] != gt[1]:
                    # print(s_idx, gt)
                    # sample = vcf_reader.samples[s_idx]
                    ps = int(variant.format('PS')[s_idx]) if 'PS' in variant.FORMAT else None
                    gt_str = f"{gt[0]}|{gt[1]}" if gt[2] == 1 else f"{gt[0]}/{gt[1]}"
                    # dp = int(variant.format('DP')[s_idx]) if 'DP' in variant.FORMAT else None
                    # ad = int(variant.format('AD')[s_idx]) if 'AD' in variant.FORMAT else None
                    var_out = [variant.CHROM,str(variant.POS),str(variant.POS),variant.REF,variant.ALT[0],gt_str,str(ps)]
                    hets.append(var_out)
    return hets

def chunkify_bed(bed_file, chunk_size):
    '''
    Divides bed file into chunks based on threads available for multiprocessing.
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

def allele_counter(curr_chunk, sites, bam):
    '''
    Count alleles across heterozygous SNP positions from phased normal bam
    '''
    print(f"    Read counting {curr_chunk} ...")   
    # open tumour bam
    with pysam.AlignmentFile(bam, "rb") as bamfile:
        list_out = []
        for site in sites:
            # print(site)
            dict_out = {}
            # if len(site[0])==len(chr_now):
            #     chr = site[0];#[3];
            # else:
            #     chr = site[0]#[3:5];
            # if chr != chr_now:
            #     continue
            # chr2    61791980    61791980    G    A    1|0    58907085
            chr = site[0]
            start = int(site[1])-1
            end = int(site[2])  # start and end cannot be the same position or otherwise the fetching does not work
            ref = site[3]
            alt = site[4]
            GT = site[5]
            PS = site[6]

            reads = [
                read for read in bamfile.fetch(
                    chr,
                    start,
                    end,
                    multiple_iterators=True
                )
                ]
            # filter reads
            reads = [
                read for read in reads if 
                read.mapping_quality >= mapq and not read.is_secondary and not read.is_supplementary
                ]
            # perform base counting 
            try:
                if len(reads) >= min_reads:
                    key_now = f"{chr}_{start}_{end}_{ref}_{alt}"
                    # key_now = "{}\t{}\t{}\t{}\t{}".format(chr, start, end, ref ,alt)
                    for read in reads:
                        read_sequence = read.seq

                        aligned_pos = read.get_reference_positions(
                            full_length=True
                        )  
                        try:
                            idx = aligned_pos.index(start)
                            #DP+=1
                            snp_read = read_sequence[idx]
                            if key_now in dict_out: #.has_key(key_now):
                                dict_out[key_now][snp_read] = dict_out[key_now][snp_read] + 1
                            else:
                                dict_out[key_now]={"A":0, "C":0, "G":0, "T":0, "N":0}
                                dict_out[key_now][snp_read] = dict_out[key_now][snp_read] + 1
                        except:
                            continue
                    DP = dict_out[key_now]["A"] + dict_out[key_now]["C"] + dict_out[key_now]["G"] + dict_out[key_now]["T"] 
                    if GT == "0|1":
                        AF_0 = dict_out[key_now][ref] / DP
                        AF_1 = dict_out[key_now][alt] / DP
                    if GT == "1|0":
                        AF_0 = dict_out[key_now][alt] / DP
                        AF_1 = dict_out[key_now][ref] / DP
                    #print("{} {}".format(AF_0, AF_1))
                    dict_out[key_now]["AF_0"] = AF_0
                    dict_out[key_now]["AF_1"] = AF_1

                    out = [chr, str(start), str(end), ref, alt, str(dict_out[key_now]["A"]), str(dict_out[key_now]["C"]), str(dict_out[key_now]["G"]), str(dict_out[key_now]["T"]), str(dict_out[key_now]["N"]) , str(dict_out[key_now]["AF_0"]), str(dict_out[key_now]["AF_1"]), GT, str(PS)]

                list_out.append(out)
            except:
                continue
    return(list_out)


#----
# 3. Prepare input and run allele counter across heterozygous SNPs
#----
def main():
    # extract heterozygous SNPs
    outfile = open(hets_bed_path, "w")
    het_snps = extract_hets(phased_vcf_path)
    for r in het_snps:
        Line = '\t'.join(r) + '\n'
        outfile.write(Line)
    outfile.close()

    ## RUN IN CHUNKS TO SPEED UP
    # Define bed_file chunks
    bed_file =  pybedtools.BedTool(hets_bed_path)
    bed_length = sum(1 for _ in bed_file)
    chunk_size = bed_length // threads + (bed_length % threads > 0)
    print(f"    splitting bed file into n = {math.ceil(bed_length / chunk_size)} chunks ...")
    # Split the BED file into chunks
    chunks = chunkify_bed(bed_file, chunk_size)

    # only use multiprocessing if more than 1 thread available/being used. 
    if threads == 1:
        # loop through chromosomes
        print("multithreading skipped.")
        allele_counts = []
        for idx,bed_chunk in enumerate(chunks):
            curr_chunk = f"chunk {idx+1}"
            allele_counts_chunk = allele_counter(curr_chunk, bed_chunk, bam_file_path)
            allele_counts.append(allele_counts_chunk)
        allele_counts = [x for xs in allele_counts for x in xs]

    else:
        print(f"multithreading using {threads} threads.")
        args_in = [[str(f"chunk {idx+1}"),bed_chunk,bam_file_path] for idx,bed_chunk in enumerate(chunks)]
        # print(args_in)
        with Pool(processes=threads) as pool:
            allele_counts = [x for xs in list(pool.starmap(allele_counter, args_in)) for x in xs]

    # print(allele_counts)
    # print(len(allele_counts))

    #----
    # 4. Get results and write out
    #----
    outfile = open(f"{outdir}/{prefix}_allele_counts_hetSNPs.tsv", "w")
    for r in allele_counts:
        Line = '\t'.join(r) + '\n'
        outfile.write(Line)
    outfile.close()           


if __name__ == "__main__":
    start_t = timeit.default_timer()
    main()
    stop = timeit.default_timer()
    Seconds = round(stop - start_t)
    print(f"Computation time (allele_counts_phased_hetSNPs): {Seconds} seconds\n") 







## BY CHROMOSOME
# with open(hets_bed_path) as bed:
#         reader = csv.reader(bed, delimiter="\t")
#         sites = list(reader)
# chr_names = list(dict.fromkeys([x[0] for x in sites]))
# ### TESTING ###
# chr_names = ['chr19', 'chr20', 'chr21', 'chr22']
# chr_in = [[chrom, bam_file_path, [x for x in sites if x[0] == chrom]] for chrom in chr_names]

# ## only use multiprocessing if more than 1 thread available/being used. 
# if threads == 1:
#     # loop through chromosomes
#     print("multithreading skipped.")
#     allele_counts = []
#     for contig in chr_in:
#         chrom, bam, chr_sites = contig[0], contig[1], contig[2]
#         # print(chrom)
#         chr_res = allele_counter(chrom, bam, chr_sites)
#         allele_counts.append(chr_res)
# else:
#     print(f"multithreading using {threads} threads.")
#     with Pool(processes=threads) as pool:
#         allele_counts = list(pool.starmap(allele_counter, chr_in))

