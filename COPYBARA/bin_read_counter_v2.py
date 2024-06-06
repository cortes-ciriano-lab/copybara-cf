from multiprocessing import Pool
import pysam
import argparse
import os
import timeit
import statistics
import math


start_t = timeit.default_timer()

##############
# TO DO
##############

# STEP 1 - multiprocessing by chromosome
# STEP 2 - normalisation [self,pon,norm] 

# THEN... smoothening | segmentation | tMAD | rascal fitting?? ploidy/purity | --> test on real samples and longitudinal!
# Also add gene mapper? 
##############

#----
# 1. Parse command line arguments
#----
parser = argparse.ArgumentParser(description="Count reads in given bam files across bins (using bin annotation bed file).")
parser.add_argument('-b', '--bam', type=str, help='Path to bam file.', required=True)
# parser.add_argument('-nb', '--normal_bam', type=str, help='Path to matched normal bam file.', required=False)
# parser.add_argument('-np', '--normal_panel', type=str, help='Path to panel-of-normals (PoN) file.', required=False)
# parser.add_argument('-nb', '--normal_bam', type=str, help='Method of normalisation (self, pon, or norm; default = self).', required=False)
parser.add_argument('-s', '--sample_prefix', default='bam', type=str, help='Sample name of prefix to use for out files.', required=False)
parser.add_argument('-a', '--bin_annotations', type=str, help='Path to bed file with bin annotations (need to be generated prior to this). ', required=True)
parser.add_argument('-q', '--mapping_quality', type=int,  default=20, help='Mapping quality threshold used for read counting', required=False)
parser.add_argument('--no_blacklist', dest='blacklisting', action='store_false')
parser.set_defaults(blacklisting=True)
parser.add_argument('-t', '--bl_threshold', type=int,  default='0', help='Percentage overlap between bin and blacklist threshold to tolerate for read counting (default = 0, i.e. no overlap tolerated). Please specify percentage threshold as integer, e.g. "-t 5" ', required=False)
parser.add_argument('--no_basesfilter', dest='bases_filter', action='store_false')
parser.set_defaults(bases_filter=True)
parser.add_argument('-bt', '--bases_threshold', type=int,  default='95', help='Percentage of known bases per bin required for read counting (default = 0, i.e. no filtering). Please specify percentage threshold as integer, e.g. "-bt 95" ', required=False)
# parser.add_argument('-sm', '--smoothing_level', type=int,  default='3', help='Number of bases (window size) used to smooth outlier bins (default = 3). Set to "-sm 0" to skip smoothening.', required=False)
parser.add_argument('-o', '--out_dir', type=str, default='.', help='Path to out directory. Default is working directory', required=False)
args = parser.parse_args()

bam_file_path = args.bam
prefix = args.sample_prefix
bed_file_path = args.bin_annotations
mapq = args.mapping_quality
blacklisting = args.blacklisting
bl_threshold = args.bl_threshold
bases_filter = args.bases_filter
bases_threshold = args.bases_threshold
smoothing_level = args.smoothing_level
outdir = args.out_dir

if outdir != '.':
    try:
        os.mkdir(outdir)
    except FileExistsError:
        pass
else:
    pass

# Open the BAM file using pysam
bam_file = pysam.AlignmentFile(bam_file_path, "rb")

outfile = open(f"{outdir}/{prefix}_read_counts.tsv", "w")

# Read the BED file to get regions and start/end positions 
with open(bed_file_path, "r") as bed_file:
    
    List_of_read_counts = []
    chunk_read_incrementer = {}
    
    for line in bed_file:
        fields = line.strip().split("\t")
        chrom, start, end, bases = fields[0], int(fields[1]), int(fields[2]), float(fields[3])
        bin_name = f"{chrom}:{start}_{end}"
        
        # Read Counting - Iterate through the BAM file and count reads in each bin
        chunk_read_count = 0
        for read in bam_file.fetch(chrom, start, end):
            if read.is_secondary or read.is_supplementary or read.mapping_quality < mapq:
                continue
            else: 
                read_start = read.reference_start
                read_end = read.reference_end
            
                if read_start >= start and read_start <= end:
                        chunk_read_count += 1
                elif read_end >= start and read_end <= end:
                        chunk_read_count += 1
            
        # Filtering of bins by assigning each bin to boolean 'use' variable
        if blacklisting == False:
            if bases_filter == False:
                use = True
            elif bases_filter == True:
                if bases >= bases_threshold:
                    use = True
                else:
                    use = False
            List_of_read_counts.append([bin_name, chrom, str(start), str(end), str(bases), str(use), str(chunk_read_count)])
        elif blacklisting == True:
            blacklist = float(fields[4])
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
            List_of_read_counts.append([bin_name, chrom, str(start), str(end), str(bases), str(blacklist), str(use), str(chunk_read_count)])
# Close the BAM file
bam_file.close()

# write output
for r in List_of_read_counts:
    Line = '\t'.join(r) + '\n'
    outfile.write(Line)        
outfile.close()            

# Remove bins with 0 reads (no sequencing data in this region). - need to be excluded for segmentation and log2 transformation
# If blacklisting == True or bases_filter == True this will also filter the results to exclude blacklisted regions.
outfile2 = open(f"{outdir}/{prefix}_read_counts_filtered.tsv", "w")
filtered_counts = [x for x in List_of_read_counts if x[-2] == 'True' and int(x[-1]) != 0]
for r in filtered_counts:
    Line = '\t'.join(r) + '\n'
    outfile2.write(Line)        
outfile2.close()

# Normalise reads by median
med = statistics.median([int(x[-1]) for x in filtered_counts])
std = statistics.stdev([int(x[-1]) for x in filtered_counts])

outfile3 = open(f"{outdir}/{prefix}_read_counts_mednorm.tsv", "w")
normalised_counts = filtered_counts
for r in normalised_counts:
    r[-1] = str(int(r[-1])/med) # median normalise readcounts
    # r[-1] = str(math.sqrt((int(r[-1])/med) + 3/8)) #anscombe sqrt transform
    # might want to change this to normalis to PoN or make that an option...
    Line = '\t'.join(r) + '\n'
    outfile3.write(Line)
outfile3.close()

# Log2 transform data (prior to smoothening) 
outfile4 = open(f"{outdir}/{prefix}_read_counts_log2r.tsv", "w")
log2r_counts = normalised_counts
for r in log2r_counts:
    r[-1] = str(math.log2(float(r[-1])))
    Line = '\t'.join(r) + '\n'
    outfile4.write(Line)
outfile4.close()


# Log transform
# outfile4 = open(f"{outdir}/{prefix}_read_counts_normalised_log.tsv", "w")
# log_counts = normalised_counts
# for r in log_counts:
#     r[-1] = str(math.log(float(r[-1])))
#     # might want to change this to normalis to PoN or make that an option...
#     Line = '\t'.join(r) + '\n'
#     outfile4.write(Line)
# outfile4.close()

# Log data if no smoothing is applied prior to CBS

# Smoothen normalised bins
# def smooth_outliers(normalised_data, smoothing_level=3):
#     if smoothing_level <= 0:
#         return normalised_data  # Skip smoothing and use original (normalised) data for CBS

#     smoothed_counts = []
#     half_window = smoothing_level // 2

#     for val in range(len(normalised_data)):
#         start = max(0, val - half_window)
#         end = min(len(normalised_data), val + half_window + 1)
#         window = normalised_data[start:end]
#         smoothed_val = sum(window) / len(window)
#         smoothed_counts.append(smoothed_val)
        
#     return smoothed_counts

# data_to_smoothen = [float(x[-1]) for x in normalised_counts]
# smoothed_counts = smooth_outliers(data_to_smoothen, smoothing_level)

# outfile4 = open(f"{outdir}/{prefix}_binned_copy_number.tsv", "w")
# binned_copy_number = normalised_counts
# for i in range(len(binned_copy_number)):
#     binned_copy_number[i][-1] = str(smoothed_counts[i])
#     # print(binned_copy_number[i][-1], smoothed_counts[i])
#     Line = '\t'.join(binned_copy_number[i]) + '\n'
#     outfile4.write(Line)
# outfile4.close()

stop = timeit.default_timer()
Seconds = round(stop - start_t)
print(f"Computation time: {Seconds} seconds\n") 
