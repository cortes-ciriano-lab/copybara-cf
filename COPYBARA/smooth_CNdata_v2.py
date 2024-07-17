"""
Script to smooth relative copy number values prior to segmentation
Created: 15/01/2024
Python 3.9.7
Carolin Sauer
"""

import argparse
import copy
from math import sqrt
from math import log2
from scipy.stats import norm
from statistics import median
import os
import timeit


#----
# 1. Parse command line arguments
#----
parser = argparse.ArgumentParser(description="smoothen copy number count data (Log2R)")
parser.add_argument('-rc', '--read_counts', type=str, help='Path to read_counts_log2r.tsv file.', required=True)
parser.add_argument('-o', '--out_dir', type=str, default='.', help='Path to out directory. Default is working directory', required=False)
parser.add_argument('-sl', '--smoothing_level', type=int, default='10', help='Size of neighbourhood for smoothing.', required=False)
parser.add_argument('-tr', '--trim', type=float, default='0.025', help='Trimming percentage to be used.', required=False)
# parser.add_argument('-s', '--sample_prefix', type=str, help='Sample name or prefix to use for out files.', required=True)
# parser.add_argument('-t', '--threads', type=int,  default=24, help='number of threads to be used for multiprocessing of chromosomes. Use threads = 1 to avoid multiprocessing.', required=False)


args = parser.parse_args()

read_counts_path = args.read_counts
outdir = args.out_dir
smoothing_level = args.smoothing_level
trim = args.trim
# prefix = args.sample_prefix
# threads = args.threads

#----
# 2. Define functions
#----

# function to calculate an influence factor (scaling factor) used to determine outliers while using a trimming percentage (value between 0-1)
def dnorm(x):
    return 1 / sqrt(2 * 3.141592653589793) * (2.718281828459045 ** (-0.5 * x ** 2))
# dnorm needed to estimate constants related to the probability density function (PDF) of the standard normal distribution
def inflfact(trim):
    # a = -1.0 * sqrt(2) * norm.ppf(trim)
    a = norm.ppf(1-trim)
    x = [-a + i * (2 * a) / 10000 for i in range(10001)]
    x1 = [(x[i] + x[i + 1]) / 2 for i in range(10000)]
    result = 1 / (sum([(x1[i] ** 2 * dnorm(x1[i])) / (1 - 2 * trim) for i in range(10000)]) * (2 * a / 10000))
    return result

# Function to calculate a variance measure used to define oSD and sSD for smoothing CN data
def trimmed_variance(genomdat, trim):
    n = len(genomdat)
    n_keep = round((1 - 2 * trim) * (n - 1))
    infl_fact = inflfact(trim)
    sorted_diffs = sorted([abs(genomdat[i] - genomdat[i + 1]) for i in range(n - 1)])
    trimmed_var = sum([sorted_diffs[i] ** 2 for i in range(n_keep)]) / (2 * n_keep) * infl_fact
    return trimmed_var

# Function to smoothen CN data (adapted from DNAcopy Fortran code)
def smoothen(chr, in_data, trim, smoothing_level):
    print(f"    Smoothening log2r read counts for chromosome {chr} ...")
    #slice out chromosome to process...
    chr_in_data = [x for x in in_data if x[1] == chr]
    if smoothing_level <= 0:
        return chr_in_data  # Skip smoothing and use original (normalised) data for CBS
    else:
        vals = []
        for x in chr_in_data:
            val = float(x[-1])
            vals.append(val)
        oSD = sqrt(trimmed_variance(vals, trim)) * 4
        # print(f"oSD = {oSD}")
        sSD = sqrt(trimmed_variance(vals, trim)) * 2
        # print(f"sSD = {sSD}")
        k = smoothing_level
        start = 0
        end = len(vals)
        # print(oSD, sSD, k, start, end)
        smoothed_counts = [0.0] * len(vals)

        for i in range(len(vals)):
            # print(f"i = {i}, {vals[i]}")
            nbh_start = max(start, i - k)
            nbh_end = min(end, i + k)
            
            above_nb = 100 * oSD #initate as high number
            below_nb = 100 * oSD #initate as high number
            
            # neighbourhood array
            nbh_i = vals[nbh_start:nbh_end+1]
            # print(nbh_i)

            for nb in range(nbh_start, nbh_end):
                if nb != i:
                    dist_inb = vals[i] - vals[nb]
                    # print(f"dist_inb = {dist_inb}")
                    if abs(dist_inb) <= oSD:
                        smoothed_counts[i] = vals[i]
                        break
                    else:
                        if dist_inb < above_nb:
                            above_nb = dist_inb
                        if -dist_inb < below_nb:
                            below_nb = -dist_inb
                    # print(f"above_nb and below_nb for {i} are {above_nb} and {below_nb}")

            #SMOOTHING
            # if all points in the neighbourhoud are > i (i.e. i is a true outlier), then above_nb < 0 and below_nb > oSD; 
            # and if all points in neighbourhoud are < i (i.e. i is a true outlier). then above_nb > oSD and below_nb < 0
            
            # However, if both above_nb and below_nb are negative, then i lies between neighbouring points. No smoothing will be done.
                    if above_nb <= 0 and below_nb <= 0: 
                        smoothed_counts[i] = vals[i]
                    else:
                        # calculate median of the neighbourhood
                        nbh_med = median(nbh_i)
                        # print(nbh_med)
                        # smooth outlier: if i > points in nbh bring it down
                        if above_nb > 0:
                            # print(nbh_med + sSD)
                            smoothed_counts[i] = nbh_med + sSD
                        # smooth outlier: if i < points in nbh bring it up
                        if below_nb > 0:
                            # print(nbh_med - sSD)
                            smoothed_counts[i] = nbh_med - sSD
        
        # prepare output
        chr_out_data = copy.deepcopy(chr_in_data)

        #check if output smoothed and input values have the same length!
        if len(chr_in_data) != len(smoothed_counts):
            print(f"length of input data ({len(chr_in_data)}) for chromosome {chr} does not equal length of smoothed data ({len(smoothed_counts)})!")
            exit
        else:
            # generate output data and print smoothing results for smoothed bins
            for i in range(len(chr_in_data)): 
                # print smoothing results to out/log   
                if float(chr_in_data[i][-1]) == float(smoothed_counts[i]):
                    continue
                else:
                    print(f"bin {i} of {chr} at position {chr_in_data[i][0]} smoothed from {chr_in_data[i][-1]} to {smoothed_counts[i]}")
                # output data
                chr_out_data[i][-1] = str(smoothed_counts[i])

        # return smoothed_counts
        return chr_out_data

# main function to run above code
def main():

    file_name = os.path.split(read_counts_path)[1].removesuffix('.tsv')
    
    print(f"Smoothening log2R read counts with for {file_name} sl: {smoothing_level},trim: {trim}...")       

    # process input data
    in_data = []
    with open(read_counts_path, "r") as file:
        for line in file:
            fields = line.strip().split("\t")
            in_data.append(fields)

    # Define contig names from input log2r read count file
    chr_names = list(dict.fromkeys([x[1] for x in in_data]))

    # This code only takes 1s to run - no multiprocessing needed. 
    # loop through chromosomes to not smooth across different chromosomes
    smoothedData = []
    for chr in chr_names:
        smoothed_chr = smoothen(chr, in_data, trim, smoothing_level)
        smoothedData.append(smoothed_chr)
    smoothedData = [x for xs in smoothedData for x in xs] 

    # Get results and write out
    outfile = open(f"{outdir}/{file_name}_smoothened_sl{smoothing_level}_t{trim}.tsv", "w")
    for r in smoothedData:
        Line = '\t'.join(r) + '\n'
        outfile.write(Line)        
    outfile.close()

#----
# 3. Run smoothening
#----
if __name__ == "__main__":
    start_t = timeit.default_timer()
    main()
    stop = timeit.default_timer()
    Seconds = round(stop - start_t)
    print(f"Computation time: {Seconds} seconds\n") 


