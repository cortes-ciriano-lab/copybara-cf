import argparse
from math import sqrt
from math import log2
from scipy.stats import norm
from statistics import median

import timeit

start_t = timeit.default_timer()

#----
# 1. Define functions
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
def smoothen(normalised_data, trim, smoothing_level):
    if smoothing_level <= 0:
        return normalised_data  # Skip smoothing and use original (normalised) data for CBS
    else:
        oSD = sqrt(trimmed_variance(normalised_data, trim)) * 4
        # print(f"oSD = {oSD}")
        sSD = sqrt(trimmed_variance(normalised_data, trim)) * 2
        # print(f"sSD = {sSD}")
        k = smoothing_level
        start = 0
        end = len(normalised_data)
        # print(oSD, sSD, k, start, end)
        smoothed_counts = [0.0] * len(normalised_data)

        for i in range(len(normalised_data)):
            # print(f"i = {i}, {normalised_data[i]}")
            nbh_start = max(start, i - k)
            nbh_end = min(end, i + k)
            
            above_nb = 100 * oSD #initate as high number
            below_nb = 100 * oSD #initate as high number
            
            # neighbourhood array
            nbh_i = normalised_data[nbh_start:nbh_end+1]
            # print(nbh_i)

            for nb in range(nbh_start, nbh_end):
                if nb != i:
                    dist_inb = normalised_data[i] - normalised_data[nb]
                    # print(f"dist_inb = {dist_inb}")
                    if abs(dist_inb) <= oSD:
                        smoothed_counts[i] = normalised_data[i]
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
                        smoothed_counts[i] = normalised_data[i]
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
        
        return smoothed_counts


#########################
# DEFINE INPUT
#########################
parser = argparse.ArgumentParser(description="smoothen copy number count data (Log2R)")
parser.add_argument('-rc', '--read_counts', type=str, help='Path to read_counts_log2r.tsv file.', required=True)
parser.add_argument('-s', '--sample_prefix', type=str, help='Sample name or prefix to use for out files.', required=True)
parser.add_argument('-o', '--out_dir', type=str, default='.', help='Path to out directory. Default is working directory', required=False)
parser.add_argument('-sl', '--smoothing_level', type=int, default='10', help='Size of neighbourhood for smoothing.', required=False)
parser.add_argument('-t', '--trim', type=float, default='0.025', help='Trimming percentage to be used.', required=False)

args = parser.parse_args()

read_counts_path = args.read_counts
prefix = args.sample_prefix
outdir = args.out_dir
smoothing_level = args.smoothing_level
trim = args.trim


### TESTING ###

# path = 'OUT/smooth_test_full_data2_read_counts_normalised.tsv'
# path = 'OUT/smooth_test_full_data_read_counts_normalised_log.tsv'
# path = 'OUT2/sample2_test_read_counts_normalised.tsv'
# path = 'DEV_20240220/out_test_20240220/CCSBEST325_read_counts_anscombe_with0s_normalised.tsv'
# path = 'DEV_20240220/out_test2_20240227/CCSBEST328_read_counts_anscombe_normalised.tsv'
# path = 'DEV_20240220/out_test3_colo829_20240228/colo829_reads0.5mio_TF0.75_read_counts_anscombe_normalised.tsv'
# path = 'DEV_20240220/out_test_20240229/CCSBEST325_read_counts_log2r.tsv'
# path = 'CNAPS_presentation/cn_out/smo0206/SMO0206_read_counts_log2r.tsv'
with open(read_counts_path, "r") as file:
    original_data = []
    input = []
    for line in file:
        fields = line.strip().split("\t")
        val = float(fields[-1])
        # print(val)
        original_data.append(fields)
        input.append(val)
# print(input)

# run 
        
# smoothing_level = 10
# trim = 0.025
print(f"sl: {smoothing_level},trim: {trim}")

smoothed = smoothen(input, trim, smoothing_level)


# for i in range(len(input)):
#     if float(original_data[i][-1]) == float(smoothed[i]):
#         continue
#     else:
#         print(f"old: {original_data[i][-1]}\tnew: {smoothed[i]}")

# out data
output = original_data

#check if output smoothed and input values are all the same length!
if len(original_data) != len(smoothed):
    print(f"length of input data ({len(original_data)} does not equal length of smoothed data ({len(smoothed)}))")
    exit
else:
    # generate output data and print smoothing results for smoothed bins
    for i in range(len(input)): 
        # print smoothing results to out    
        if float(original_data[i][-1]) == float(smoothed[i]):
            continue
        else:
            print(f"bin {i} at position {original_data[i][0]} smoothed from {original_data[i][-1]} to {smoothed[i]}")
        # output data
        output[i][-1] = str(smoothed[i])

# outfile = open(f"DEV_20240220/out_test_20240220/CCSBEST325_read_counts_log2r_normalised_smoothened_level10_sd4-2_t0.025.tsv", "w")
# outfile = open(f"DEV_20240220/out_test2_20240227/CCSBEST328_read_counts_anscombe_normalised_smoothened_level10_sd4-2_t0.025.tsv", "w")
# outfile = open(f"DEV_20240220/out_test_20240229/CCSBEST325_read_counts_log2r_smoothened_level10_sd4-2_t0.025.tsv", "w")

# print(f"{outdir}/{prefix}_read_counts_log2r_smoothened_sl{smoothing_level}_t{trim}.tsv")
outfile = open(f"{outdir}/{prefix}_read_counts_log2r_smoothened_sl10_t0.025.tsv", "w")
for r in output:
    Line = '\t'.join(r) + '\n'
    # Line = str(r) + '\n'
    # print(Line)
    outfile.write(Line)        
outfile.close()





# anscombe sqrt transform
# outfile2 = open(f"DEV_20240220/out_test_20240229/CCSBEST325_read_counts_log2r_smoothened_level10_sd4-2_t0.025_invlog2_anscombe.tsv", "w")

# outfile2 = open(f"CNAPS_presentation/cn_out/smo0206/SMO0206_read_counts_log2r_smoothened_level10_sd4-2_t0.025_invlog2_anscombe.tsv", "w")
# # anscombe = smoothed
# for r in smoothed:
#     # print(r, log(r))
#     # test = sqrt(float(r) + 3/8)
#     r = str(sqrt((2**r) + 3/8)) # inverse log2 then anscombe transform
#     Line = str(r) + '\n'
#     # print(Line)
#     outfile2.write(Line)        
# outfile2.close()

# anscombe sqrt transform followed by log2
# outfile3 = open(f"DEV_20240220/out_test_20240229/CCSBEST325_read_counts_log2r_smoothened_level10_sd4-2_t0.025_invlog2_anscombe_log2.tsv", "w")

# outfile3 = open(f"CNAPS_presentation/cn_out/smo0206/SMO0206_read_counts_log2r_smoothened_level10_sd4-2_t0.025_invlog2_anscombe_log2.tsv", "w")
# for r in smoothed:
#     # print(r)
#     r = str(log2(sqrt((2**r) + 3/8))) # log2 anscombe transformed data
#     Line = str(r) + '\n'
#     # print(Line)
#     outfile3.write(Line)        
# outfile3.close()

stop = timeit.default_timer()
Seconds = round(stop - start_t)
print(f"Computation time: {Seconds} seconds\n") 