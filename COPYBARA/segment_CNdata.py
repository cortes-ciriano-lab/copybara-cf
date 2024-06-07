
import argparse
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import logging
from statistics import median
from statistics import mean
from math import sqrt
from math import log2

log = logging.getLogger()
logging.basicConfig(level=logging.WARN)

import timeit
start_t = timeit.default_timer()

#----
# 1. Define functions
#----

# Code adapted from: https://github.com/jeremy9959/cbs
def cbs_stat(x):
    '''Given x, Compute the subinterval x[i0:i1] with the maximal segmentation statistic t. 
    Returns t, i0, i1'''
    
    x0 = x - np.mean(x)
    n = len(x0)
    y = np.cumsum(x0)
    e0, e1 = np.argmin(y), np.argmax(y)
    i0, i1 = min(e0, e1), max(e0, e1)
    s0, s1 = y[i0], y[i1]
    return (s1-s0)**2*n/(i1-i0+1)/(n+1-i1+i0), i0, i1+1


def tstat(x, i):
    '''Return the segmentation statistic t testing if i is a (one-sided)  breakpoint in x'''
    n = len(x)
    s0 = np.mean(x[:i])
    s1 = np.mean(x[i:])
    return (n-i)*i/n*(s0-s1)**2

# def cbs(x, shuffles=1000, p=.05):
# def cbs(x, min_segment_size=5, shuffles=10000, p=.05):
def cbs(x, min_segment_size=5, shuffles=1000, p=.05):
    '''Given x, find the interval x[i0:i1] with maximal segmentation statistic t. Test that statistic against
    given (shuffles) number of random permutations with significance p.  Return True/False, t, i0, i1; True if
    interval is significant, false otherwise.'''

    max_t, max_start, max_end = cbs_stat(x)
    # Check for segment size
    if max_end - max_start < min_segment_size:
        return False, max_t, max_start, max_end
    
    if max_end-max_start == len(x):
        return False, max_t, max_start, max_end
    if max_start < 5:
        max_start = 0
    if len(x)-max_end < 5:
        max_end = len(x)
    thresh_count = 0
    alpha = shuffles*p
    xt = x.copy()
    for i in range(shuffles):
        np.random.shuffle(xt)
        threshold, s0, e0 = cbs_stat(xt)
        if threshold >= max_t:
            thresh_count += 1
        if thresh_count > alpha:
            return False, max_t, max_start, max_end
    return True, max_t, max_start, max_end


def rsegment(x, start, end, L=[], min_segment_size=5, shuffles=1000, p=.05):
    '''Recursively segment the interval x[start:end] returning a list L of pairs (i,j) where each (i,j) is a significant segment.
    '''
    threshold, t, s, e = cbs(x[start:end], min_segment_size=min_segment_size, shuffles=shuffles, p=p)
    log.info('Proposed partition of {} to {} from {} to {} with t value {} is {}'.format(start, end, start+s, start+e, t, threshold))
    if (not threshold) | (e-s < min_segment_size) | (e-s == end-start):
        L.append((start, end))
    else:
        if s > 0:
            rsegment(x, start, start+s, L, min_segment_size=min_segment_size)
        if e-s > 0:
            rsegment(x, start+s, start+e, L, min_segment_size=min_segment_size)
        if start+e < end:
            rsegment(x, start+e, end, L, min_segment_size=min_segment_size)
    return L


def segment(x, min_segment_size=5, shuffles=1000, p=.05):
    '''Segment the array x, using significance test based on shuffles rearrangements and significance level p
    '''
    start = 0
    end = len(x)
    L = []
    rsegment(x, start, end, L, min_segment_size=min_segment_size, shuffles=shuffles, p=p)
    return L


def validate(x, L, shuffles=1000, p=.01):
    S = [x[0] for x in L]+[len(x)]
    SV = [0]
    SV_se = []
    left = 0
    for test, s in enumerate(S[1:-1]):
        t = tstat(x[S[left]:S[test+2]], S[test+1]-S[left])
        log.info('Testing validity of {} in interval from {} to {} yields statistic {}'.format(S[test+1], S[left], S[test+2], t))
        threshold = 0
        thresh_count = 0
        site = S[test+1]-S[left]
        xt = x[S[left]:S[test+2]].copy()
        flag = True
        for k in range(shuffles):
            np.random.shuffle(xt)
            threshold = tstat(xt, site)
            if threshold > t:
                thresh_count += 1
            if thresh_count >= p*shuffles:
                flag = False
                log.info('Breakpoint {} rejected'.format(S[test+1]))
                break
        if flag:
            log.info('Breakpoint {} accepted'.format(S[test+1]))
            SV.append(S[test+1])
            left += 1
    SV.append(S[-1])
    
    for pos, start in enumerate(SV[:-1]):
        end = SV[pos+1]
        # if pos+1 == len(SV):
        #     end = len(x)
        # else:
        #     end = SV[pos+1]
        SV_se.append((start, end))

    return SV, SV_se

# ### Merge segments based on SD threshold and difference of segment means
# def merge_segments(x, Sse, threshold_sd=0.5):
#     '''Merges adjacent segments if their means are not at least threshold_sd standard deviations apart.
    
#     :param x: The original data array.
#     :param L: List of tuples representing the segments (start, end).
#     :param threshold_sd: The number of standard deviations to use as the threshold for merging.
#     :return: A new list of segments after merging.'''
#     if not Sse:
#         return []
    
#     # Calculate the overall standard deviation
#     overall_sd = np.std(x)
#     # Initialize the new list of segments with the first segment
#     new_segments = [Sse[0]]

#     for current_start, current_end in Sse[1:]:
#         # Get the last segment in the new list
#         last_start, last_end = new_segments[-1]
        
#         # Calculate means of the current and last segments
#         last_mean = np.mean(x[last_start:last_end])
#         current_mean = np.mean(x[current_start:current_end])
        
#         # Calculate the difference between the means
#         mean_diff = abs(current_mean - last_mean)
        
#         # If the segments are not at least threshold_sd standard deviations apart, merge them
#         if mean_diff < threshold_sd * overall_sd:
#             # Merge the current segment with the last one in the new list
#             new_segments[-1] = (last_start, current_end)
#         else:
#             # Otherwise, add the current segment as a new entry in the list
#             new_segments.append((current_start, current_end))
    
#     return new_segments

# ### Merge segments based on Xth percentile of changepoints between segments
# def merge_segments(x, Sse, quantile = 0.5):
#     '''Merges adjacent segments if their means are not at least threshold_sd standard deviations apart.
    
#     :param x: The original data array.
#     :param L: List of tuples representing the segments (start, end).
#     :param threshold_sd: The number of standard deviations to use as the threshold for merging.
#     :return: A new list of segments after merging.'''
#     if not Sse:
#         return []
    
#     # Calculate the overall standard deviation
#     # overall_sd = np.std(x)
    
#     med_changepoints = []
    
#     # initiate segment slider with first segment
#     seg_slider = [Sse[0]]

#     for current_start, current_end in Sse[1:]:
#         # Get the last segment in the new list
#         last_start, last_end = seg_slider[-1]
#         # Calculate the difference between the medians
#         median_diff = abs(np.median(x[current_start:current_end]) - np.median(x[last_start:last_end]))
#         # append list collecting median change sizes between segments
#         med_changepoints.append(median_diff)
#         # update seg_slider to move to next segment
#         seg_slider = [(current_start, current_end)]

#     # calculate threshold
#     threshold = np.quantile(med_changepoints, quantile)

#     # Initialize the new list of segments with the first segment
#     new_segments = [Sse[0]]

#     for current_start, current_end in Sse[1:]:
#         # Get the last segment in the new list
#         last_start, last_end = new_segments[-1]
        
#         # Calculate medians of the current and last segments
#         last_median = np.median(x[last_start:last_end])
#         current_median = np.median(x[current_start:current_end])
        
#         # Calculate the difference between the medians
#         median_diff = abs(current_median - last_median)
      
#         # If the segments are not at least threshold apart, merge them
#         if median_diff <= threshold:
#             # Merge the current segment with the last one in the new list
#             new_segments[-1] = (last_start, current_end)
#         else:
#             # Otherwise, add the current segment as a new entry in the list
#             new_segments.append((current_start, current_end))
    
#     return new_segments


### Merge segments based on Xth percentile of changepoints between segments (comparing all to all other segments)
### initially set 0.25 as default quantile; now lowered to 0.2 (20240503)
def merge_segments(x, Sse, quantile = 0.2):
    '''Merges adjacent segments if their medians are not at least threshold_sd standard deviations apart.
    
    :param x: The original data array.
    :param L: List of tuples representing the segments (start, end).
    :param threshold_sd: The number of standard deviations to use as the threshold for merging.
    :return: A new list of segments after merging.'''
    if not Sse:
        return []
    
    # Calculate the overall standard deviation
    # overall_sd = np.std(x)
    # initialise list to collect absolute differences between segment medians
    med_changepoints_diffs = []

    # Iterate through all pairs of segments to calculate absolute median differences
    for i, (start_i, end_i) in enumerate(Sse):
        for j, (start_j, end_j) in enumerate(Sse):
            if j > i:  # Exclude comparing a segment to itself and ensure only following segments are being compared to to avoid duplication of abs differences
                # Calculate the absolute difference between segment medians
                median_diff = abs(np.median(x[start_i:end_i]) - np.median(x[start_j:end_j]))
                med_changepoints_diffs.append(median_diff)
    
    # Calculate the merging threshold based on the specified quantile
    threshold = np.quantile(med_changepoints_diffs, quantile)
    
    # Initialize the new list of segments with the first segment
    new_segments = [Sse[0]]

    for current_start, current_end in Sse[1:]:
        # Get the last segment in the new list
        last_start, last_end = new_segments[-1]
        
        # Calculate medians of the current and last segments
        last_median = np.median(x[last_start:last_end])
        current_median = np.median(x[current_start:current_end])
        
        # Calculate the difference between the medians
        median_diff = abs(current_median - last_median)
      
        # If the segments are not at least threshold apart, merge them
        if median_diff <= threshold:
            # Merge the current segment with the last one in the new list
            new_segments[-1] = (last_start, current_end)
        else:
            # Otherwise, add the current segment as a new entry in the list
            new_segments.append((current_start, current_end))
    
    return new_segments


def draw_segmented_data(data, M, chr_positions, title=None):
    '''Draw a scatterplot of the data with vertical lines at segment boundaries and horizontal lines at medians of 
    the segments. S is a list of segment boundaries.'''
    breaks = [x[0] for x in M]
    breaks.append(M[-1][1])

    chr_breaks = [0]
    for x in chr_positions:
        chr_breaks.append(x[2])

    ticks = [x[1] for x in chr_positions]
    labels = [x[0] for x in chr_positions]

    # sns.set_context("paper", rc={"font.size":6,"axes.titlesize":5,"axes.labelsize":5,"xtick.labelsize": 5,"ytick.labelsize": 5})   

    j=sns.scatterplot(x=range(len(data)),y=data,color='black',s=1.5,legend=None)
    # for x in breaks:
    #     j.axvline(x, linewidth=1)
    for x in chr_breaks:
        j.axvline(x, linewidth=0.5, linestyle='-', color='black')

    for i in range(1,len(breaks)):
        # j.hlines(np.mean(data[breaks[i-1]:breaks[i]]),breaks[i-1],breaks[i],color='green')
        j.hlines(np.median(data[breaks[i-1]:breaks[i]]),breaks[i-1],breaks[i],color='#CC79A7',linewidth=2)
    
    #change axis ticks to chromosomes
    j.set_xticks(ticks)
    j.set_xticklabels(labels, rotation=45)
    

    j.set_title(title)
    j.set_ylim(-3, 3)
    # j.set_ylim(0, 2.5)
    # j.set_ylim(-1, 2)
    j.get_figure().set_size_inches(16,4)
    # j.get_figure().set_size_inches(9,2)

    return j


#########################
# DEFINE INPUT
#########################
parser = argparse.ArgumentParser(description="segment copy number count data... edit (add anscombe optional etc). Also add parameters (optional) for permutation number, pvalues, min_seg size and merging threshold. Also add functionality to include and exclude y and sex chr in general. Also plotting title with info and stats..")
parser.add_argument('-i', '--input_path', type=str, help='Path to prepared (transformed and smoothened) copy number count data file.', required=True)
parser.add_argument('-s', '--sample_prefix', type=str, help='Sample name or prefix to use for out files.', required=True)
parser.add_argument('-o', '--out_dir', type=str, default='.', help='Path to out directory. Default is working directory', required=False)
# parser.add_argument('-sl', '--smoothing_level', type=int, default='10', help='Size of neighbourhood for smoothing.', required=False)
# parser.add_argument('-t', '--trim', type=float, default='0.025', help='Trimming percentage to be used.', required=False)

args = parser.parse_args()

input_path = args.input_path
prefix = args.sample_prefix
outdir = args.out_dir





if __name__ == '__main__':

    log.setLevel(logging.INFO)
    # sample = generate_normal_time_series(5)
    # sample = np.array([-0.663, -0.725, -0.666, -0.752, -0.688, -0.598, -0.721, 0.491, 0.468, 0.512, 0.517, 0.463])
    # path = 'OUT2/smooth_test_s2_20231019_full_level10_sd4-2_t0.025_anscombe.tsv'
    # path = 'DEV_20240220/out_test_20240220/CCSBEST325_read_counts_log2r_normalised_smoothened_level10_sd4-2_t0.025.tsv'
    
    # path = 'DEV_20240220/out_test_20240220/CCSBEST325_read_counts_anscombe_normalised_smoothened_level10_sd4-2_t0.025.tsv'
    # path = 'DEV_20240220/out_test2_20240227/CCSBEST328_read_counts_anscombe_normalised_smoothened_level10_sd4-2_t0.025.tsv'
    # path = 'DEV_20240220/out_test3_colo829_20240228/colo829_reads0.5mio_TF0.75_read_counts_anscombe_normalised_smoothened_level10_sd4-2_t0.025.tsv'
    # path = 'DEV_20240220/out_test_20240229/CCSBEST325_read_counts_log2r_smoothened_level10_sd4-2_t0.025_invlog2_anscombe.tsv'
    
    # path = 'CNAPS_presentation/cn_out/smo0206/SMO0206_read_counts_log2r_smoothened_level10_sd4-2_t0.025_v2.tsv'
    with open(input_path, "r") as file:
        original_data = []
        input_ansc = []
        input_log2 = []
        for line in file:

            fields = line.strip().split("\t")
            original_data.append(fields)

            val = float(fields[-1])
            input_log2.append(val)
            # anscombe transform to stabilise variance!
            val_ansc = sqrt((2**val) + 3/8) # inverse log2 then anscombe transform
            input_ansc.append(val_ansc)
           
            # print(fields)
            # print(val_ansc)
    # print(original_data)
    # print(type(original_data))
    # print(input)

    # shuf = 10000

    ### VARIANCE STABILISATION --> ANSCOMBE

    L = segment(input_ansc, shuffles=1000, p = 0.1)
    # print(L)
    
    S, Sse = validate(input_ansc, L, shuffles=1000, p = 0.05)
    # print(S)

    # # M = merge_segments(sample, Sse, threshold_sd=0.5)
    M = merge_segments(input_ansc, Sse, quantile = 0.2)
    # M = merge_segments(input_ansc, Sse, quantile = 0)

    # print(M)

    # BUILD OUTPUT DATA
    out_data_full = []
    # out_data_segmented = [] ### ADD THIS IN LATER!
    for i, (s_start,s_end) in enumerate(M):
        # print(i, s_start, s_end)
        seg=str(f"seg{i+1}")
        # for row in original_data[s_start,s_end]:
        #     val = float(row[-1])
        cur_seg = original_data[s_start:s_end]
        seg_vals = [float(row[-1]) for row in cur_seg]
        seg_median = median(seg_vals)
        
        for row in cur_seg:
            # print(type(row))
            row.append(seg)
            row.append(str(seg_median))
            out_data_full.append(row)

    
    ### WRITE OUT FILE ###
    outfile = open(f"{outdir}/{prefix}_segmented_copy_number_log2r.tsv", "w")
    for r in out_data_full:
        Line = '\t'.join(r) + '\n'
        outfile.write(Line)
    outfile.close()

    
    # define chromosomes and positions for plotting
    # chr_unique = np.unique([x[1] for x in original_data if x[1] != "chrY"])
    chr_unique = np.unique([x[1] for x in original_data])
    indeces = [(index, sublist[1]) for index, sublist in enumerate(original_data)]
    # print(chr_unique)
    # print(indeces)
    chr_positions = []
    for chr in chr_unique:
        tmp = [x for x in indeces if x[1] == chr]
        chr_break = max(x[0] for x in tmp)
        midpoint = mean([x[0] for x in tmp])
        chr_positions.append([chr, midpoint, chr_break])
    
    # print(chr_positions)
    title=str(f"copy number (log2R) of sample {prefix} (No. segments = {len(M)})")
    ax = draw_segmented_data(input_log2,  M, chr_positions, title=title)
    # ax.tick_params(axis='x', rotation=45)


    # # # ax.get_figure().savefig('DEV_20240220/out_test_20240220/CCSBEST325_read_counts_log2r_normalised_smoothened_level10_sd4-2_t0.025_p0.1valp0.05_t2.png')
    
    # ax.get_figure().savefig('DEV_20240220/out_test_20240220/CCSBEST325_read_counts_anscombe_normalised_smoothened_level10_sd4-2_t0.025_p0.1valp0.05_minseg5_merge_changept_ALL_quant0.2_median.v2.png')
    # ax.get_figure().savefig('DEV_20240220/out_test2_20240227/CCSBEST328_read_counts_anscombe_normalised_smoothened_level10_sd4-2_t0.025_p0.1valp0.05_minseg5_merge_changept_ALL_quant0.25_median.v2.png')
    # ax.get_figure().savefig('DEV_20240220/out_test3_colo829_20240228/colo829_reads0.5mio_TF0.75_read_counts_anscombe_normalised_smoothened_level10_sd4-2_t0.025_p0.1valp0.05_minseg5_merge_changept_ALL_quant0.25_median.v2.png')
    # ax.get_figure().savefig('DEV_20240220/out_test_20240229/CCSBEST325_read_counts_log2r_smoothened_level10_sd4-2_t0.025_invlog2_anscombe_p0.1valp0.05_minseg5_merge_changept_ALL_quant0.25_median.v2.png')
    
    # ax.get_figure().savefig('CNAPS_presentation/cn_out/smp0003/SMP0003_read_counts_log2r_smoothened_level10_sd4-2_t0.025_invlog2_anscombe_p0.1valp0.05_minseg5_merge_changept_ALL_quant0.25_median.v2.png')

    ax.get_figure().savefig(f'{outdir}/{prefix}_segmented_copy_number_log2r_plot.png')



stop = timeit.default_timer()
Seconds = round(stop - start_t)
print(f"Computation time: {Seconds} seconds\n") 