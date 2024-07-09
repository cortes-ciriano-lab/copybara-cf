# temporary script to build ploidy code
##########################################

import argparse
import statistics
import math
import copy

parser = argparse.ArgumentParser(description="Estimate purity and ploidy and fit absolute copy number.")
parser.add_argument('-i', '--input', type=str, help='Path to segmented log2r copy number tsv file.', required=True)
parser.add_argument('-pur', '--purity', type=float, help='estimated purity from step before.', required=True)
parser.add_argument('-s', '--sample_prefix', type=str, help='Sample name of prefix to use for out files.', required=False)
parser.add_argument('-o', '--out_dir', type=str, default='.', help='Path to out directory. Default is working directory', required=False)
args = parser.parse_args()

input_file = args.input
cellularity = args.purity
prefix = args.sample_prefix
outdir = args.out_dir

## PARAMETERS ##
min_ploidy = 1.5
max_ploidy = 7
ploidy_step = 0.01
min_cellularity = max(0,cellularity - 0.1)
max_cellularity = min(1,cellularity + 0.1) 
cellularity_step = 0.01
distance_function = "RMSD"

distance_filter_scale_factor = 1.25
max_proportion_zero = 0.1
min_proportion_close_to_whole_number = 0.5
max_distance_from_whole_number = 0.25

####
# FUNCTIONS
####

# def relative_to_absolute_CN(relative_CN, purity, ploidy):
#     acn = ploidy + (relative_CN - 1)*(ploidy+(2/purity)-2) 
#     return acn

# def acn_distance(relative_CN, purity, ploidy, distance_function="RMSD", weights=None): # function needs list of relative seg CNs and list of weights
#     acn = [relative_to_absolute_CN(x, purity, ploidy) for x in relative_CN]
#     differences = [abs(x - round(x)) for x in acn]
#     if weights == None:
#         weights = [1]*len(relative_CN)
#     if distance_function == "MAD":
#         # estimate weighted distances using mean absolute deviation (MAD)
#         distance = sum(differences[i] * weights[i] for i in range(len(differences))) / sum(weights)
#     if distance_function == "RMSD":
#         distance = math.sqrt(sum(differences[i]**2 * weights[i] for i in range(len(differences))) / sum(weights))
#     return distance


# def define_search_space(min, max, by=0.01):
#     digs = len(str(by))-2 if isinstance(by,int) != True else 1
#     # print(by,digs)
#     if by <= 0:
#         raise ValueError("by must be >= zero")
#     result = []
#     current = min
#     if (min < max and by > 0):
#         while (current <= max):
#             result.append(current)
#             current += by
#             current = round(current,digs)
#     return result


# def build_search_grid(purity_seq, ploidy_seq):
#     unique_combinations = []
#     for i in range(len(purity_seq)):
#         for j in range(len(ploidy_seq)):
#             unique_combinations.append([purity_seq[i], ploidy_seq[j], i+1, j+1])
#     return unique_combinations


# def estimate_grid_distances(min_cellularity, max_cellularity, cellularity_step, min_ploidy, max_ploidy, ploidy_step, relative_CN, weights=None, distance_function="RMSD"):
#     purs = define_search_space(min_cellularity,max_cellularity,by=cellularity_step)
#     plois = define_search_space(min_ploidy,max_ploidy,by=ploidy_step)
#     grid = build_search_grid(purs,plois)
#     fits = []
#     for pair in grid:
#         cur_pur = pair[0]
#         cur_ploi = pair[1]
#         d = acn_distance(relative_CN, cur_pur, cur_ploi, weights=weights, distance_function=distance_function)
#         pair.append(d)
#         fits.append(pair)
#     return fits


# def reduce_grid(fits,distance_filter_scale_factor = 1.25):
#     reduced_grid = copy.deepcopy(fits)
#     for xdelta in [-1, 0, 1]:
#         for ydelta in [-1, 0, 1]:
#             if xdelta != 0 or ydelta != 0:
#                 keep_idxs = []
#                 for idx,sol in enumerate(reduced_grid):
#                     xc = sol[-2] + xdelta
#                     yc = sol[-3] + ydelta
#                     # Check if there's a corresponding (xc, yc) in fits
#                     dc = next((d[-1] for d in fits if d[-2] == xc and d[-3] == yc),None)
#                     if dc is None or sol[-1] <= dc:
#                         keep_idxs.append(idx)
#                 reduced_grid = [reduced_grid[i] for i in keep_idxs]
#     if isinstance(distance_filter_scale_factor, (int,float)):
#         scaler = min([x[-1] for x in reduced_grid]) * distance_filter_scale_factor
#         reduced_grid = [x for x in reduced_grid if x[-1] < scaler]
#     return reduced_grid


# def is_acceptable_fit(purity, ploidy, relative_CN, weights, max_proportion_zero = 0.1, min_proportion_close_to_whole_number = 0.5, max_distance_from_whole_number = 0.25):
#     acn = [relative_to_absolute_CN(x, purity, ploidy) for x in relative_CN]
#     acn_int = [round(x) for x in acn]
#     differences = [abs(x - round(x)) for x in acn]
#     if weights == None:
#         weights = [1]*len(relative_CN)
#     # Filter fits based on fraction genome fitted to ACN of >= 0
#     zeros_idxs =  [i for i in range(len(acn_int)) if acn_int[i] <= 0]
#     prop_zero = sum([weights[i] for i in zeros_idxs]) / sum(weights)
#     if prop_zero > max_proportion_zero:
#         return False
#     # Filter fits based on fraction genome that is fitted to < min_proportion_close_to_whole_number
#     state_idxs = [i for i in range(len(differences)) if differences[i] < max_distance_from_whole_number] 
#     prop_state = sum([weights[i] for i in state_idxs]) / sum(weights)
#     if prop_state < min_proportion_close_to_whole_number:
#         return False
#     # Remove fits which result in overstretchin/overfitting of copy number states (i.e. copy number state skipping) normally caused by oversegmentation
#     most_common = statistics.mode(acn_int)
#     second_most_common = statistics.mode([x for x in acn_int if x != most_common])
#     if abs(most_common - second_most_common) >= 2:
#         return False
#     else:
#         return True


# def viable_solutions(fits_r, relative_CN, weights, max_proportion_zero = 0.1, min_proportion_close_to_whole_number = 0.5, max_distance_from_whole_number = 0.25):
#     solutions = []
#     for sol in fits_r:
#         purity,ploidy = sol[0],sol[1]
#         if is_acceptable_fit(purity, ploidy, relative_CN, weights, 
#                         max_proportion_zero = max_proportion_zero, 
#                         min_proportion_close_to_whole_number = min_proportion_close_to_whole_number, 
#                         max_distance_from_whole_number = max_distance_from_whole_number) == True:
#             solutions.append(sol)
#     # sort solutions by distance function
#     return solutions


# def rank_solutions(solutions,distance_precision=3):
#     ranked = copy.deepcopy(solutions)
#     # round distance function
#     for x in ranked:
#         x[-1] = round(x[-1],distance_precision)
#     ranked = sorted(ranked,key=lambda x: (x[-1],x[1])) # sort by rounded distance, and second by ploidy
#     # add rank position to output
#     for idx,x in enumerate(ranked):
#         x.append(idx+1)
#     return ranked


################### 
# READ IN AND PREPARE DATA
###################

rel_copy_number = []
with open(input_file, "r") as file:
    for line in file:
        fields = line.strip().split("\t")
        chrom,start,end,seg_id = fields[1],int(fields[2]),int(fields[3]),fields[-2]
        copy_number=2**float(fields[-1]) # revert log2 of relatitve copy number
        bin_length = end-start+1 
        # print(chrom,start,end,seg_id,copy_number)
        rel_copy_number.append([chrom,start,end,bin_length,seg_id,copy_number])

# Collapse copy number to segments
segs = list(dict.fromkeys([x[-2] for x in rel_copy_number]))

# estimate median length of all bins to determine weight of segment
med_length = statistics.median([x[3] for x in rel_copy_number])

rel_cn_segs = []
for s in segs:
    cn_seg = [x for x in rel_copy_number if x[-2] == s]
    bin_count = len(cn_seg)
    sum_of_bin_lengths = sum([x[3] for x in cn_seg])
    weight = sum_of_bin_lengths/med_length
    min_start = min([x[1] for x in cn_seg])
    max_end = max([x[2] for x in cn_seg])
    # sanity check copy number values for each segment are all the same!
    if all(i[-1] == cn_seg[0][-1] for i in cn_seg) == True:
        rel_cn_segs.append([cn_seg[0][0], min_start, max_end, s, bin_count, sum_of_bin_lengths, weight, cn_seg[0][-1]])
    else:
        print(f"    ERROR: segment {s} contains multiple copy number values.")
        break

# perform grid search and estimate distances for each purity-ploidy solution in grid
## functions to translate
## DONE 1. relative_to_absolute_copy_number
## DONE 2. weighted.mean / rmsd...
## DONE 3. absolute copy number distance
## 4. absolute copy number distance grid
## 5. grid reduction
## 6. viable solution
## 7. convert segmented RCN to ACN and write all out...

relative_CN = [x[-1] for x in rel_cn_segs]
weights = [x[-2] for x in rel_cn_segs]
fits = estimate_grid_distances(min_cellularity, max_cellularity, cellularity_step, min_ploidy, max_ploidy, ploidy_step, relative_CN, weights=weights, distance_function=distance_function)
fits_r = reduce_grid(fits,distance_filter_scale_factor = distance_filter_scale_factor)
solutions = viable_solutions(fits_r, relative_CN, weights, 
                      max_proportion_zero = max_proportion_zero, 
                      min_proportion_close_to_whole_number = min_proportion_close_to_whole_number, 
                      max_distance_from_whole_number = max_distance_from_whole_number)
solutions_ranked = rank_solutions(solutions,distance_precision=3)
final_fit = solutions_ranked[0]


