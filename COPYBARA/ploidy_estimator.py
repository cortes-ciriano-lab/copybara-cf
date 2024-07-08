# temporary script to build ploidy code
##########################################

import argparse
import statistics
import math

parser = argparse.ArgumentParser(description="Estimate purity and ploidy and fit absolute copy number.")
parser.add_argument('-i', '--input', type=str, help='Path to segmented log2r copy number tsv file.', required=True)
parser.add_argument('-pur', '--purity', type=float, help='estimated purity from step before.', required=True)
parser.add_argument('-s', '--sample_prefix', type=str, help='Sample name of prefix to use for out files.', required=False)
parser.add_argument('-o', '--out_dir', type=str, default='.', help='Path to out directory. Default is working directory', required=False)
args = parser.parse_args()

input_file = args.input
purity = args.purity
prefix = args.sample_prefix
outdir = args.out_dir

## PARAMETERS ##
min_ploidy = 1.5
max_ploidy = 7
ploidy_step = 0.01
min_cellularity = max(0,purity - 0.1)
max_cellularity = min(1,purity + 0.1) 
cellularity_step = 0.01
distance_function = "RMSD"

distance_filter_scale_factor = 1.25
max_proportion_zero = 0.1
min_proportion_close_to_whole_number = 0.5
max_distance_from_whole_number = 0.25
solution_proximity_threshold = 5


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

####
# FUNCTIONS
####

def relative_to_absolute_CN(relative_CN, purity, ploidy):
    acn = ploidy + (relative_CN - 1)*(ploidy+(2/purity)-2) 
    return acn

def acn_distance(relative_CN, purity, ploidy, distance_function="RMSD", weights=None): # function needs list of relative seg CNs and list of weights
    acn = [relative_to_absolute_CN(x, purity, ploidy) for x in relative_CN]
    differences = [abs(x - round(x)) for x in acn]
    if weights == None:
        weights = [1]*len(relative_CN)
    if distance_function == "MAD":
        # estimate weighted distances using mean absolute deviation (MAD)
        distance = sum(differences[i] * weights[i] for i in range(len(differences))) / sum(weights)
    if distance_function == "RMSD":
        distance = math.sqrt(sum(differences[i]**2 * weights[i] for i in range(len(differences))) / sum(weights))
    return distance

def distance_grid(relative_CN, weights, min_ploidy, max_ploidy, ploidy_step, min_cellularity, max_cellularity, cellularity_step, distance_function="RMSD"):
    distances = []
    return distances
