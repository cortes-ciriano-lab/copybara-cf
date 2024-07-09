import numpy as np
from scipy.stats import gaussian_kde
import statistics
import math
import copy
import argparse

import CN_functions as cnfitter

#----
# 1. Parse command line arguments
#----
parser = argparse.ArgumentParser(description="Estimate purity and ploidy and fit absolute copy number.")
parser.add_argument('-i', '--input', type=str, help='Path to segmented log2r copy number tsv file.', required=True)
parser.add_argument('-a', '--allele_counts', type=str, default=None, help='Path to hetSNPs allele counts to estimate purity.', required=False)
parser.add_argument('-s', '--sample_prefix', type=str, help='Sample name of prefix to use for out files.', required=True)
parser.add_argument('-o', '--out_dir', type=str, default='.', help='Path to out directory. Default is working directory', required=False)
# Fitting parameters
parser.add_argument('--min_ploidy', type=float, default=1.5, help='Minimum ploidy to be considered for copy number fitting.', required=False)
parser.add_argument('--max_ploidy', type=float, default=6.5, help='Maximum ploidy to be considered for copy number fitting.', required=False)
parser.add_argument('--ploidy_step', type=float, default=0.01, help='Ploidy step size for grid search space used during for copy number fitting.', required=False)
parser.add_argument('--min_cellularity', type=float, default=0.2, help='Minimum cellularity to be considered for copy number fitting. If hetSNPs allele counts are provided, this is estimated during copy number fitting. Alternatively, a purity value can be provided if the purity of the sample is already known.', required=False)
parser.add_argument('--max_cellularity', type=float, default=1, help='Maximum cellularity to be considered for copy number fitting. If hetSNPs allele counts are provided, this is estimated during copy number fitting. Alternatively, a purity value can be provided if the purity of the sample is already known.', required=False)
parser.add_argument('--cellularity_step', type=float, default=0.01, help='Ploidy step size for grid search space used during for copy number fitting.', required=False)
parser.add_argument('--distance_function', type=str, default='RMSD', help='Distance function to be used for copy number fitting.', choices=['RMSD', 'MAD'], required=False)
parser.add_argument('--distance_filter_scale_factor', type=float, default=1.25, help='Distance filter scale factor to only include solutions with distances < scale factor * min(distance).', required=False)
parser.add_argument('--distance_precision', type=int, default=3, help='Number of digits to round distance functions to', required=False)
parser.add_argument('--max_proportion_zero', type=float, default=0.1, help='Maximum proportion of fitted copy numbers to be tolerated in the zero or negative copy number state.', required=False)
parser.add_argument('--min_proportion_close_to_whole_number', type=float, default=0.5, help='Minimum proportion of fitted copy numbers sufficiently close to whole number to be tolerated for a given fit.', required=False)
parser.add_argument('--max_distance_from_whole_number', type=float, default=0.25, help='Distance from whole number for fitted value to be considered sufficiently close to nearest copy number integer.', required=False)
parser.add_argument('--min_ps_size', type=int, default=10, help='Minimum size (number of SNPs) for phaseset to be considered for purity estimation.', required=False)
parser.add_argument('--min_ps_length', type=int, default=500000, help='Minimum length (bps) for phaseset to be considered for purity estimation.', required=False)

args = parser.parse_args()

min_ps_size = 10
###############
# ADD IN FITTING PARAMETERS... i.e. PS parameters for purity and ploidy estimation parameters...
###############

input_file = args.input
allele_counts_file = args.allele_counts
prefix = args.sample_prefix
outdir = args.out_dir

min_ploidy = args.min_ploidy
max_ploidy = args.max_ploidy
ploidy_step = args.ploidy_step
min_cellularity = args.min_cellularity
max_cellularity = args.max_cellularity
cellularity_step = args.cellularity_step

distance_function = args.distance_function
distance_filter_scale_factor = args.distance_filter_scale_factor
distance_precision = args.distance_precision
max_proportion_zero = args.max_proportion_zero
min_proportion_close_to_whole_number = args.min_proportion_close_to_whole_number
max_distance_from_whole_number = args.max_distance_from_whole_number

min_ps_size = args.min_ps_size
min_ps_length = args.min_ps_length

#----
# 2. Define functions
#----
def process_allele_counts(allele_counts_file):
    allele_counts = []
    with open(allele_counts_file, "r") as file:
        for line in file:
            fields = line.strip().split("\t")
            A,C,G,T,N = int(fields[5]),int(fields[6]),int(fields[7]),int(fields[8]),int(fields[9])
            ps_id = f"{fields[0]}_{fields[-1]}" # assign phase sets IDs as same phase set might accur in different chromosomes
            fields.append(ps_id)
            DP = A+C+G+T+N
            if DP > 20: # only include het SNPs allele counts with a total depth of > 20.
                fields.append(DP)
                # print(fields.append(str(DP)))
                allele_counts.append(fields)
    # columns:"chr", "start", "end","REF", "ALT", "A","C","G","T","N","AF_0","AF_1","GT","PS","ps_id","DP"
    # get unique phase sets to iterate (or multiprocess through)
    phasesets = list(dict.fromkeys([x[-2] for x in allele_counts]))
    return allele_counts, phasesets


def estimate_cellularity_phased_hetSNPs(phasesets,allele_counts,dp_cutoff,min_ps_size=10,min_ps_length=500000):
    '''
    Estimate sample cellularity/purity using allele counts of phased hetSNPs.
    '''
    phasesets_dict = {}
    ps_summary = []
    for ps in phasesets:
        if "None" not in ps: 
            print(ps)
            # subset allele count data to ps and remove SNPs with AF0/AF1 of 0 or 1
            ac_ps = [x for x in allele_counts if x[-2] == ps and float(x[10]) != 0 and float(x[10]) != 1 and float(x[11]) != 0 and float(x[11]) != 1]
            if len(ac_ps) < min_ps_size: # skip ps with less than min_ps_size SNPs
                continue
            if (max([int(x[2]) for x in ac_ps]) - min([int(x[1]) for x in ac_ps])) < min_ps_length: # skip ps less than min_ps_length bps in size
                continue
            else:
                ps_length = max([int(x[2]) for x in ac_ps]) - min([int(x[1]) for x in ac_ps]) # genomic lenght of ps_length
                ps_snps = len(ac_ps) # number of SNPs in ps
                ps_depth = statistics.mean([int(x[-1]) for x in ac_ps]) # mean depth of ps
                # ps_weight = ps_snps/no_hetSNPs # estimate weight of PS based on number of hetSNPs present in PS compared to all hetSNPs across sample
                # all AFs at all het SNP position within ps
                af = [float(x[10]) for x in ac_ps] + [float(x[11]) for x in ac_ps]
                # check if distribution is not unimodal and estimate bimodality coefficient
                if cnfitter.is_unimodal(af) == False:
                    ps_bc = cnfitter.bimodality_coefficient(af)
                    if ps_depth < dp_cutoff:
                        phasesets_dict[ps] = ac_ps
                        ps_summary.append([ps, ps_depth, ps_length, ps_snps, ps_bc])
                        # ps_summary.append([ps, ps_depth, ps_length, ps_snps, ps_weight, ps_bc])
                
    ps_summary = sorted(ps_summary,key=lambda x: x[-1]) # sort ps by bimodality coefficient
    # ps_summary = sorted(sorted(ps_summary,key=lambda x: x[-2]), key=lambda x: round(x[-1],3)) # sort by ps_bc and weight
    # select top 10% PSs
    # # ps_top = [x[0] for x in ps_summary[-(round(len(ps_summary)*0.01)):]]
    # ps_summary_long = [x for x in ps_summary if x[2] > min_ps_length]
    # ps_top = [x[0] for x in ps_summary_long[-10:]]
    ps_top = [x[0] for x in ps_summary[-10:]]
    # pull out data for ps_top
    ps_top_acs = []
    for pst in ps_top:
        ps_top_acs += phasesets_dict[pst]
    # Estimate cellularity
    af_cutoff = 0.5
    pur0 = statistics.median([float(x[10]) for x in ps_top_acs if float(x[10]) > af_cutoff] + [float(x[11]) for x in ps_top_acs if float(x[11]) > af_cutoff])
    pur1 = statistics.median([float(x[10]) for x in ps_top_acs if float(x[10]) < (1-af_cutoff)] + [float(x[11]) for x in ps_top_acs if float(x[11]) < (1-af_cutoff)])
    cellularity = statistics.mean([1-(1-max(pur0,pur1))*2]+[1-(min(pur0,pur1))*2])
    return cellularity

def process_log2r_input(input_file):
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
    segs = list(dict.fromkeys([x[-2] for x in rel_copy_number])) # unique segment ids  
    med_length = statistics.median([x[3] for x in rel_copy_number]) # Estimate median length of all bins to determine weight of segment
    
    rel_copy_number_segments = []
    for s in segs:
        cn_seg = [x for x in rel_copy_number if x[-2] == s]
        bin_count = len(cn_seg)
        sum_of_bin_lengths = sum([x[3] for x in cn_seg])
        weight = sum_of_bin_lengths/med_length
        min_start = min([x[1] for x in cn_seg])
        max_end = max([x[2] for x in cn_seg])
        # sanity check copy number values for each segment are all the same!
        if all(i[-1] == cn_seg[0][-1] for i in cn_seg) == True:
            rel_copy_number_segments.append([cn_seg[0][0], min_start, max_end, s, bin_count, sum_of_bin_lengths, weight, cn_seg[0][-1]])
        else:
            print(f"    ERROR: segment {s} contains multiple copy number values.")
            break

    return rel_copy_number_segments

#----
# 3. Estimate purity using phased heterozygous SNPs
### Method/principles based on: https://doi.org/10.1371/journal.pone.0045835
#----

# Open and prepare data to estimate cellularity if phased hetSNPs allele counts were provided.
if allele_counts_file == None:
    print("     ... No phased hetSNPs allele counts were provided. Cellularity is estimated using grid search. Note that this is less accurate. Please provide hetSNPs allele counts if possible. See documentation for instructions on how to generate these.")
elif allele_counts_file != None:
    print("     ... Allele counts for phased hetSNPs provided and being processed to estimate sample purity ...")
    
    allele_counts, phasesets = process_allele_counts(allele_counts_file)
    # Estimate depth cutoff to exclude potentially amplified/gained regions as those will impact BAF distribution
    ac_mean_depth = statistics.mean([x[-1] for x in allele_counts])
    dp_cutoff = 2 * ac_mean_depth 
    # no_hetSNPs = len(allele_counts)

    cellularity = estimate_cellularity_phased_hetSNPs(phasesets, allele_counts, dp_cutoff, min_ps_size=min_ps_size, min_ps_length=min_ps_size)

    min_cellularity = max(0,cellularity - 0.1)
    max_cellularity = min(1,cellularity + 0.1) 

#----
# 4. Estimate ploidy and fit ACN using estimated sample purity
### Method/principles based on rascal R package: https://www.biorxiv.org/content/10.1101/2021.07.19.452658v1
#----

rel_copy_number_segments = process_log2r_input(input_file)

# Prepare input for copy number fitting
relative_CN = [x[-1] for x in rel_copy_number_segments]
weights = [x[-2] for x in rel_copy_number_segments]

# Copy number fitting
fits = cnfitter.estimate_grid_distances(min_cellularity, max_cellularity, cellularity_step, min_ploidy, max_ploidy, ploidy_step, relative_CN, weights=weights, distance_function=distance_function)
fits_r = cnfitter.reduce_grid(fits,distance_filter_scale_factor = distance_filter_scale_factor)
solutions = cnfitter.viable_solutions(fits_r, relative_CN, weights, 
                      max_proportion_zero = max_proportion_zero, 
                      min_proportion_close_to_whole_number = min_proportion_close_to_whole_number, 
                      max_distance_from_whole_number = max_distance_from_whole_number)
solutions_ranked = cnfitter.rank_solutions(solutions,distance_precision=3)
final_fit = solutions_ranked[0]

# Convert segmented relative to absolute copy number
fitted_purity, fitted_ploidy = final_fit[0], final_fit[1]
abs_copy_number_segments = copy.deepcopy(rel_copy_number_segments)
for x in abs_copy_number_segments:
    acn = cnfitter.relative_to_absolute_CN(x[-1], fitted_purity, fitted_ploidy)
    acn_int = round(acn)
    x[-1] = acn
    x.append(acn_int)

# Prepare and write out results
## ranked solutions, final fit, and converted segmented absolute copy number

outfile1 = open(f"{outdir}/{prefix}_ranked_solutions.tsv", "w")
for r in solutions_ranked:
    Line = '\t'.join(r) + '\n'
    outfile1.write(Line)
outfile1.close()           

outfile2 = open(f"{outdir}/{prefix}_fitted_purity_ploidy.tsv", "w")
for r in solutions_ranked:
    Line = '\t'.join(r) + '\n'
    outfile2.write(Line)
outfile2.close()       

outfile3 = open(f"{outdir}/{prefix}_segmented_absolute_copy_number.tsv", "w")
for r in solutions_ranked:
    Line = '\t'.join(r) + '\n'
    outfile3.write(Line)
outfile3.close() 

###############
#!!!!! ADD IN FUNCTION FOR MAJOR AND MINOR COPY NUMBER!!!

# Add in code to redefine segment breaks based on rounded acn?
###############







############################################################################################################################
# filter snps with DP<20
# get mean DP for each phasing set
# remove SNPs with AF_0 or AF_1 = 1 or = 0


# data = [1, 2, 2, 3, 3, 3, 4, 4, 10, 20, 20, 30, 30, 30, 40, 40]
# test = [0.3068182, 0.3522727, 0.2873563, 0.2926829, 0.4230769, 0.2676056, 0.3, 0.15625, 0.2459016, 0.3703704, 0.3157895, 0.2258065, 0.2419355, 0.265625, 0.2686567, 0.234375, 0.2537313, 0.2461538, 0.1343284, 0.3015873, 0.2205882, 0.2173913, 0.2878788, 0.203125, 0.2537313, 0.1940299, 0.2461538, 0.2083333, 0.2191781, 0.2266667, 0.2222222, 0.2575758, 0.1666667, 0.2368421, 0.1025641, 0.2463768, 0.2266667, 0.3670886, 0.225, 0.4, 0.1627907, 0.2289157, 0.3085106, 0.2857143, 0.2926829, 0.2, 0.278481, 0.4727273, 0.3214286, 0.2183908, 0.2592593, 0.3544304, 0.3291139, 0.5633803, 0.4025974, 0.3026316, 0.2962963, 0.2191781, 0.3066667, 0.304878, 0.2771084, 0.2597403, 0.3658537, 0.3076923, 0.2643678, 0.3295455, 0.3555556, 0.3222222, 0.3095238, 0.3181818, 0.3209877, 0.28125, 0.2717391, 0.2527473, 0.3870968, 0.2666667, 0.25, 0.255814, 0.2674419, 0.2467532, 0.297619, 0.3037975, 0.3493976, 0.443038, 0.3846154, 0.4558824, 0.375, 0.28, 0.2739726, 0.3, 0.3478261, 0.3484848, 0.3428571, 0.40625, 0.2786885, 0.6896552, 0.3823529, 0.3043478, 0.28125, 0.4067797, 0.6590909, 0.6363636, 0.7126437, 0.6585366, 0.5769231, 0.7183099, 0.7, 0.828125, 0.7540984, 0.6111111, 0.6666667, 0.7741935, 0.7258065, 0.71875, 0.7313433, 0.6875, 0.7462687, 0.6461538, 0.8507463, 0.6984127, 0.7647059, 0.7391304, 0.6818182, 0.796875, 0.7164179, 0.7910448, 0.7230769, 0.7777778, 0.739726, 0.7066667, 0.7638889, 0.7272727, 0.8333333, 0.7368421, 0.8717949, 0.7246377, 0.7733333, 0.6202532, 0.75, 0.5882353, 0.8255814, 0.7710843, 0.6914894, 0.7142857, 0.695122, 0.7866667, 0.7088608, 0.5272727, 0.6547619, 0.7816092, 0.7283951, 0.6455696, 0.6455696, 0.4366197, 0.5844156, 0.6973684, 0.691358, 0.7671233, 0.68, 0.6585366, 0.686747, 0.7012987, 0.6219512, 0.6794872, 0.7241379, 0.6590909, 0.6222222, 0.6666667, 0.6904762, 0.6590909, 0.6666667, 0.7083333, 0.7065217, 0.7032967, 0.5913978, 0.6666667, 0.7291667, 0.7325581, 0.7209302, 0.7402597, 0.702381, 0.6835443, 0.6506024, 0.556962, 0.6025641, 0.5441176, 0.5972222, 0.7066667, 0.7260274, 0.6714286, 0.6521739, 0.6363636, 0.6, 0.59375, 0.6557377, 0.3103448, 0.5735294, 0.6521739, 0.671875, 0.5932203]

# print(cnfitter.is_unimodal(test))
# print(cnfitter.bimodality_coefficient(test))


# print(cnfitter.is_unimodal(data))
# print(cnfitter.bimodality_coefficient(data))