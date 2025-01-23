"""
Script to fit absolute copy number
Created: 27/06/2023
Python 3.9.7
Carolin Sauer
"""

from multiprocessing import Pool
from multiprocessing import cpu_count
import numpy as np
from scipy.stats import gaussian_kde
import statistics
import copy
import sys


import copybara.cn_functions as cnfitter

def process_input_for_purity_estimation(log2r_cn_path):
    rel_cn = []
    with open(log2r_cn_path, "r") as file:
        next(file)
        for line in file:
            fields = line.strip().split("\t")
            chrom,start,end,cn,seg_id,seg_cn=fields[1],int(fields[2]),int(fields[3]),float(fields[-3]),fields[-2],float(fields[-1])
            rel_cn.append([chrom,start,end,cn,seg_id,seg_cn])
    return rel_cn


def define_purity_search_space(rel_cn,bc_thres,dens_thres,min_copy_number,max_copy_number,lower_threshold):
    chr_names = list(dict.fromkeys([x[0] for x in rel_cn]))
    bc_out = []
    chr_count = 0
    for CHROM in chr_names:
        if CHROM == 'chrX' or CHROM == 'chrY':
            continue
        # print(CHROM)
        chr_rel_cn = [x[-3] for x in rel_cn if x[0] == CHROM]
        # apply min and max copy number thresholds for purity estimation
        if min_copy_number is not None:
            assert isinstance(min_copy_number, (int, float)) and np.isscalar(min_copy_number)
            chr_rel_cn = [cn for cn in chr_rel_cn if cn >= min_copy_number]
        if max_copy_number is not None:
            assert isinstance(max_copy_number, (int, float)) and np.isscalar(max_copy_number)
            chr_rel_cn = [cn for cn in chr_rel_cn if cn <= max_copy_number]
        number_chr_segs=len(list(dict.fromkeys([x[-2] for x in rel_cn if x[0] == CHROM])))
        if number_chr_segs >= 2:
            try:
                # modes_result = Modes(chr_rel_cn)
                # print(CHROM, is_unimodal(chr_rel_cn))
                if  cnfitter.is_unimodal(chr_rel_cn) == False:
                    print(CHROM)
                    chr_bc = cnfitter.bimodality_coefficient(chr_rel_cn)
                    # print(CHROM, chr_bc, cnfitter.is_unimodal(chr_rel_cn))
                    if chr_bc >= bc_thres:
                        bc_out.append(chr_rel_cn)
                        chr_count += 1
            except:
                continue
        else:
            continue
    ####### Building site #######
    if len(bc_out) == 0 or chr_count == 1:
        # seg_names = list(dict.fromkeys([x[-2] for x in rel_cn]))
        # bc_segs = []
        # for i in range(len(seg_names)):
        #     s1,s2 = seg_names[i], seg_names[min(i+1,len(seg_names)-1)]
        #     # print(f"testing {s1}_{s2}")
        #     segs_rel_cn = [x[-3] for x in rel_cn if x[-2] == s1 or x[-2] == s2]
        #     if min_copy_number is not None:
        #         assert isinstance(min_copy_number, (int, float)) and np.isscalar(min_copy_number)
        #         segs_rel_cn = [cn for cn in segs_rel_cn if cn >= min_copy_number]
        #     if max_copy_number is not None:
        #         assert isinstance(max_copy_number, (int, float)) and np.isscalar(max_copy_number)
        #         segs_rel_cn = [cn for cn in segs_rel_cn if cn <= max_copy_number]
        #     try:
        #         # modes_result = Modes(chr_rel_cn)
        #         # print(CHROM, is_unimodal(chr_rel_cn))
        #         if cnfitter.is_unimodal(segs_rel_cn) == False:
        #             print(f'{s1}-{s2}')
        #             seg_bc = cnfitter.bimodality_coefficient(segs_rel_cn)
        #             # print(CHROM, chr_bc, cnfitter.is_unimodal(chr_rel_cn))
        #             if seg_bc >= bc_thres:
        #                 # bc_out.append(segs_rel_cn)
        #                 for s in s1,s2:
        #                     if s not in bc_segs:
        #                         bc_segs.append(s) 
        #     except:
        #         continue
        # bc_out = []
        # for S in bc_segs:
        #     out = [x[-3] for x in rel_cn if x[-2] == S]
        #     bc_out.append(out)
        ##########
        # go by half chromosomes
        bc_out = [] # empty bc_out from previous
        splits = []
        for CHROM in chr_names:
            c_min,c_max = min([x[1] for x in rel_cn if x[0] == CHROM]),max([x[2] for x in rel_cn if x[0] == CHROM])
            med=int((c_min+c_max-1)/2)
            splits.append([CHROM, med])
        for i in range(len(splits)):
            if i+1 < len(splits):
                # print(f'{splits[i][0]}-{splits[i+1][0]}')
                split_rel_cn = [x[-3] for x in rel_cn if (x[0] == splits[i][0] and x[1] > splits[i][1]) or (x[0] == splits[i+1][0] and x[1] <= splits[i+1][1])]
                if min_copy_number is not None:
                    assert isinstance(min_copy_number, (int, float)) and np.isscalar(min_copy_number)
                    split_rel_cn = [cn for cn in split_rel_cn if cn >= min_copy_number]
                if max_copy_number is not None:
                    assert isinstance(max_copy_number, (int, float)) and np.isscalar(max_copy_number)
                    split_rel_cn = [cn for cn in split_rel_cn if cn <= max_copy_number]
                try:
                    if cnfitter.is_unimodal(split_rel_cn) == False:
                        print(f'{splits[i][0]}-{splits[i+1][0]}')
                        split_bc = cnfitter.bimodality_coefficient(split_rel_cn)
                        # print(CHROM, chr_bc, cnfitter.is_unimodal(chr_rel_cn))
                        if split_bc >= bc_thres:
                            bc_out.append(split_rel_cn)
                except:
                    continue
        ##########
    # based on bc_out estimate purity centre... 
    bc_out = [x for xs in bc_out for x in xs]
    # test if all chromosomes unimodal and define out:
    if len(bc_out) == 0:
        print("All chromosomes show unimodal distribution. Minimum purity = 0; maximum purity = 0.1")
        pur_centre = 0
    elif len(bc_out) != 0:
        # if multimodal chromosomes identified:
        # look at distribution and find maxima... using python implementation of density default function from R
        dens_x,dens_y = cnfitter.r_density_default(bc_out, n=512)
        # Filter density by min/max copy number thresholds
        filtered_density = [(x, y) for x, y in zip(dens_x,dens_y) if (min_copy_number is None or x >= min_copy_number) and (max_copy_number is None or x <= max_copy_number)]
        # Finding maxima
        maxima = []
        for i in range(1, len(filtered_density) - 1):
            if filtered_density[i][1] > filtered_density[i - 1][1] and filtered_density[i][1] > filtered_density[i + 1][1]:
                maxima.append(filtered_density[i])
        # Filter maxima by lower_threshold
        maxima = [(x, d) for x, d in maxima if d >= lower_threshold * max(dens_y)]
        # Select maxima with density > 0.2 and sort by density
        maxima_select = sorted([m for m in maxima if m[1] > dens_thres], key=lambda m: -m[1])
        # # PLOTTING
        # copy_number_array = np.array(bc_out)
        # plt.figure(figsize=(8, 5))
        # plt.hist(copy_number_array, bins=30, density=True, alpha=0.5, label="Histogram")
        # plt.plot(dens_x, dens_y, label="Density", color='blue')
        # for x, _ in maxima_select:
        #     plt.axvline(x=x, color='red', linestyle='--', label=f"Maxima at {x:.2f}")
        # plt.xlim(min_copy_number, max_copy_number)
        # plt.legend()
        # plt.show()
        # Compute pur_centre and estimate min and max purity search space
        if len(maxima_select) == 1:
            print("unimodal distribution. Minimum purity = 0; maximum purity = 0.1")
            pur_centre = 0
        if len(maxima_select) == 2:
            pur_centre = abs(maxima_select[0][0] - maxima_select[1][0])
        if len(maxima_select) > 2:
            maxval = maxima_select[0][0]
            try:
                # lowerval = max([m[0] for m in maxima_select[1:] if m[0] < maxval])
                lowerval = statistics.mean([m[0] for m in maxima_select[1:] if m[0] < maxval])
            except:
                lowerval = None
            try:
                # upperval = min([m[0] for m in maxima_select[1:] if m[0] > maxval])
                upperval = statistics.mean([m[0] for m in maxima_select[1:] if m[0] > maxval])
            except:
                upperval = None
            if upperval == None and lowerval != None:
                pur_centre = abs(maxval-lowerval)
            elif upperval != None and lowerval == None:
                pur_centre = abs(maxval-upperval)
            elif upperval != None and lowerval != None:
                pur_centre = max([abs(maxval-lowerval),abs(maxval-upperval)])
    pur_centre = 1 if pur_centre >= 1 else pur_centre
    minp = round(max(pur_centre-0.15,0),2)
    min_purity = 0 if minp <= 0.1 else minp
    max_purity = round(min(pur_centre+0.15,1),2) if pur_centre >= 0.15 else round(min(pur_centre+pur_centre,1),2)
    max_purity = 0.1 if pur_centre == 0 else max_purity
    return(pur_centre,min_purity,max_purity)


def process_log2r_input(log2r_cn_path):
    rel_copy_number = []
    with open(log2r_cn_path, "r") as file:
        next(file)
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

def fit_absolute_cn(outdir, log2r_cn_path, sample,
    bc_thres,dens_thres,min_copy_number,max_copy_number,lower_threshold,
    min_ploidy, max_ploidy, ploidy_step, 
    min_cellularity, max_cellularity, cellularity_step,
    distance_function, distance_filter_scale_factor, distance_precision,
    max_proportion_zero, min_proportion_close_to_whole_number, max_distance_from_whole_number, main_cn_step_change,
    threads):
    '''
    Fit absolute copy number and estimate purity and ploidy fit.
    '''
    # check and define threads
    new_threads = min(threads, cpu_count())
    print(f"... CN fitting will use threads = {new_threads}. (threads = {threads} defined; threads = {cpu_count()} available) ...")
    threads = new_threads

    #----
    # 3. Determine purity centre and search space based on chromosomes with multimodal CN data
    #----
    rel_cn = process_input_for_purity_estimation(log2r_cn_path)
    cellularity,min_cellularity,max_cellularity=define_purity_search_space(rel_cn,bc_thres,dens_thres,min_copy_number,max_copy_number,lower_threshold)
    print(cellularity,min_cellularity,max_cellularity)
    print(f"        Purity centre for copy number fitting = {cellularity}. Purity search space: {min_cellularity}-{max_cellularity}.")

    #----
    # 4. Estimate ploidy and fit ACN using estimated sample purity
    ### Method/principles based on rascal R package: https://www.biorxiv.org/content/10.1101/2021.07.19.452658v1
    #----
    rel_copy_number_segments = process_log2r_input(log2r_cn_path)

    # Prepare input for copy number fitting
    relative_CN = [x[-1] for x in rel_copy_number_segments]
    weights = [x[-2] for x in rel_copy_number_segments]

    # Copy number fitting
    fits = cnfitter.estimate_grid_distances(min_cellularity, max_cellularity, cellularity_step, min_ploidy, max_ploidy, ploidy_step, relative_CN, weights=weights, distance_function=distance_function)
    fits_r = cnfitter.reduce_grid(fits,distance_filter_scale_factor = distance_filter_scale_factor)
    solutions = cnfitter.viable_solutions(fits_r, relative_CN, weights,
                        max_proportion_zero = max_proportion_zero,
                        min_proportion_close_to_whole_number = min_proportion_close_to_whole_number,
                        max_distance_from_whole_number = max_distance_from_whole_number, main_cn_step_change = main_cn_step_change)

    # check if viable solutions were found. If not, terminate script and write out error message and arguments to file for inspection and adjustment
    if len(solutions) == 0:
        print("No fits found. See No_fit_found_PARAMS_out.tsv in output")
        with open(f"{outdir}/No_fit_found_PARAMS.tsv", 'w') as params_out:
            params_out.write(f'No viable solution fullfilling set paramaters was found. \nPerform QC, review parameters and rerun if required/appropriate with adjusted parameters.\n')
            # TODO: add this back in (either manually or one step up)
            """
            params_out.write(f'\nPARAMETERS:\n')
            for key, value in vars(args).items():
                    params_out.write(f'{key}: {value}\n')
            """
            params_out.write(f'\nCandidate fits found prior to checking viability/acceptability using set parameters:\n')
            header=['purity','ploidy','x','y','distance']
            params_out.write('\t'.join(header)+'\n')
            for r in fits_r:
                Line = '\t'.join(str(e) for e in r) + '\n'
                params_out.write(Line)
        sys.exit(1) # Exit the script with a status code of 1 (indicating an error)

    solutions_ranked = cnfitter.rank_solutions(solutions,distance_precision)
    final_fit = solutions_ranked[0]
    print(f"        Data fitted to purity = {final_fit[0]} and ploidy = {final_fit[1]}.")


    # Convert relative to absolute copy number and obtain major and minor copy number values if allele counts provided
    # if allele_counts_bed_path == None:
        # print("     ... No phased hetSNPs allele counts provided. Only total absolute copy number will be estimated. Consider providing hetSNPs allele counts if possible. See documentation for instructions on how to generate these.")
    fitted_purity, fitted_ploidy = final_fit[0], final_fit[1]
    abs_copy_number_segments = copy.deepcopy(rel_copy_number_segments)
    for x in abs_copy_number_segments:
        acn = round(cnfitter.relative_to_absolute_CN(x[-1], fitted_purity, fitted_ploidy),4)
        x[-1] = acn if acn > 0 else 0
    # set negative values to 0
    # for x in abs_copy_number_segments:
    #     if x[-1] < 0:
    #         x[-1] = 0

    #----
    # 5. Prepare and write out results
    ### ranked solutions, final fit, and converted segmented absolute copy number
    #----
    outfile1 = open(f"{outdir}/{sample}_ranked_solutions.tsv", "w")
    header=['purity','ploidy','distance','rank']
    outfile1.write('\t'.join(header)+'\n')
    for r in solutions_ranked:
        Line = '\t'.join(str(e) for e in r) + '\n'
        outfile1.write(Line)
    outfile1.close()

    outfile2 = open(f"{outdir}/{sample}_fitted_purity_ploidy.tsv", "w")
    header=['purity','ploidy','distance','rank', 'purity_centre', 'min_purity', 'max_purity']
    outfile2.write('\t'.join(header)+'\n')
    Line = '\t'.join(str(e) for e in final_fit) + f'\t{cellularity}' + f'\t{min_cellularity}' + f'\t{max_cellularity}' + '\n'
    outfile2.write(Line)
    outfile2.close()

    outfile3 = open(f"{outdir}/{sample}_segmented_absolute_copy_number.tsv", "w")
    # if allele_counts_bed_path == None:
    header=['chromosome','start','end','segment_id', 'bin_count', 'sum_of_bin_lengths', 'weight', 'copyNumber']
    # elif allele_counts_bed_path != None:
    #     header=['chromosome','start','end','segment_id', 'bin_count', 'sum_of_bin_lengths', 'weight', 'copyNumber', 'minorAlleleCopyNumber', 'meanBAF', 'no_hetSNPs']
    outfile3.write('\t'.join(header)+'\n')
    for r in abs_copy_number_segments:
        Line = '\t'.join(str(e) for e in r) + '\n'
        outfile3.write(Line)
    outfile3.close()

    # TODO: add this back in - either manually or one step up
    """
    with open(f"{outdir}/PARAMS_out.tsv", 'w') as params_out:
            params_out.write('PARAMETERS:\n')
            for key, value in vars(args).items():
                    params_out.write(f'{key}: {value}\n')
    """

