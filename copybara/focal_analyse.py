"""
Script to analyse focal log2 read counts for putative ecDNA detection
Created: 24/08/2025
Python 3.9.7
Carolin Sauer
"""

import os
from scipy import stats
import numpy as np

import copybara.cn_functions as cnfitter

def process_read_counts(read_counts_path, roi_annot, blacklisting):
    # filename for later output
    # file_name = os.path.split(read_counts_path)[1].removesuffix('.tsv')
    # process input data
    in_data = []
    with open(read_counts_path, "r") as file:
        header_line=next(file)
        for line in file:
            fields = line.strip().split("\t")
            fields[-2] = roi_annot[1] if fields[0] == roi_annot[0] else 'background'
            in_data.append(fields)
    # write out and replace previous file
    outfile = open(read_counts_path, "w")
    if blacklisting == True:
        header=['bin', 'chromosome','start','end', 'gc_content', 'known_bases', 'overlap_blacklist', 'region', 'log2r_copynumber']
    else: 
        header=['bin', 'chromosome','start','end', 'gc_content', 'known_bases', 'region', 'log2r_copynumber']
    outfile.write('\t'.join(header)+'\n')
    for r in in_data:
        Line = '\t'.join(r) + '\n'
        outfile.write(Line)
    outfile.close()
    return in_data

def process_rois(roi):
    # function to read in regions of interest for focal amplification/ecDNA analysis
    bin, roi_name = f"{roi[0]}:{roi[1]}_{roi[2]}", roi[3]
    roi_annot = [bin, roi_name]
    return roi_annot, roi_name

def flatten(xss):
    return [x for xs in xss for x in xs]

def roi_prediction_test(roi_cn, bg_cn, lower_threshold, dens_thres, p_thres = 0.05 ,alternative="greater"):
    bg = np.asarray(bg_cn, dtype=float)
    n_bg = bg.size
    mean_bg = bg.mean()
    sd = bg.std(ddof=1)
    T = (roi_cn - mean_bg) / (sd * np.sqrt(1 + 1/n_bg))
    # p-values
    if alternative == "two-sided":
        p = 2 * stats.t.sf(abs(T), df=n_bg-1)
    elif alternative == "greater":   # enrichment
        p = stats.t.sf(T, df=n_bg-1)
    else:  # "less"
        p = stats.t.cdf(T, df=n_bg-1)
    # estimate CN change separation if p_val is significant
    cnsep = estimate_CN_change_sep(roi_cn, bg_cn, lower_threshold, dens_thres) if p <= p_thres else None
    return dict(mean_bg=mean_bg, sd_bg=sd, df=n_bg-1, T=T, p=p, cnsep=cnsep)

def estimate_CN_change_sep(roi_cn, bg_cn, lower_threshold, dens_thres):
    dens_x,dens_y = cnfitter.r_density_default(bg_cn, n=512)
    density = [(x, y) for x, y in zip(dens_x,dens_y)]
    # Finding maxima
    maxima = []
    for i in range(1, len(density) - 1):
        if density[i][1] > density[i - 1][1] and density[i][1] > density[i + 1][1]:
            maxima.append(density[i])
    # Filter maxima by lower_threshold
    maxima = [(x, d) for x, d in maxima if d >= lower_threshold * max(dens_y)]
    maxima_select = sorted([m for m in maxima if m[1] > dens_thres], key=lambda m: -m[1])
    # estimate cn change sep
    cnsep = abs(maxima_select[0][0] - roi_cn)
    return cnsep

def analyse_focal(outdir, sample, read_counts_path, blacklisting, roi, lower_threshold, dens_thres):
    ''' main function to analyse focal read counts '''
    # extract roi info for annotation
    roi_annot, roi_name = process_rois(roi)
    # process input data and annotate with region name
    in_data = process_read_counts(read_counts_path, roi_annot, blacklisting)
    roi_read_count = flatten([x for x in in_data if x[-2] == roi_name])
    # compute statistics and summary values and estimate cn change separation
    bg_cn = [float(x[-1]) for x in in_data if x[-2] == "background"]
    roi_cn = float(roi_read_count[-1])
    roi_stats = roi_prediction_test(roi_cn, bg_cn, lower_threshold, dens_thres, p_thres=0.05, alternative="greater")
    # prepare main out
    select = [1,2,3,-2,-1]
    out = [roi_read_count[i] for i in select] + [str(x) for x in roi_stats.values()]
    # write output
    outfile = open(f"{outdir}/{sample}_focal_analysis_stats_{roi_name}.tsv", "w")
    header=['chromosome','start','end', 'region', 'log2r_copynumber', 'background_mean', 'background_sd', 'df', 'T', 'p-val', 'cn_change_sep']
    outfile.write('\t'.join(header)+'\n')
    outfile.write('\t'.join(out) + '\n')
    outfile.close()
    return in_data, out
  

#### builing site code dump ####

