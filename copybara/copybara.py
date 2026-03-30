"""
COPYBARA copy number analysis of long-read cfDNA data - main program
Created: 23/10/2024
Python 3.9.7
Carolin Sauer
"""
#!/usr/bin/env python3

import sys
import os
import argparse
import math

from time import time
from multiprocessing import cpu_count

import pysam

import copybara.bin_generator as bin_generator
import copybara.focal_bin_generator as focal_bin_generator
import copybara.read_counter as read_counter
import copybara.focal_analyse as focal_analyse
import copybara.pon_generator as pon_generator
import copybara.smooth as smooth
import copybara.segment as segment
import copybara.fit_absolute as fit_absolute
import copybara.plotting as plotting
# import copybara.tmad as tmad
import copybara.helper as helper

logo = """
 ▗▄▄▖ ▗▄▖ ▗▄▄▖▗▖  ▗▖▗▄▄▖  ▗▄▖ ▗▄▄▖  ▗▄▖     
▐▌   ▐▌ ▐▌▐▌ ▐▌▝▚▞▘ ▐▌ ▐▌▐▌ ▐▌▐▌ ▐▌▐▌ ▐▌      ▐▌▀▀▘   
▐▌   ▐▌ ▐▌▐▛▀▘  ▐▌  ▐▛▀▚▖▐▛▀▜▌▐▛▀▚▖▐▛▀▜▌▀▀▗▞▀▘▐▛▀▘ 
▝▚▄▄▖▝▚▄▞▘▐▌    ▐▌  ▐▙▄▞▘▐▌ ▐▌▐▌ ▐▌▐▌ ▐▌  ▝▚▄▖▐▌  
"""

#####
# COPYBARA pon
#####
def copybara_pon(args):
    """ main function for PoN generation """
    # print(f'blacklisting = {args.blacklisting}\nblacklist: {args.blacklist}')
    # print(f'minploidy={args.min_ploidy}, maxploidy={args.max_ploidy}')
    if not args.pon_name:
        # set sample name to default
        args.pon_name = "PoN"
    print(f'Running as panel of normal generator with prefix {args.pon_name}')
    outdir = helper.check_outdir(args.outdir, args.overwrite, illegal=None)

    # set number of threads to cpu count if none set
    if not args.threads:
        args.threads = cpu_count()

    # confirm pon_list exists
    if not os.path.exists(args.pon_list):
        sys.exit(f'Provided PoN list of bam files: "{args.pon_list}" does not exist. Please provide full path')
    
    # # check if files are bam or cram (must have indices)
    # if args.bam.endswith('bam'): #and args.normal_bam.endswith('bam'):
    #     args.is_cram = False
    #     aln_files = {
    #         'tbam': pysam.AlignmentFile(args.bam, "rb"),
    #         # 'nbam': pysam.AlignmentFile(args.normal_bam, "rb")
    #     }
    # elif args.bam.endswith('cram'): #and args.normal_bam.endswith('cram'):
    #     args.is_cram = True
    #     aln_files = {
    #         'tbam': pysam.AlignmentFile(args.bam, "rc"),
    #         # 'nbam': pysam.AlignmentFile(args.normal_bam, "rc")
    #     }
    # else:
    #     sys.exit('Unrecognized file extension. Input files must be BAM/CRAM.')

    # confirm ref and ref fasta index exist
    if not os.path.exists(args.ref):
        sys.exit(f'Provided reference: "{args.ref}" does not exist. Please provide full path')
    elif args.ref_index and not os.path.exists(args.ref_index):
        sys.exit(f'Provided reference fasta index: "{args.ref_index}" does not exist. Please provide full path')
    elif not os.path.exists(f'{args.ref}.fai'):
        sys.exit(f'Default reference fasta index: "{args.ref}.fai" does not exist. Please provide full path')
    else:
        args.ref_index = f'{args.ref}.fai' if not args.ref_index else args.ref_index
        print(f'Found {args.ref_index} to use as reference fasta index')
    
    # # set 'missing' parameters for bin generator
    args.blacklist = None
    args.blacklist_buffer = None

    # initialize timing
    checkpoints = [time()]
    time_str = []

    # 1. generate bins
    bin_annotations_path = bin_generator.generate_bins(outdir, args.pon_name, args.ref, args.chromosomes, args.cn_binsize, args.blacklist, args.blacklist_buffer, args.threads)
    helper.time_function("Binned reference genome", checkpoints, time_str)
    # 2. perform read counting across bins
    # pon_path = pon_generator.count_reads(outdir, args.pon_list, args.pon_name, bin_annotations_path, args.mapq, args.blacklisting, args.bl_threshold, args.bases_filter, args.bases_threshold, args.threads)
    pon_path = pon_generator.count_reads(outdir, args.pon_list, args.pon_name, bin_annotations_path, args.mapq, args.threads)
    print(f'PoN generated: {pon_path}')
    helper.time_function("Performed read counting", checkpoints, time_str)
    helper.time_function("Total time to generate panel of normals", checkpoints, time_str, final=True)

#####
# COPYBARA focal
#####
def copybara_focal(args):
    """ main function for focal analysis and ecDNA detection in regions of interest (ROIs) """
    # define function to process roi input
    def process_rois(roi_path):
    # function to read in regions of interest for focal amplification/ecDNA analysis
        roi_list = []
        with open(roi_path, "r") as file:
            for line in file:
                fields = line.strip().split("\t")
                roi_length = int(fields[2])-int(fields[1])+1
                roi_list.append([fields[0],int(fields[1]),int(fields[2]),fields[3], roi_length])
        return roi_list
    
    def output_roi_summary(roi_summary, outdir, sample, cnfit):
        outfile = open(f"{outdir}/{sample}_focal_analysis_summary_stats.tsv", "w")
        if cnfit != None:
            header=['chromosome','start','end', 'region', 'log2r_copynumber', 'background_mean', 'background_sd', 'df', 'T', 'p-val', 'cn_change_sep', 'copyNumber', 'category']
        else:
            header=['chromosome','start','end', 'region', 'log2r_copynumber', 'background_mean', 'background_sd', 'df', 'T', 'p-val', 'cn_change_sep']
        outfile.write('\t'.join(header)+'\n')
        for r in roi_summary:
            Line = '\t'.join(r) + '\n'
            outfile.write(Line)
        outfile.close()

    def output_focal_out(outdir, sample, focal_out, roi_name, cnfit):
        outfile = open(f"{outdir}/{sample}_focal_analysis_stats_{roi_name}.tsv", "w")
        if cnfit != None:
            header=['chromosome','start','end', 'region', 'log2r_copynumber', 'background_mean', 'background_sd', 'df', 'T', 'p-val', 'cn_change_sep', 'copyNumber', 'category']
        else:
            header=['chromosome','start','end', 'region', 'log2r_copynumber', 'background_mean', 'background_sd', 'df', 'T', 'p-val', 'cn_change_sep']
        outfile.write('\t'.join(header)+'\n')
        outfile.write('\t'.join(focal_out) + '\n')
        outfile.close()
    
    # check arguments/parameters and define all required input
    print(f'blacklisting = {args.blacklisting}\nblacklist: {args.blacklist}')
    if not args.sample:
        # set sample name to default
        args.sample = os.path.splitext(os.path.basename(args.bam))[0]
    print(f'Running as sample {args.sample}')
    outdir = helper.check_outdir(args.outdir, args.overwrite, illegal=None)

    # set number of threads to cpu count if none set
    if not args.threads:
        args.threads = cpu_count()

    # set blacklisting parameter if not provided
    if not args.blacklist:
        args.blacklisting = False

    # check if files are bam or cram (must have indices)
    if args.bam.endswith('bam'): #and args.normal_bam.endswith('bam'):
        args.is_cram = False
        aln_files = {
            'tbam': pysam.AlignmentFile(args.bam, "rb"),
        }
    elif args.bam.endswith('cram'): #and args.normal_bam.endswith('cram'):
        args.is_cram = True
        aln_files = {
            'tbam': pysam.AlignmentFile(args.bam, "rc"),
        }
    else:
        sys.exit('Unrecognized file extension. Input files must be BAM/CRAM.')

    # confirm ref and ref fasta index exist
    if not os.path.exists(args.ref):
        sys.exit(f'Provided reference: "{args.ref}" does not exist. Please provide full path')
    elif args.ref_index and not os.path.exists(args.ref_index):
        sys.exit(f'Provided reference fasta index: "{args.ref_index}" does not exist. Please provide full path')
    elif not os.path.exists(f'{args.ref}.fai'):
        sys.exit(f'Default reference fasta index: "{args.ref}.fai" does not exist. Please provide full path')
    else:
        args.ref_index = f'{args.ref}.fai' if not args.ref_index else args.ref_index
        print(f'Found {args.ref_index} to use as reference fasta index')

    # initialize timing
    checkpoints = [time()]
    time_str = []
    ###

    # go through rois 
    roi_list = process_rois(args.roi)
    roi_summary = []
    for roi in roi_list:
        roi_name=roi[3]
        print(f'*** Focal analysis for region: {roi_name} ***')
        # outdir_reg=f"{outdir}/{roi_name}"
        outdir_reg=helper.check_outdir(f"{outdir}/{roi_name}", args.overwrite, illegal=None)
        # 1. generate annotated bin bed file 
        bin_annotations_path = focal_bin_generator.generate_bins(outdir_reg, args.sample, args.ref, roi, args.roi_buffer, args.n_regions, args.chromosomes, args.blacklist)
        helper.time_function("Binned reference genome", checkpoints, time_str)
        # 2. count reads and normalise
        try:
            read_counts_path,nmode,coverage = read_counter.count_reads(outdir_reg, args.bam, None, None ,args.sample, bin_annotations_path, args.ref, args.mapq, args.blacklisting, args.bl_threshold, args.bases_filter, args.bases_threshold, False, None, None, args.threads)
            helper.time_function("Performed read counting", checkpoints, time_str)
        except: 
            print(f"read counter normalisation and gc correction failed likely due to insufficient coverage. moving on to next region...")
            focal_out = [str(x) for x in roi[0:4]] + ['NA','NA','NA','NA','NA','NA','NA']
            if args.cnfit != None:
                focal_out += ['NA','NA']
            roi_summary.append(focal_out)
            output_focal_out(outdir_reg, args.sample, focal_out, roi_name, args.cnfit)
            continue
        # 3. analyse focal readcounts
        dens_thres=0.2
        lower_threshold=0   
        focal_cn,focal_out,make_plot =focal_analyse.analyse_focal(outdir_reg, args.sample, read_counts_path, args.blacklisting, roi, args.cnfit, args.alternative, args.p_thres, lower_threshold, dens_thres)
        output_focal_out(outdir_reg, args.sample, focal_out, roi_name, args.cnfit)
        roi_summary.append(focal_out)
        # 4. plotting output
        if make_plot == True:
            plotting.plot_focal_results(focal_cn, focal_out, args.sample, outdir_reg)
            helper.time_function("Performed focal analysis and visualisation", checkpoints, time_str)
    # prepare and output all roi out data
    output_roi_summary(roi_summary, outdir, args.sample, args.cnfit)
    helper.time_function("Total time to perform focal analysis", checkpoints, time_str, final=True) 


#####
# COPYBARA (main)
#####
def copybara_main(args):
    """ main function for copy number analysis """
    print(f'blacklisting = {args.blacklisting}\nblacklist: {args.blacklist}')
    # print(f'minploidy={args.min_ploidy}, maxploidy={args.max_ploidy}')
    if not args.sample:
        # set sample name to default
        args.sample = os.path.splitext(os.path.basename(args.bam))[0]
    print(f'Running as sample {args.sample}')
    outdir = helper.check_outdir(args.outdir, args.overwrite, illegal=None)
    
    # set number of threads to cpu count if none set
    if not args.threads:
        args.threads = cpu_count()

    # set blacklisting parameter if not provided
    if not args.blacklist:
        args.blacklisting = False

    # check if files are bam or cram (must have indices)
    if args.bam.endswith('bam'): #and args.normal_bam.endswith('bam'):
        args.is_cram = False
        aln_files = {
            'tbam': pysam.AlignmentFile(args.bam, "rb"),
            # 'nbam': pysam.AlignmentFile(args.normal_bam, "rb")
        }
    elif args.bam.endswith('cram'): #and args.normal_bam.endswith('cram'):
        args.is_cram = True
        aln_files = {
            'tbam': pysam.AlignmentFile(args.bam, "rc"),
            # 'nbam': pysam.AlignmentFile(args.normal_bam, "rc")
        }
    else:
        sys.exit('Unrecognized file extension. Input files must be BAM/CRAM.')

    # confirm ref and ref fasta index exist
    if not os.path.exists(args.ref):
        sys.exit(f'Provided reference: "{args.ref}" does not exist. Please provide full path')
    elif args.ref_index and not os.path.exists(args.ref_index):
        sys.exit(f'Provided reference fasta index: "{args.ref_index}" does not exist. Please provide full path')
    elif not os.path.exists(f'{args.ref}.fai'):
        sys.exit(f'Default reference fasta index: "{args.ref}.fai" does not exist. Please provide full path')
    else:
        args.ref_index = f'{args.ref}.fai' if not args.ref_index else args.ref_index
        print(f'Found {args.ref_index} to use as reference fasta index')

    # if not args.size_select:
    #     args.size_select = False
    if args.size_select == True and (args.min_read_size is None or args.max_read_size is None):
        sys.exit(f'-size_select requires --min_read_size and --max_read_size.')
        if len(args.min_read_size) != len(args.max_read_size):
            sys.exit(f'the length of --min_read_size and --max_read_size must be identical.')
    if args.size_select == True:
        print(f'Performing size selection prior to copy number analysis. (min_read_size: {args.min_read_size} and max_read_size: {args.max_read_size}')

    # set minimum number of bins for segment size based on min_segment_size argument
    bin_size=args.cn_binsize * int(1000)
    args.min_segment_size = math.ceil(args.min_segment_size / bin_size)
    print(f'Min segment size = {args.min_segment_size} bins.')

    # initialize timing
    checkpoints = [time()]
    time_str = []

    # print(f'size select: {args.size_select} | min read size: {args.min_read_size} | {type(args.min_read_size)}')

    # 1. generate bins
    bin_annotations_path = bin_generator.generate_bins(outdir, args.sample, args.ref, args.chromosomes, args.cn_binsize, args.blacklist, args.blacklist_buffer, args.threads)
    helper.time_function("Binned reference genome", checkpoints, time_str)
    # 2. perform read counting across bins
    read_counts_path,nmode,coverage = read_counter.count_reads(outdir, args.bam, args.normal_bam, args.panel_of_normal ,args.sample, bin_annotations_path, args.ref, args.mapq, args.blacklisting, args.bl_threshold, args.bases_filter, args.bases_threshold, args.size_select, args.min_read_size, args.max_read_size, args.threads)
    helper.time_function("Performed read counting", checkpoints, time_str)
    # smooth the copy number data
    smoothened_cn_path = smooth.smooth_copy_number(outdir, read_counts_path, args.smoothing_level, args.trim)
    helper.time_function("Performed smoothing", checkpoints, time_str)
    # segment the copy number data
    log2r_cn_path = segment.segment_copy_number(outdir, smoothened_cn_path, args.min_segment_size, args.shuffles, args.p_seg, args.p_val, args.quantile, args.quantile_low_cov, args.min_coverage, coverage, args.threads)
    helper.time_function("Performed CBS", checkpoints, time_str)
    # fit absolute copy number
    bc_thres=0.25
    dens_thres=0.2
    min_copy_number=-2
    max_copy_number=2
    lower_threshold=0
    absolute_cn_path,cn_fit_path = fit_absolute.fit_absolute_cn(outdir, nmode, log2r_cn_path, args.sample, coverage, args.goi,
        bc_thres,dens_thres,min_copy_number,max_copy_number,lower_threshold,
        args.min_ploidy, args.max_ploidy, args.ploidy_step, args.min_cellularity, args.max_cellularity, args.cellularity_step, args.cellularity_buffer, 
        args.distance_function, args.distance_filter_scale_factor, args.distance_precision,
        args.max_proportion_zero, args.min_proportion_close_to_whole_number, args.max_distance_from_whole_number, args.main_cn_step_change,
        args.threads)
    plotting.plot_copy_number(absolute_cn_path, log2r_cn_path, cn_fit_path, args.sample, args.no_plot_points, outdir)
    helper.time_function("Fit absolute copy number", checkpoints, time_str)
    helper.time_function("Total time to perform copy number calling", checkpoints, time_str, final=True)

####
# bam/cram checker into pon_generator when opening lists! 
# above values into parameters and defaults! 
# improve segmentation and decrease merging?
# fix no_blacklist flag... 
####

def parse_args(args):
    """ parse arguments - separated into subcommands """
    global_parser = argparse.ArgumentParser(prog="copybara", description="COPYBARA - Copy number analysis of long-read cfDNA sequencing data")
    global_parser.add_argument('--version', action='version', version=f'COPYBARA {helper.__version__}')
    subparsers = global_parser.add_subparsers(title="subcommands", help='COPYBARA sub-commands', dest='command')
    subparsers.required = False
    
    ####
    # arguments for copybara pon
    ####
    pon_parser = subparsers.add_parser("pon", help="generate panel of normals from bam files for COPYBARA copy number analysis")
    pon_parser.add_argument('-i','--pon_list', nargs='?', type=str, default=None, required=True, help='Path to text file containing a list of bam files from which to build the panel of normals')
    pon_parser.add_argument('--pon_name', nargs='?', type=str, default='PoN', required=False, help='Name to prepend to PoN output file (default="PoN")')
    pon_parser.add_argument('--ref', nargs='?', type=str, required=True, help='Full path to reference genome')
    pon_parser.add_argument('--ref_index', nargs='?', type=str, required=False, help='Full path to reference genome fasta index (ref path + ".fai" by default)')
    pon_parser.add_argument('--mapq', nargs='?', type=int, default=5, help='Minimum MAPQ to consider a read counting (default=5)')
    pon_parser.add_argument('--threads', nargs='?', type=int, const=0, help='Number of threads to use (default=max)')
    pon_parser.add_argument('--outdir', nargs='?', required=True, help='Output directory (can exist but must be empty)')
    pon_parser.add_argument('--overwrite', action='store_true', help='Use this flag to write to output directory even if files are present')
    pon_parser.add_argument('-w', '--cn_binsize', type=int, default=500, help='Bin window size in kbp (default=500)', required=False)
    pon_parser.add_argument('-c', '--chromosomes', nargs='+', default='all', help='Contigs/chromosomes to consider. (optional, default=all). To run only a subset of chromosomes, specify the chromosome numbers separated by spaces. For x and y chromosomes, use 23 and 24, respectively.  E.g. use "-c 1 4 23 24" to run chromosomes 1, 4, X and Y', required=False)
    pon_parser.set_defaults(func=copybara_pon)

    ####
    # arguments for copybara focal
    ####
    focal_parser = subparsers.add_parser("focal", help="detect putative ecDNA/focal amplifications across regions of interest")
    focal_parser.add_argument('-b', '--bam', nargs='?', type=str, required=True, help='BAM file (must have index)')
    focal_parser.add_argument('--sample', nargs='?', type=str, help='Name to prepend to output files (default=tumour BAM filename without extension)')
    focal_parser.add_argument('--roi', nargs='?', type=str, required=True, help='Full path to region list of interest for focal analysis (must be tsv or bed file in bed file format)')
    focal_parser.add_argument('--roi_buffer', nargs='?', type=int, default=500000, required=False, help='Length of region (in bp) flanking the a given region of interest to be excluded from random background region sampling')
    focal_parser.add_argument('--n_regions', nargs='?', type=int, default = 1000, required=False, help='Number of random background regions to be sampled for computing of ROIs')    
    focal_parser.add_argument('--ref', nargs='?', type=str, required=True, help='Full path to reference genome')
    focal_parser.add_argument('--ref_index', nargs='?', type=str, required=False, help='Full path to reference genome fasta index (ref path + ".fai" by default)')
    focal_parser.add_argument('--cnfit', nargs='?', type=str, required=False, default = None, help='Full path to "fitted_purity_ploidy.tsv" COPYBARA output file')
    focal_parser.add_argument('--alternative', nargs='?', type=str, required=False, default = 'greater', help='Alternative to test ["greater", "two-sided", "less"]')
    focal_parser.add_argument('--p_thres', nargs='?', type=float, required=False, default = 0.05, help='p value threshold to consider significance')
    focal_parser.add_argument('--mapq', nargs='?', type=int, default=5, help='Minimum MAPQ to consider a read counting (default=5)')
    focal_parser.add_argument('--threads', nargs='?', type=int, const=0, help='Number of threads to use (default=max)')
    focal_parser.add_argument('--outdir', nargs='?', required=True, help='Output directory (can exist but must be empty)')
    focal_parser.add_argument('--overwrite', action='store_true', help='Use this flag to write to output directory even if files are present')
    focal_parser.add_argument('-c', '--chromosomes', nargs='+', default='all', help='Contigs/chromosomes to consider. (optional, default=all). To run only a subset of chromosomes, specify the chromosome numbers separated by spaces. For x and y chromosomes, use 23 and 24, respectively.  E.g. use "-c 1 4 23 24" to run chromosomes 1, 4, X and Y', required=False)
    focal_parser.add_argument('-bl', '--blacklist', type=str, help='Path to the blacklist file', required=False)
    focal_parser.add_argument('--no_blacklist', dest='blacklisting', action='store_false')
    focal_parser.set_defaults(blacklisting=True)
    focal_parser.add_argument('-blt', '--bl_threshold', type=int,  default='1', help='Percentage overlap between bin and blacklist threshold to tolerate for read counting (default = 0, i.e. no overlap tolerated). Please specify percentage threshold as integer, e.g. "-t 5" ', required=False)
    focal_parser.add_argument('--no_basesfilter', dest='bases_filter', action='store_false')
    focal_parser.set_defaults(bases_filter=True)
    focal_parser.add_argument('-bt', '--bases_threshold', type=int,  default='75', help='Percentage of known bases per bin required for read counting (default = 0, i.e. no filtering). Please specify percentage threshold as integer, e.g. "-bt 95" ', required=False)
    focal_parser.set_defaults(func=copybara_focal)

    try:
        global_parser.exit_on_error = False
        subparser = global_parser.parse_args().command if not args else global_parser.parse_args(args).command
    except argparse.ArgumentError as _:
        # unable to parse args, set args to None
        subparser = None

    ####
    # arguments for default/main copybara run
    ####
    if not subparser:
        global_parser.add_argument('-b', '--bam', nargs='?', type=str, required=True, help='BAM file (must have index)')
        control_group = global_parser.add_mutually_exclusive_group()
        control_group.add_argument('-nb','--normal_bam', nargs='?', type=str, default=None, required=False, help='Matched normal BAM file if available (must have index)')
        control_group.add_argument('-pon','--panel_of_normal', nargs='?', type=str, default=None, required=False, help='Path to panel of normal read count file')
        global_parser.add_argument('--ref', nargs='?', type=str, required=True, help='Full path to reference genome')
        global_parser.add_argument('--ref_index', nargs='?', type=str, required=False, help='Full path to reference genome fasta index (ref path + ".fai" by default)')
        global_parser.add_argument('--goi', nargs='?', type=str, required=False, help='Full path to gene list of interest (must be tsv or bed file in bed file format)')
        global_parser.add_argument('--mapq', nargs='?', type=int, default=5, help='Minimum MAPQ to consider a read counting (default=5)')
        global_parser.add_argument('--threads', nargs='?', type=int, const=0, help='Number of threads to use (default=max)')
        global_parser.add_argument('--outdir', nargs='?', required=True, help='Output directory (can exist but must be empty)')
        global_parser.add_argument('--overwrite', action='store_true', help='Use this flag to write to output directory even if files are present')
        global_parser.add_argument('--sample', nargs='?', type=str, help='Name to prepend to output files (default=tumour BAM filename without extension)')
        global_parser.add_argument('-w', '--cn_binsize', type=int, default=500, help='Bin window size in kbp', required=False)
        global_parser.add_argument('-c', '--chromosomes', nargs='+', default='all', help='Contigs/chromosomes to consider. (optional, default=all). To run only a subset of chromosomes, specify the chromosome numbers separated by spaces. For x and y chromosomes, use 23 and 24, respectively.  E.g. use "-c 1 4 23 24" to run chromosomes 1, 4, X and Y', required=False)
        global_parser.add_argument('--blacklist', type=str, help='Path to the blacklist file', required=False)
        global_parser.add_argument('--blacklist_buffer', type=int, default=0, help='Length of region (in bp) flanking the blacklisted region to be excluded. (default = 0)', required=False)
        global_parser.add_argument('--no_blacklist', dest='blacklisting', action='store_false')
        global_parser.set_defaults(blacklisting=True)
        global_parser.add_argument('--bl_threshold', type=int,  default='1', help='Percentage overlap between bin and blacklist threshold to tolerate for read counting (default = 0, i.e. no overlap tolerated). Please specify percentage threshold as integer, e.g. "-t 5" ', required=False)
        global_parser.add_argument('--no_basesfilter', dest='bases_filter', action='store_false')
        global_parser.set_defaults(bases_filter=True)
        global_parser.add_argument('--bases_threshold', type=int,  default='75', help='Percentage of known bases per bin required for read counting (default = 0, i.e. no filtering). Please specify percentage threshold as integer, e.g. "-bt 95" ', required=False)
        global_parser.add_argument('--smoothing_level', type=int, default='10', help='Size of neighbourhood for smoothing.', required=False)
        global_parser.add_argument('--trim', type=float, default='0.025', help='Trimming percentage to be used.', required=False)
        global_parser.add_argument('--min_segment_size', type=int,  default=2500000, help='Minimum size for a segement to be considered a segment (default = 2500000 (2.5Mb)).', required=False)
        global_parser.add_argument('--shuffles', type=int,  default=1000, help='Number of permutations (shuffles) to be performed during CBS (default = 1000).', required=False)
        global_parser.add_argument('--p_seg', type=float,  default=0.05, help='p-value used to test segmentation statistic for a given interval during CBS using (shuffles) number of permutations (default = 0.05).', required=False)
        global_parser.add_argument('--p_val', type=float,  default=0.01, help='p-value used to test validity of candidate segments from CBS using (shuffles) number of permutations (default = 0.01).', required=False)
        global_parser.add_argument('--quantile', type=float,  default=0, help='Quantile of changepoints (absolute median differences across all segments) used to estimate threshold for segment merging (default = 0, i.e. no segment merging).', required=False)
        global_parser.add_argument('--quantile_low_cov', type=float,  default=0.2, help='Quantile of changepoints used to estimate threshold for segment merging if sample < min coverage (default = 0.2).', required=False)
        global_parser.add_argument('--min_coverage', type=float,  default=0.01, help='Quantile of changepoint (absolute median differences across all segments) used to estimate threshold for segment merging (default = 0, i.e. no segment merging).', required=False)
        global_parser.add_argument('--min_ploidy', type=float, default=1.7, help='Minimum ploidy to be considered for copy number fitting.', required=False)
        global_parser.add_argument('--max_ploidy', type=float, default=3.7, help='Maximum ploidy to be considered for copy number fitting.', required=False)
        global_parser.add_argument('--ploidy_step', type=float, default=0.1, help='Ploidy step size for grid search space used during for copy number fitting.', required=False)
        # global_parser.add_argument('--set_ploidy', type=float, default=None, help='Set to sample`s ploidy if known.', required=False)   
        # global_parser.add_argument('--ploidy_buffer', type=float, default=0.3, help='Ploidy buffer to define ploidy grid search space during copy number fitting when --set_ploidy is provided (default = 0.3).', required=False)
        global_parser.add_argument('--min_cellularity', type=float, default=0, help='Minimum cellularity to be considered for copy number fitting. If hetSNPs allele counts are provided, this is estimated during copy number fitting. Alternatively, a purity value can be provided if the purity of the sample is already known.', required=False)
        global_parser.add_argument('--max_cellularity', type=float, default=1, help='Maximum cellularity to be considered for copy number fitting. If hetSNPs allele counts are provided, this is estimated during copy number fitting. Alternatively, a purity value can be provided if the purity of the sample is already known.', required=False)
        global_parser.add_argument('--cellularity_step', type=float, default=0.01, help='Cellularity step size for grid search space used during for copy number fitting.', required=False)
        global_parser.add_argument('--cellularity_buffer', type=float, default=0.05, help='Cellularity buffer to define purity grid search space during copy number fitting (default = 0.05).', required=False)
        # global_parser.add_argument('--overrule_cellularity', type=float, default=None, help='Set to sample`s purity if known. This value will overrule the cellularity estimated using hetSNP allele counts (not used by default).', required=False)   
        global_parser.add_argument('--distance_function', type=str, default='RMSD', help='Distance function to be used for copy number fitting.', choices=['RMSD', 'MAD'], required=False)
        global_parser.add_argument('--distance_filter_scale_factor', type=float, default=1.25, help='Distance filter scale factor to only include solutions with distances < scale factor * min(distance).', required=False)
        global_parser.add_argument('--distance_precision', type=int, default=3, help='Number of digits to round distance functions to', required=False)
        global_parser.add_argument('--max_proportion_zero', type=float, default=0.1, help='Maximum proportion of fitted copy numbers to be tolerated in the zero or negative copy number state.', required=False)
        global_parser.add_argument('--min_proportion_close_to_whole_number', type=float, default=0.5, help='Minimum proportion of fitted copy numbers sufficiently close to whole number to be tolerated for a given fit.', required=False)
        global_parser.add_argument('--max_distance_from_whole_number', type=float, default=0.25, help='Distance from whole number for fitted value to be considered sufficiently close to nearest copy number integer.', required=False)
        global_parser.add_argument('--main_cn_step_change', type=int, help='Max main copy number step change across genome to be considered for a given solution.', required=False)
        global_parser.add_argument('--size_select', dest='size_select', action='store_true', help='add --size_select to turn on size_selection.')
        global_parser.set_defaults(size_select=False)
        global_parser.add_argument('--min_read_size', type=int, nargs='+', required=False, help='minimum read size for size selection prior to copy number analysis. Provide multiple values if looking to combine multiple sizes.')
        global_parser.add_argument('--max_read_size', type=int, nargs='+', required=False, help='maximum read size for size selection prior to copy number analysis. Provide multiple values if looking to combine multiple sizes.')
        global_parser.add_argument('--no_plot_points', type=int, default=7500 ,required=False, help='Number of points to plot, e.g. downsampling when small binsize is used for cleaner visuals. (default = 10000).')
        global_parser.set_defaults(func=copybara_main)
        parsed_args = global_parser.parse_args() if not args else global_parser.parse_args(args)
    else:
        global_parser.exit_on_error = True
        parsed_args = global_parser.parse_args() if not args else global_parser.parse_args(args)
        
    return parsed_args

def main(args=None):
    """ main function for COPYBARA - collects command line arguments and executes algorithm """
    print(logo)
    print(f'Version {helper.__version__}')
    src_location = __file__
    print(f'Source: {src_location}\n')
    args = parse_args(args)
    args.func(args)

if __name__ == "__main__":
    main()
