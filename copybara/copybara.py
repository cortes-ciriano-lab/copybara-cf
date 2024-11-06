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

from time import time
from multiprocessing import cpu_count

import pysam

import copybara.bin_generator as bin_generator
import copybara.read_counter as read_counter
import copybara.smooth as smooth
import copybara.segment as segment
import copybara.fit_absolute as fit_absolute
# import copybara.tmad as tmad
import copybara.helper as helper

logo = """
 ▗▄▄▖ ▗▄▖ ▗▄▄▖▗▖  ▗▖▗▄▄▖  ▗▄▖ ▗▄▄▖  ▗▄▖     
▐▌   ▐▌ ▐▌▐▌ ▐▌▝▚▞▘ ▐▌ ▐▌▐▌ ▐▌▐▌ ▐▌▐▌ ▐▌      ▐▌▀▀▘   
▐▌   ▐▌ ▐▌▐▛▀▘  ▐▌  ▐▛▀▚▖▐▛▀▜▌▐▛▀▚▖▐▛▀▜▌▀▀▗▞▀▘▐▛▀▘ 
▝▚▄▄▖▝▚▄▞▘▐▌    ▐▌  ▐▙▄▞▘▐▌ ▐▌▐▌ ▐▌▐▌ ▐▌  ▝▚▄▖▐▌  
"""

"""
Add in matched normal and panel of normal to args define in copybara run function and arugments of parser
Add parameter to overrule and set ploidy if known from tumour
Write code to estimate tMAD score!
Gene mapping for genes of interest and categorise if gain amp loss del present 
Correct tumour purity if purity 1 and ploidy 0 (or if completely flat profile, i.e. no changes present...)
"""

def copybara_main(args):
    """ main function for copy number analysis """
    if not args.sample:
        # set sample name to default
        args.sample = os.path.splitext(os.path.basename(args.bam))[0]
    print(f'Running as sample {args.sample}')
    outdir = helper.check_outdir(args.outdir, args.overwrite, illegal=None)
    # set number of threads to cpu count if none set
    if not args.threads:
        args.threads = cpu_count()

    # check if files are bam or cram (must have indices)
    if args.bam.endswith('bam') and args.nbam.endswith('bam'):
        args.is_cram = False
        aln_files = {
            'tbam': pysam.AlignmentFile(args.bam, "rb"),
            'nbam': pysam.AlignmentFile(args.nbam, "rb")
        }
    elif args.bam.endswith('cram') and args.nbam.endswith('cram'):
        args.is_cram = True
        aln_files = {
            'tbam': pysam.AlignmentFile(args.bam, "rc"),
            'nbam': pysam.AlignmentFile(args.nbam, "rc")
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
    
     # threads
    if not args.threads:
        args.threads = cpu_count()

    # initialize timing
    checkpoints = [time()]
    time_str = []
   
    # 1. generate bins
    bin_annotations_path = bin_generator.generate_bins(outdir, args.sample, args.ref, args.chromosomes, args.cn_binsize, args.blacklist, args.threads)
    helper.time_function("Binned reference genome", checkpoints, time_str)
    # 2. perform read counting across bins
    read_counts_path = read_counter.count_reads(outdir, args.tumour, args.normal, args.sample, bin_annotations_path, args.readcount_mapq, args.blacklisting, args.bl_threshold, args.bases_filter, args.bases_threshold, args.threads)
    helper.time_function("Performed read counting", checkpoints, time_str)
    # smooth the copy number data
    smoothened_cn_path = smooth.smooth_copy_number(outdir, read_counts_path, args.smoothing_level, args.trim)
    helper.time_function("Performed smoothing", checkpoints, time_str)
    # segment the copy number data
    log2r_cn_path = segment.segment_copy_number(outdir, smoothened_cn_path, args.min_segment_size, args.shuffles, args.p_seg, args.p_val, args.quantile, args.threads)
    helper.time_function("Performed CBS", checkpoints, time_str)
    # fit absolute copy number
    fit_absolute.fit_absolute_cn(outdir, log2r_cn_path, args.sample,
        args.min_ploidy, args.max_ploidy, args.ploidy_step, args.min_cellularity, args.max_cellularity, args.cellularity_step, args.cellularity_buffer, args.overrule_cellularity,
        args.distance_function, args.distance_filter_scale_factor, args.distance_precision,
        args.max_proportion_zero, args.min_proportion_close_to_whole_number, args.max_distance_from_whole_number, args.main_cn_step_change,
        args.min_ps_size, args.min_ps_length, args.threads)
    helper.time_function("Fit absolute copy number", checkpoints, time_str)
    # cleanup tmpdir
    # helper.clean_tmpdir(args.tmpdir, outdir)
    helper.time_function("Total time to perform copy number calling", checkpoints, time_str, final=True)


def parse_args(args):
    """ parse arguments - separated into subcommands """
    global_parser = argparse.ArgumentParser(prog="copybara", description="COPYBARA - Copy number analysis of long-read cfDNA sequencing data")
    global_parser.add_argument('--version', action='version', version=f'COPYBARA {helper.__version__}')
    
    # arguments for default copybara run
    global_parser.add_argument('-bam', '--tumour', nargs='?', type=str, required=True, help='BAM file (must have index)')
    control_group = global_parser.add_mutually_exclusive_group()
    control_group.add_argument('-nbam','--normal', nargs='?', type=str, required=False, help='Matched normal BAM file if available (must have index)')
    control_group.add_argument('-pon','--panel_of_normal', nargs='?', type=str, required=False, help='Panel of normal if available')
    
    global_parser.add_argument('--ref', nargs='?', type=str, required=True, help='Full path to reference genome')
    global_parser.add_argument('--ref_index', nargs='?', type=str, required=False, help='Full path to reference genome fasta index (ref path + ".fai" by default)')
    global_parser.add_argument('--contigs', nargs='?', type=str, help="Contigs/chromosomes to consider (optional, default=All)")

    global_parser.add_argument('--mapq', nargs='?', type=int, default=5, help='Minimum MAPQ to consider a read counting (default=5)')

    global_parser.add_argument('--threads', nargs='?', type=int, const=0, help='Number of threads to use (default=max)')
    global_parser.add_argument('--outdir', nargs='?', required=True, help='Output directory (can exist but must be empty)')
    global_parser.add_argument('--overwrite', action='store_true', help='Use this flag to write to output directory even if files are present')
    
    global_parser.add_argument('--sample', nargs='?', type=str, help='Name to prepend to output files (default=tumour BAM filename without extension)')
    
    global_parser.add_argument('-w', '--cn_binsize', type=int, default=10, help='Bin window size in kbp', required=False)
    global_parser.add_argument('-b', '--blacklist', type=str, help='Path to the blacklist file', required=False)
    global_parser.add_argument('-c', '--chromosomes', nargs='+', default='all', help='Contigs/chromosomes to consider. (optional, default=all). To run only a subset of chromosomes, specify the chromosome numbers separated by spaces. For x and y chromosomes, use 23 and 24, respectively.  E.g. use "-c 1 4 23 24" to run chromosomes 1, 4, X and Y', required=False)
    global_parser.add_argument('--no_blacklist', dest='blacklisting', action='store_false')
    global_parser.set_defaults(blacklisting=True)
    global_parser.add_argument('-blt', '--bl_threshold', type=int,  default='5', help='Percentage overlap between bin and blacklist threshold to tolerate for read counting (default = 0, i.e. no overlap tolerated). Please specify percentage threshold as integer, e.g. "-t 5" ', required=False)
    global_parser.add_argument('--no_basesfilter', dest='bases_filter', action='store_false')
    global_parser.set_defaults(bases_filter=True)
    global_parser.add_argument('-bt', '--bases_threshold', type=int,  default='75', help='Percentage of known bases per bin required for read counting (default = 0, i.e. no filtering). Please specify percentage threshold as integer, e.g. "-bt 95" ', required=False)
    global_parser.add_argument('-sl', '--smoothing_level', type=int, default='10', help='Size of neighbourhood for smoothing.', required=False)
    global_parser.add_argument('-tr', '--trim', type=float, default='0.025', help='Trimming percentage to be used.', required=False)
    global_parser.add_argument('-ms', '--min_segment_size', type=int,  default=5, help='Minimum size for a segement to be considered a segment (default = 5).', required=False)
    global_parser.add_argument('-sf', '--shuffles', type=int,  default=1000, help='Number of permutations (shuffles) to be performed during CBS (default = 1000).', required=False)
    global_parser.add_argument('-ps', '--p_seg', type=float,  default=0.05, help='p-value used to test segmentation statistic for a given interval during CBS using (shuffles) number of permutations (default = 0.05).', required=False)
    global_parser.add_argument('-pv', '--p_val', type=float,  default=0.01, help='p-value used to test validity of candidate segments from CBS using (shuffles) number of permutations (default = 0.01).', required=False)
    global_parser.add_argument('-qt', '--quantile', type=float,  default=0.2, help='Quantile of changepoint (absolute median differences across all segments) used to estimate threshold for segment merging (default = 0.2; set to 0 to avoid segment merging).', required=False)
    global_parser.add_argument('--min_ploidy', type=float, default=1.5, help='Minimum ploidy to be considered for copy number fitting.', required=False)
    global_parser.add_argument('--max_ploidy', type=float, default=5, help='Maximum ploidy to be considered for copy number fitting.', required=False)
    global_parser.add_argument('--ploidy_step', type=float, default=0.01, help='Ploidy step size for grid search space used during for copy number fitting.', required=False)
    
    global_parser.add_argument('--set_ploidy', type=float, default=None, help='Set to sample`s ploidy if known.', required=False)   
    
    global_parser.add_argument('--min_cellularity', type=float, default=0.2, help='Minimum cellularity to be considered for copy number fitting. If hetSNPs allele counts are provided, this is estimated during copy number fitting. Alternatively, a purity value can be provided if the purity of the sample is already known.', required=False)
    global_parser.add_argument('--max_cellularity', type=float, default=1, help='Maximum cellularity to be considered for copy number fitting. If hetSNPs allele counts are provided, this is estimated during copy number fitting. Alternatively, a purity value can be provided if the purity of the sample is already known.', required=False)
    global_parser.add_argument('--cellularity_step', type=float, default=0.01, help='Cellularity step size for grid search space used during for copy number fitting.', required=False)
    global_parser.add_argument('--cellularity_buffer', type=float, default=0.1, help='Cellularity buffer to define purity grid search space during copy number fitting (default = 0.1).', required=False)
    global_parser.add_argument('--overrule_cellularity', type=float, default=None, help='Set to sample`s purity if known. This value will overrule the cellularity estimated using hetSNP allele counts (not used by default).', required=False)   
    global_parser.add_argument('--distance_function', type=str, default='RMSD', help='Distance function to be used for copy number fitting.', choices=['RMSD', 'MAD'], required=False)
    global_parser.add_argument('--distance_filter_scale_factor', type=float, default=1.25, help='Distance filter scale factor to only include solutions with distances < scale factor * min(distance).', required=False)
    global_parser.add_argument('--distance_precision', type=int, default=3, help='Number of digits to round distance functions to', required=False)
    global_parser.add_argument('--max_proportion_zero', type=float, default=0.1, help='Maximum proportion of fitted copy numbers to be tolerated in the zero or negative copy number state.', required=False)
    global_parser.add_argument('--min_proportion_close_to_whole_number', type=float, default=0.5, help='Minimum proportion of fitted copy numbers sufficiently close to whole number to be tolerated for a given fit.', required=False)
    global_parser.add_argument('--max_distance_from_whole_number', type=float, default=0.25, help='Distance from whole number for fitted value to be considered sufficiently close to nearest copy number integer.', required=False)
    global_parser.add_argument('--main_cn_step_change', type=int, default=1, help='Max main copy number step change across genome to be considered for a given solution.', required=False)
    global_parser.add_argument('--min_ps_size', type=int, default=10, help='Minimum size (number of SNPs) for phaseset to be considered for purity estimation.', required=False)
    global_parser.add_argument('--min_ps_length', type=int, default=500000, help='Minimum length (bps) for phaseset to be considered for purity estimation.', required=False)
    global_parser.set_defaults(func=copybara_main)
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
