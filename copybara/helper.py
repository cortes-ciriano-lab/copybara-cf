"""
Module containing misc. useful functions for COPYBARA
Created: 23/10/2024
Python 3.9.7
Carolin Sauer - Adapted from Hillary Elrick
"""

import re
import os
import csv
import sys

from time import time
from datetime import datetime
import argparse

__version__ = "1.0.2_focal_dev"

# def get_contigs(contig_file, ref_index):
# 	""" use the contigs file to return contigs and lengths - otherwise use index """
# 	if contig_file:
# 		with open(contig_file, encoding="utf-8") as f:
# 			contigs = f.readlines()
# 			contigs = [contig.rstrip() for contig in contigs]
# 			return contigs
# 	# otherwise, use the fai to get the contig names
# 	contigs = []
# 	with open(ref_index, encoding="utf-8") as f:
# 		tab_reader = csv.reader(f, delimiter='\t')
# 		for line in tab_reader:
# 			contig = line[0]
# 			contigs.append(contig)
# 	return contigs

# def get_contig_lengths(ref_index):
# 	""" get the contig lengths from the reference """
# 	contig_lengths = {}
# 	with open(ref_index, encoding="utf-8") as f:
# 		tab_reader = csv.reader(f, delimiter='\t')
# 		for line in tab_reader:
# 			contig = line[0]
# 			length = line[1]
# 			contig_lengths[contig] = int(length)
# 	return contig_lengths

# adapted from Hillary Elrick
def time_function(desc, checkpoints, time_str, final=False):
	""" prints the number of seconds elapsed compared to previous checkpoint """
	checkpoints.append(time())
	if not final:
		formatted_time = f'{desc:<60}{round(checkpoints[-1] - checkpoints[-2], 3)} seconds'
	else:
		formatted_time = f'{desc:<60}{round(checkpoints[-1] - checkpoints[0], 3)} seconds\n'
	time_str.append(formatted_time)
	print(formatted_time)
	return

def check_outdir(args_outdir, args_overwrite, illegal=None):
	# create output dir if it doesn't exist
	outdir = os.path.join(os.getcwd(), args_outdir)
	if not os.path.exists(outdir):
		print(f'Creating directory {outdir} to store results')
		os.mkdir(outdir)
	if args_overwrite:
		# don't check for files, overwrite them
		return outdir
	if not illegal:
		# throw error if ANY files present
		if os.listdir(outdir):
			sys.exit(f'Output directory "{outdir}" already exists and contains files. Please remove the files or supply a different directory name.')
	else:
		# check if illegal files exist in outdir
		for f in os.listdir(outdir):
			if f.endswith(illegal):
				sys.exit(f'Output directory "{outdir}" already exists and contains {illegal} files which may be overwritten. Please remove the files or supply a different directory name.')

	return outdir


if __name__ == "__main__":
	print("Helper functions for SAVANA")
