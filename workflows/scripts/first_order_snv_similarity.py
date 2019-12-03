import sys
import argparse
import pandas as pd
import re
import numpy as np

############################################
# Compute first-order similarity between 
# primary-derived pairs based on CNVs
# Inputs:
#   See args
# Outputs:
# mcadosch 8/28/2018
############################################


#######################################################################################################
# 0. Parse command line inputs
#######################################################################################################
parser = argparse.ArgumentParser(description='Compute first-order CNVs similarity between primary-derived pairs')
parser.add_argument('--pair_id', type=str, required=True,
					help='pair_id')
parser.add_argument('--case_sample_id', type=str, required=True,
					help='case sample id in primary-derived pairs')
parser.add_argument('--control_sample_id', type=str, required=True,
					help='control sample id in primary-derived pairs')
parser.add_argument('--case_sample_snvs_path', type=str, required=True,
					help='case sample snvs path')
parser.add_argument('--control_sample_snvs_path', type=str, required=True,
					help='control sample snvs path')

####################################################################################################
# HELPER FUNCTIONS
####################################################################################################
def add_log2_cnvs(sample_seg):
	"""Add log2 of relative coverage
	"""
	sample_seg["segment_mean_log2"] = np.log2(sample_seg["Segment_Mean"])
	return sample_seg

def discretize_cnv_events(sample_seg, amp_thresh=0.3, del_thresh=-0.7):
	"""Discretize AMP and DEL events
	"""
	sample_seg["CNV_event"] = 0
	sample_seg.loc[sample_seg["segment_mean_log2"]>=amp_thresh, "CNV_event"] = 1
	sample_seg.loc[sample_seg["segment_mean_log2"]<=del_thresh, "CNV_event"] = -1
	return sample_seg

def overlap(ival1, ival2, bp_buffer=200):
	"""Returns True if two intervals overlap
	Parmas:
		- ival1, ival2: (list) intervals 
		- bp_buffer: deprecated
	"""
	a = ival1[0]; c = ival1[1]; b = ival2[0]; d = ival2[1]
	e = max(a,b); f = min(c,d)
	if (a <= d and c >= b):
		return (True, [e, f], f-e)
	else:
		return (False, 0, 0)

def main(pair_id, case_sample_id, control_sample_id, sample_1_seg_path, sample_2_seg_path):
	# Chromosome list for iteration
	chrom_list = [str(i) for i in np.arange(1, 23)]

	## Read Seg files
	sample_1_seg = pd.read_table(sample_1_seg_path, comment="#")
	sample_2_seg = pd.read_table(sample_2_seg_path, comment="#")

	## Add log2 copy number
	sample_1_seg = add_log2_cnvs(sample_1_seg)
	sample_2_seg = add_log2_cnvs(sample_2_seg)

	## Discretize copy number events
	sample_1_seg = discretize_cnv_events(sample_1_seg)
	sample_2_seg = discretize_cnv_events(sample_2_seg)

	## Only keep CNV events
	sample_1_cnv_events = sample_1_seg[sample_1_seg["CNV_event"]!=0]
	sample_2_cnv_events = sample_2_seg[sample_2_seg["CNV_event"]!=0]

	## Overlap intervals and overlap interval lengths
	overlap_ivals = []; overlap_ivals_lens = []; overlap_chroms = [];
	for chrom in chrom_list:
		# Sample 1 and 2 CNVs in chromosome @chrom
		sample_1_chrom_cnvs = sample_1_cnv_events[sample_1_cnv_events["Chromosome"]==chrom]
		sample_2_chrom_cnvs = sample_2_cnv_events[sample_2_cnv_events["Chromosome"]==chrom]
		# Iterate through both CNV events in chromosome, and count overlaps
		for idx, cnv1 in sample_1_chrom_cnvs.iterrows():
			for jdx, cnv2 in sample_2_chrom_cnvs.iterrows():
				ival1 = cnv1[["Start", "End"]].tolist()
				ival2 = cnv2[["Start", "End"]].tolist()
				is_overlap, overlap_ival, overlap_len = overlap(ival1, ival2)
				if is_overlap: 
					overlap_chroms.append(chrom)
					overlap_ivals.append(overlap_ival) 
					overlap_ivals_lens.append(overlap_len);

	## Prepare dataframe for saving csv file
	overlaps = pd.DataFrame({"overlap_start": [i[0] for i in overlap_ivals], \
							 "overlap_end": [i[1] for i in overlap_ivals], \
							 "overlap_range": overlap_ivals_lens, \
							 "overlap_chrom": overlap_chroms})

	# Add sample ids
	overlaps["case_sample_id"] = case_sample_id
	overlaps["control_sample_id"] = control_sample_id
							 
	overlaps[["case_sample_id", "control_sample_id", "overlap_chrom", "overlap_start", "overlap_end", "overlap_range"]] \
			.to_csv("%s.cnv_first_order_similarities.txt"%pair_id, sep="\t", index=None)

####################################################################################################
# MAIN
####################################################################################################
args = parser.parse_args()
main(args.pair_id, args.case_sample_id, args.control_sample_id, args.case_sample_snvs_path, args.control_sample_snvs_path)

	

	

