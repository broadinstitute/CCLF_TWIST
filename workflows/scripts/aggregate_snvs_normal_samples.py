import argparse
import pandas as pd
import re

############################################
# Aggregate SNVs from all normal samples in sample set
# Inputs:
#   See args
# Outputs:
# mcadosch 9/13/2017
############################################


#######################################################################################################
# 0. Parse command line inputs
#######################################################################################################
parser = argparse.ArgumentParser(description='Aggregate SNVs from all samples in sample set')
parser.add_argument('--tumor_variant_files', type=str, required=True, nargs='+',
                    help='Paths to files with somatic variant calls')
parser.add_argument('--external_ids', type=str, required=True, nargs='+',
                    help='external_ids')
parser.add_argument('--sample_ids', type=str, required=True, nargs='+',
                    help='sample_ids')
parser.add_argument('--tsca_id', type=str, required=True,
                    help='tsca_id')
parser.add_argument('--sample_types', type=str, required=True, nargs='+', 
					help="sample_types")

####################################################################################################
# 0. Ensure inputs are correct
# Problems arise especially if a sample does not have a file, 
# then the arrays `sample_ids`, `tumor_variant_files` will not correspond to each other.
# As a safety measure, the sample_ids should be part of the filenames
####################################################################################################
def assert_inputs(sample_ids, external_ids, files, sample_types):
	"""Ensure inputs are correct, and that the external_ids, sample_ids, files arrays correspond to each other
	"""
	print( "Number of external_ids: %s / Number of files: %s" %(len(external_ids), len(files)) )
	
	# > Mismatch between number of sample_ids and files: some file must be missing
	if len(external_ids) != len(files):
		clean_sample_ids = []; clean_external_ids = []; clean_files = []; clean_sample_types = [];
		for idx, sid in enumerate(sample_ids):
			# Find file associated with sid (sid must appear in filename)
		    associated_file = [f for f in files if sid in f]
		    # There is file associated with sample id
		    if len(associated_file) > 0: 
			    clean_sample_ids.append(sid)
			    clean_external_ids.append(external_ids[idx])
			    clean_sample_types.append(sample_types[idx])
			    clean_files.append(associated_file[0])			
		return clean_sample_ids, clean_external_ids, clean_files, clean_sample_types
	
	# > Inputs are correct
	return sample_ids, external_ids, files, sample_types

####################################################################################################
# 1. Aggregate somatic SNVs from all tumor samples
####################################################################################################
def aggregate_snvs(external_ids, tumor_variant_files):
	
	dfs = []
	for external_id, tumor_variants_file in zip(external_ids, tumor_variant_files):
		df = pd.read_table(tumor_variants_file)
		df['external_id'] = external_id
		dfs.append(df)

	aggregate_snvs = pd.concat(dfs, axis=0)
	aggregate_snvs.to_csv("%s.aggregate_somatic_SNVs.txt"%args.tsca_id, sep="\t", index=None)

	return


####################################################################################################
# MAIN
####################################################################################################
args = parser.parse_args()
sample_ids, external_ids, tumor_variant_files, sample_types = \
	assert_inputs(args.sample_ids, args.external_ids, args.tumor_variant_files, args.sample_types)
print(sample_ids, tumor_variant_files)
aggregate_snvs(external_ids, tumor_variant_files)
