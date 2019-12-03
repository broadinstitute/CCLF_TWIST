import argparse
import sys
import os
import pandas as pd

############################################
# Aggregate data depth of coverage data for all samples
# Inputs:
#  See args
# Outputs:
#	Combined {gene/interval}-level coverage for all samples in {tsca_id}
#	{tsca_id}.depth_of_cov_by_{gene/interval}_cvg.txt
# mcadosch 8/17
############################################

### Parse command line inputs
parser = argparse.ArgumentParser(description='Aggregate depth of coverage statistics for many samples')
parser.add_argument('-s', '--statistic', type=str, required=True,
                    help='Type of statistic to look at: mean or total')
parser.add_argument('-i', '--interval', type=str, required=True,
                    help='Interval to look at: gene or target')
parser.add_argument('-f', '--files', type=str, required=True, nargs='+',
					help='Depth of coverage files (*.sample_interval_summary or *.sample_gene_summary)')
parser.add_argument('--tsca_id', type=str, required=True,
					help='TSCA ID')
parser.add_argument('--sample_ids', type=str, required=True, nargs='+', 
					help='sample ids')
args = parser.parse_args()


####################################################################################################
# 1. Aggregate data depth of coverage data for all samples
# Data is in files end in *.sample_interval_summary
# Aggregate based on statistic (mean/total coverage) and interval (target/gene)
####################################################################################################

# The statistic to aggregate depends on what was passed in the command line
if args.statistic == "mean":
	statistic = "average_coverage"
elif args.statistic == "total":
	statistic = "total_coverage"
# The statistic to aggregate depends on what was passed in the command line
if args.interval == "target":
	interval = "Target"
elif args.interval == "gene":
	interval = "Gene"


# Read in data first sample
df0 = pd.read_table(args.files[0], index_col=None, header=0, usecols=[interval, statistic]) \
        .rename(columns={statistic:args.sample_ids[0]})

dfs = [df0]

# Append data for remaining samples
for f, sid in zip(args.files[1:], args.sample_ids[1:]):
    df_to_aggregate = pd.read_table(f, index_col=None, header=0, usecols=[statistic]) \
                        .rename(columns={statistic:sid})
    dfs.append(df_to_aggregate)

df = pd.concat(dfs, axis=1)


filename = "%s.depth_of_cov_by_%s.%s_cvg.txt" %(args.tsca_id, args.interval, args.statistic)
df.to_csv(filename, sep="\t", index=None)