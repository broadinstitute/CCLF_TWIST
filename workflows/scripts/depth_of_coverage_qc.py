import pandas as pd
import argparse

############################################
# Compute whether a sample has enough coverage to pass TSCA
# Inputs:
#  See args
# Outputs:
#	depth_of_cov.txt: constains string pass or fail
# mcadosch 8/17
############################################

parser = argparse.ArgumentParser(description='Depth of coverage QC')
parser.add_argument('--gene_summary_file', type=str, required=True,
                    help='*.sample_gene_summary file')
parser.add_argument('--interval_summary_file', type=str, required=True,
                    help='*.sample_interval_summary file')
parser.add_argument('--min_depth', type=int, required=True, 
					help='minimum coverage to pass QC')

args = parser.parse_args()
gene_summary = pd.read_table(args.gene_summary_file)
interval_summary = pd.read_table(args.interval_summary_file)

# Pass the QC if either gene-level or interval-level mean coverage above threshold
pass_qc  = ( gene_summary['average_coverage'].mean() >= args.min_depth) \
				or ( interval_summary['average_coverage'].mean() >= args.min_depth )

out_file = open("depth_of_cov_result.txt", "w")
if pass_qc:
	out_file.write("pass")
else:
	out_file.write("fail")

mean_gene_coverage 				= gene_summary['average_coverage'].mean()
mean_interval_coverage 			= interval_summary['average_coverage'].mean()

# Write out results
mean_gene_coverage_file 		= open("mean_gene_cvg.txt", "w")
mean_interval_coverage_file 	= open("mean_interval_cvg.txt", "w")

mean_gene_coverage_file.write( "%s"%mean_gene_coverage )
mean_interval_coverage_file.write( "%s"%mean_interval_coverage )

out_file.close()
mean_gene_coverage_file.close()
mean_interval_coverage_file.close()