import matplotlib
matplotlib.use('Agg')
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
import matplotlib.pyplot as plt

############################################
# Estimate gender from depth of coverage data
# Inputs:
#   See args
# Outputs:
#   ./gender.png: distribution of coverage across chromosomes
#   ./gender.txt: gender M/F/NA
# mcadosch 10/17
############################################


#######################################################################################################
# 0. Parse command line inputs
#######################################################################################################
parser = argparse.ArgumentParser(description='Estimate gender from depth of coverage data')
parser.add_argument('--sample_id', type=str, required=True,
                    help='sample_id')
parser.add_argument('--external_id', type=str, required=True,
                    help='external_id')
parser.add_argument('--depth_of_cov_by_interval', type=str, required=True,
                    help='Depth of coverage by interval (sample_interval_summary)')

args = parser.parse_args()
sample_id                       = args.sample_id
external_id                     = args.external_id
depth_of_cov_by_interval_path   = args.depth_of_cov_by_interval

#######################################################################################################
# 1. Analysis
#######################################################################################################
def estimate_gender(data, sample_id, external_id):
    """Estimate gender
    Args:
        - data: depth of coverage by interval
    """
    
    # Non-gender chromosomes
    NG = data.loc[~(data['Target'].str.contains('X')) & ~(data['Target'].str.contains('Y')), :]
    NG.loc[:, 'chr'] = NG['Target'].apply(lambda x: x.split(':')[0])
    # X and Y chromosomes
    X = data.loc[data['Target'].str.contains('X'), :]
    Y = data.loc[data['Target'].str.contains('Y'), :]
    
    # Mean coverage for non-gender chromosomes
    NG_mean_coverage = NG.groupby('chr').mean()
    
    # Compute p-value
    mu = NG_mean_coverage['average_coverage'].mean()
    sigma = NG_mean_coverage['average_coverage'].std()
    X_mean_coverage = X['average_coverage'].mean()
    Y_mean_coverage = Y['average_coverage'].mean()
    p_val_x = stats.norm(mu, sigma).cdf(X_mean_coverage)
    p_val_y = stats.norm(mu, sigma).cdf(Y_mean_coverage)
    
    # Plot distribution of chromosome coverages
    fig, ax = plt.subplots(figsize=(16, 8))
    sns.distplot(NG_mean_coverage, ax=ax)
    ax.axvline(X_mean_coverage, c='r', ls='--', label='mean coverage X chromosome (p-val: %s)'%p_val_x)
    ax.axvline(Y_mean_coverage, c='g', ls='--', label='mean coverage Y chromosome (p-val: %s)'%p_val_y)
    ax.set_title('Chromosome mean coverage distribution for sample %s'%external_id)
    ax.legend(loc='upper right')

    fname = "%s.chromosome_cov_distribution"% (sample_id)
    fig.savefig("%s.png"%fname)
    
    # File with gender estiamte
    out_file = open("gender_estimate.txt", "w")
    alpha = 0.05
    
    # If p_val_x < 0.05, reject hypothesis that chromosome X comes from same distribution of coverages as non-gender chromosomes 
    # Since males have only one copy of X chromosome (vs. 2 copies of non-gender chromosomes), sample is most likely male.
    if p_val_x <= alpha:
        out_file.write("M")
    else:
        out_file.write("F")

    out_file.close()
    return

####################################################################################################
# MAIN
####################################################################################################
if __name__ == '__main__':
    # Read data
    depth_of_cov = pd.read_table(depth_of_cov_by_interval_path, usecols=['Target', 'average_coverage'])
    estimate_gender(depth_of_cov, sample_id, external_id)
