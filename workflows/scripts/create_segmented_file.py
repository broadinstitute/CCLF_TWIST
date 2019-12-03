import matplotlib
matplotlib.use('Agg')
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


############################################
# Save segmented CNV data to file
# Inputs:
#   See args
# Outputs:
#   ./sample_id.complete_seg.txt: Raw data for heatmap (segmented with all intervals)
# mcadosch 10/17
############################################

#######################################################################################################
# 0. Parse command line inputs
#######################################################################################################
parser = argparse.ArgumentParser(description='Plot Results from Somatic CNV Calls')
parser.add_argument('--tn_file', type=str, required=True, nargs=1,
                    help='*.tn.tsv file')
parser.add_argument('--seg_file', type=str, required=True, nargs=1,
                    help='*.called file')
parser.add_argument('--sample_id', type=str, required=True, nargs=1, 
                    help='sample id')
parser.add_argument('--external_id', type=str, required=True, nargs=1, 
                    help='external validation id')
parser.add_argument('--out_dir', type=str, required=False, default=".", nargs=1,
                    help='Output directory')

args = parser.parse_args()

tn_file            	= args.tn_file[0]
seg_file           	= args.seg_file[0]
sample_id          	= args.sample_id[0]
external_id        	= args.external_id[0]
out_dir				= args.out_dir[0]

####################################################################################################
# 0. Add Segment Mean to Non-Segmented File
####################################################################################################
def copy_segmented_mean(sample_id, tn_file, seg_file):
    """For the non-segmented file (*.tn.tsv), add a 'segment_mean' column with the segment_mean from the segmented file (*.called) for every target
    """
    non_segmented = pd.read_table(tn_file, comment="#")
    segmented = pd.read_table(seg_file)
    chromosomes = non_segmented.contig.unique().tolist()

    # Rename sample_id column (containing tangent normalized values) to tn_values
    non_segmented.rename(columns={sample_id: 'tn_values'})

    # Iterate over chromosomes
    for chrm in chromosomes:
        # Targets in chromosome
        chrm_targets = non_segmented[non_segmented.contig==chrm]
        for idx, target in chrm_targets.iterrows():
            # Find segment_mean for target in non-segmented df, and set it to the non-segmented df target value
            seg_mean = segmented.loc[(segmented.Chromosome==chrm) & (segmented.Start <= target['start']) \
                     & (segmented.End >= target['stop']), 'Segment_Mean']
            non_segmented.loc[idx, sample_id] = seg_mean.item()

    return non_segmented

def save_segmented_cnv_calls(sample_id, external_id, tn_file, seg_file):
    """Save CNV data
    """
    df = copy_segmented_mean(sample_id, tn_file, seg_file)
    df.drop('name', axis=1, inplace=True)
    fname = "%s/%s_complete_seg.txt" %(out_dir, sample_id)
    df.to_csv(fname, sep="\t", index=None)
    print("Saving CNV Segmented data for %s to %s" %(sample_id, out_dir))
    return

####################################################################################################
# MAIN
####################################################################################################
if __name__ == '__main__':
    merged_df = save_segmented_cnv_calls(sample_id, external_id, tn_file, seg_file)