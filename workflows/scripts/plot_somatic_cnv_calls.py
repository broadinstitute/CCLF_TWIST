import matplotlib
matplotlib.use('Agg')
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import math

############################################
# Plot somatic CNV heatmap
# Edited from the mcadosch/TSCA docker image from 20191008
# Inputs:
#   See args
# Outputs:
#   ./tscaXX.cnv_calls.png: CNV calls for all intervals
#   ./tscaXX.cnv_calls.txt: Raw data for heatmap
# mcadosch 8/17
############################################


#######################################################################################################
# 0. Parse command line inputs
#######################################################################################################
parser = argparse.ArgumentParser(description='Plot Results from Somatic CNV Calls')
parser.add_argument('--tsca_id', type=str, required=True,
                    help='TSCA ID')
parser.add_argument('--tn_files', type=str, required=True, nargs='+',
                    help='List of *tn.tsv files')
parser.add_argument('--sample_ids', type=str, required=True, nargs='+',
                    help='List of sample ids')
parser.add_argument('--external_ids', type=str, required=True, nargs='+',
                    help='List of external validation ids')
parser.add_argument('--participant_ids', type=str, required=True, nargs='+',
                    help='List of patient ids')
parser.add_argument('--depth_of_cov_qcs', type=str, required=True, nargs='+',
                    help='List of results from depth of coverage QC (pass/fail)')
parser.add_argument('--sample_types', type=str, required=True, nargs='+',
                    help='List of sample types (Tumor vs. Normal)')
args = parser.parse_args()

####################################################################################################
# Remove normal samples
####################################################################################################

# tumor_indices   = [i for i, x in enumerate(args.sample_types) if x=='Tumor']
# sample_ids      = [args.sample_ids[i] for i in tumor_indices]
# external_ids    = [args.external_ids[i] for i in tumor_indices]
# tn_files        = [args.tn_files[i] for i in tumor_indices]

tsca_id = args.tsca_id
tn_files = args.tn_files
sample_ids = args.sample_ids
external_ids = args.external_ids
participant_ids = args.participant_ids
depth_of_cov_qcs = args.depth_of_cov_qcs
sample_types = args.sample_types


def argsort(seq):
    # http://stackoverflow.com/questions/3071415/efficient-method-to-calculate-the-rank-vector-of-a-list-in-python
    return sorted(range(len(seq)), key=seq.__getitem__)

####################################################################################################
# 1. Plot the CNV calls for all intervals
####################################################################################################


def assert_inputs(sample_ids, external_ids, sample_types, participant_ids, files, depth_of_cov_qcs):
    """Ensure inputs are correct, and that the external_ids, sample_ids, sample_types,
    participant_ids, files arrays correspond to each other
    Returns:
        - True / False if inputs are correct
    """
    print("Asserting inputs...")
    print("Number of external_ids: %s \n Number of sample_ids: %s \n Number of sample_types: %s \n \
        Number of participant_ids: %s \n Number of segmented files: %s \
                     \n Number of depth of coverage QCs: %s" % (len(external_ids), len(sample_ids),
                                                                len(sample_types), len(participant_ids), len(files), len(depth_of_cov_qcs)))
    # check: will the normals have a patient ID? That is, do we expect any samples to NOT have participant_ids?
    if len(external_ids) == len(sample_ids) == len(sample_types) == len(participant_ids) == len(files) == len(depth_of_cov_qcs):
        print("Same number of sample_ids, external_ids, sample_types, participant_ids, files, \
            and depth_of_cov_qc results passed. Inputs are correct...")
        return True
    else:
        print("Inputs have different lengths. One of the inputs is wrong...")
        return False


def plot_raw_cnv_calls(sample_ids, external_ids, sample_types, participant_ids, files):
    """Plot raw cnv calls
    Args:
        - files: files to use, either TN (tumor-normalized) and PTN (pre-tangent-normalized)
    """
    ################################################
    # Create data df
    # Replace sample_id with external_id, as external_ids are more informative
    # check: I think I also need to include the participant_ids as a column for eventual sorting
    ################################################
    print("Plotting unsegmented CNV calls...")
    print("There are %d sample_ids, %d external_ids, %d  sample_types, %d participant_ids, %d files" %
          (len(sample_ids), len(external_ids), len(sample_types), len(participant_ids), len(files)))

    # First DF (necessary to establish)
    # check: need to add participant_ids / pid
    dfs = [pd.read_table(files[0],
                         index_col=None, header=0, comment='#',
                         usecols=['contig', 'start', 'stop', 'name', sample_ids[0]]
                         ).rename(columns={sample_ids[0]: external_ids[0]})]

    # Append data for remaining samples
    for f, sid, eid, in zip(files[1:], sample_ids[1:], external_ids[1:]):
        print(sid, eid)
        df_to_add = pd.read_table(f, index_col=None, header=0, comment='#',
                                  usecols=[sid]).rename(columns={sid: eid})
        dfs.append(df_to_add)

    # NOTE: In order to concatente samples into a single table,
    # we need to shave off the extra rows in the samples with extra rows
    # check: why would there ever be any extra rows?
    n_rows = min([df.shape[0] for df in dfs])

    # Merge dfs
    df = pd.concat([df.iloc[0:n_rows] for df in dfs], axis=1)

    ################################################
    # Creating chromosome labels for plot
    ################################################
    # Chromosome as strings 1, 2, ..., X, Y
    df = df.rename(columns={'contig': 'chrm_str'})
    # Chromosome as ints 1, 2, ..., 23, 24 (Used for sorting)
    df['chrm_int'] = df['chrm_str']
    # Set X, Y to 23, 24
    df.loc[df['chrm_int'] == 'X', 'chrm_int'] = 23
    df.loc[df['chrm_int'] == 'Y', 'chrm_int'] = 24
    df['chrm_int'] = df['chrm_int'].astype(int)
    # Sort values
    df.sort_values(by='chrm_int', ascending=False, inplace=True)

    ################################################
    # Sort datatable by patient_id, external_id
    # want nested sorting: sorted by chrm values, then patient_id, then external_id.
    # also want all the normals on the left
    # why sort? aids interpretation of the CNV plots
    ################################################
    # make boolean column identifying the normal samples; want to eventually have these on the far left of the CNV plot
    df['is_normal'] = [i == "Normal" for i in sample_types]  # check:/fix
    df['participant_id'] = participant_ids
    df['external_id'] = external_ids
    df[df.columns[4:]] = df[df.columns[4:]].T.sort_values(by=['is_normal', 'participant_id', 'external_id']).T
    df.drop(['is_normal', 'participant_id'], axis=1, inplace=True)
    # even if have patient specific normal, I think I need to include it on the left hand side
    # thus put boolean column before the patient_id

    # perform nested sort in the dataframe so that the CNV plot we create will be sorted by normals
    # on left, then sorted by patient id (and sorted by ext_id within pt_id)
    df.sort_values(by=['chrm_int'], ascending=False, inplace=True)
    # get rid of the column we added for sorting

    ################################################
    # Save raw data to file
    ################################################
    print("Saving raw data to file...")
    fname = "%s.cnv_calls_unsegmented" % (args.tsca_id)
    df.to_csv("./%s.txt" % fname, sep="\t", index=False)

    ################################################
    # Save raw data to file with sample ids
    ################################################
    df_sample_ids = df.rename(columns=dict(zip(external_ids, sample_ids)))
    df_sample_ids.to_csv("./%s.sample_ids.txt" % fname, sep="\t", index=False)

    # ################################################
    # ## Save raw data to file with external ids
    # ################################################
    # df_external_ids = df.rename(columns=dict(zip(sample_ids, external_ids)))
    # df_external_ids.to_csv("%s.cnv_calls_unsegmented_with_external_id.txt"%args.tsca_id, sep="\t", index=False)

    ################################################
    # Plot and save figure
    ################################################
    # list of chromosomes
    chromosomes = df['chrm_str'].unique().tolist()
    n_intervals_per_chromosome = [df['chrm_str'].value_counts()[i] for i in chromosomes]

    # Create tick arrays for plot
    tick_positions = np.cumsum(n_intervals_per_chromosome)

    # Divide samples into multiple figures if too many samples
    samples_per_fig = 80
    num_figs = int(math.ceil(float(len(external_ids)) / samples_per_fig))
    print("num_figs: ", num_figs)
    fig_height = int(30 * num_figs)
    fig, axs = plt.subplots(num_figs, 1, figsize=(25, fig_height))
    print(axs)
    if num_figs > 1:
        axs = axs.ravel()
    else:
        axs = [axs]
    print("num_figs: ", num_figs)
    print("axs: ", axs)

    # Draw figure for each sample set
    # check: want to have all normals on the left, then sorted by participant_id and external_id
    for fig_num in np.arange(num_figs):
        fig_external_ids = external_ids[fig_num * samples_per_fig: (fig_num + 1) * samples_per_fig]
        axs[fig_num].pcolor(df[fig_external_ids].values, cmap=plt.cm.RdBu_r, vmin=-2, vmax=2)
        axs[fig_num].set_yticklabels(chromosomes, minor=False)
        axs[fig_num].set_yticks(tick_positions)
        axs[fig_num].set_xticklabels(df[fig_external_ids].columns.tolist(), rotation=90, ha="left")
        axs[fig_num].set_xticks(range(len(df[fig_external_ids].columns.tolist())))

    fig.subplots_adjust(bottom=0.55)
    fig.savefig("%s.png" % fname)

    return


def remove_samples_low_coverage(sample_ids, external_ids, sample_types, participant_ids, files, depth_of_cov_qcs):
    """Remove samples with low coverage
    """
    print("Removing samples with low coverage...")
    indices_samples_to_exclude = [idx for (idx, qc) in enumerate(depth_of_cov_qcs) if qc == 'fail']
    samples_excluded = [sample_ids[i] for i in indices_samples_to_exclude]
    print("Excluding samples: %s" % (samples_excluded))

    indices_samples_to_keep = [idx for idx, (sid, eid, stype, pid, f, qc) in enumerate(zip(sample_ids,
                                                                                    external_ids, sample_types, participant_ids, files, depth_of_cov_qcs)) if qc == 'pass']
    sample_ids = [sample_ids[i] for i in indices_samples_to_keep]
    external_ids = [external_ids[i] for i in indices_samples_to_keep]
    files = [files[i] for i in indices_samples_to_keep]
    participant_ids = [participant_ids[i] for i in indices_samples_to_keep]  # check: do these match up?
    sample_types = [sample_types[i] for i in indices_samples_to_keep]  # check: do these match up?

    return sample_ids, external_ids, sample_types, participant_ids, files


####################################################################################################
# MAIN
####################################################################################################
if __name__ == '__main__':
    print("Running plot_somatic_cnv_calls.py...\n")
    if not assert_inputs(sample_ids, external_ids, sample_types, participant_ids, tn_files, depth_of_cov_qcs):
        import sys
        sys.exit()
    sample_ids, external_ids, sample_types, participant_ids, tn_files = remove_samples_low_coverage(sample_ids,
                                                                                                    external_ids,
                                                                                                    sample_types,
                                                                                                    participant_ids,
                                                                                                    tn_files,
                                                                                                    depth_of_cov_qcs)
    plot_raw_cnv_calls(sample_ids, external_ids, sample_types, participant_ids, tn_files)
