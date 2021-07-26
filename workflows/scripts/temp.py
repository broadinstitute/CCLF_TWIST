def plot_raw_cnv_calls(sample_ids, external_ids, sample_types, participant_ids, files):
    """Plot raw cnv calls
    Args:
        - files: files to use, either TN (tumor-normalized) and PTN (pre-tangent-normalized)
    """
    tsca_id = "test_batch"
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
    df = pd.read_csv(files[0],
                         index_col=None, header=0, comment='#',
                         usecols=['contig', 'start', 'stop', 'name', sample_ids[0]],
                         sep = '\t'
                         ).rename(columns={sample_ids[0]: external_ids[0]})

    # Append data for remaining samples by taking the intersection of intervals found in all samples
    for f, sid, eid, in zip(files[1:], sample_ids[1:], external_ids[1:]):
        df1 = pd.read_csv(f,comment="#", sep='\t').rename(columns={sid: eid})
        df = pd.merge(df, df1, on=["contig", "start", "stop", "name"], how="inner")

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

    # sort rows by chrm
    df.sort_values(by='chrm_int', ascending=False, inplace=True)
    df = df.drop('chrm_int', axis=1)

    ################################################
    # Save raw data to file
    ################################################
    print("Saving raw data to file...")
    fname = tsca_id + ".cnv_calls_unsegmented"
    df.to_csv(fname + ".txt", sep="\t", index=False)
    print(fname + ".txt")

    ################################################
    # Save raw data to file with sample ids
    ################################################
    df_sample_ids = df.rename(columns=dict(zip(external_ids, sample_ids)))
    df_sample_ids.to_csv(fname + ".sample_ids.txt", sep="\t", index=False)
    print(fname + ".sample_ids.txt")
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

    # need to select external ids in correct order here; above sorting of the df doesn't matter at all
    subdf = df[df.columns[4:]].T
    subdf['is_normal'] = [i == "Normal" for i in sample_types]
    subdf['participant_id'] = participant_ids
    subdf['external_id'] = external_ids

    # sorted columns
    external_ids_sorted = subdf.sort_values(by=['is_normal', 'participant_id', 'external_id'], ascending = [False, True, True]).loc[:,'external_id']
    participant_ids_sorted = subdf.sort_values(by=['is_normal', 'participant_id', 'external_id'], ascending = [False, True, True]).loc[:,'participant_id']

    # Create tick arrays for plot
    # Note that these ticks will be evenly spaced even if they represent intervals of differing size
    # TODO: space the ticks to match the size of the interval
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
    # with all normals on the left of the first figure, then the remainder sorted by participant_id and external_id
    # and split into multiple figures if there are more than samples_per_fig samples
    for fig_num in np.arange(num_figs):

        # make x-axis labels with format "participant_id: external_id"
        fig_external_ids = external_ids_sorted[fig_num * samples_per_fig: (fig_num + 1) * samples_per_fig]
        fig_participant_ids = participant_ids_sorted[fig_num * samples_per_fig: (fig_num + 1) * samples_per_fig]
        labels_for_xaxis = [str(m)+': '+str(n) for m,n in zip(fig_participant_ids, fig_external_ids)]

        axs[fig_num].pcolor(df[fig_external_ids].values, cmap=plt.cm.RdBu_r, vmin=-2, vmax=2)
        axs[fig_num].set_yticklabels(chromosomes, minor=False)
        axs[fig_num].set_yticks(tick_positions)
        axs[fig_num].set_xticklabels(labels_for_xaxis, rotation=90, ha="left")
        axs[fig_num].set_xticks(range(len(df[fig_external_ids].columns.tolist())))

    fig.subplots_adjust(hspace = 0.15)
    fig.savefig(fname + ".pdf")
    print(fname+".pdf")
#     os.system('ls -alh')


    return(fig)


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
    participant_ids = [participant_ids[i] for i in indices_samples_to_keep]
    sample_types = [sample_types[i] for i in indices_samples_to_keep]

    return sample_ids, external_ids, sample_types, participant_ids, files
