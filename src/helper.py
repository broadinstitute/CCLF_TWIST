 # Helper functions for running CCLF QC pipeline
import pandas as pd
import numpy as np
import dalmatian as dm
from gsheets import Sheets
import sys
pathtoJK_parent = "../"
sys.path.append(pathtoJK_parent)
from JKBio import terra


sheets = Sheets.from_files('~/.client_secret.json', '~/.storage.json')


def getCohortAbbreviations(df, cohorts2id_url="https://docs.google.com/spreadsheets/d/1R97pgzoX0YClGDr5nmQYQwimnKXxDBGnGzg7YPlhZJU"):
    # mapping abbreviations to full names/descriptions
    cohorts = sheets.get(cohorts2id_url).sheets[0].to_frame()
    # match collection data and error out
    cohortlist = []
    no_cohort_match = set()
    for k, val in df['Collection'].iteritems():
        res = cohorts[cohorts['Name'] == val]
        if len(res) == 0:
            print(res)
            no_cohort_match.add(str(val))
            cohortlist.append('nan')
        else:
            cohortlist.append(res['ID'].values[0])
    print("we do not have a corresponding cohort abbreviation for these cohorts: [{0}]".format(", ".join(str(i) for i in no_cohort_match)))

    df['cohorts'] = cohortlist
    return df


def regenerate_variables(wm, samplesetnames, cohorts2id_url):
    # Run function when need to re-generate key variables: sample_info, all_pairsets, the cohorts_per_batch dictionary, cohort_pairsets
    # Inputs:
    # - wm: Dalmation workspace manager for the data processing workspace

    # wm = dm.WorkspaceManager(workspace)
    # will use "cohort_pairsets" to create cohort-specific SNV tsv and CN heat map
    cohorts = sheets.get(cohorts2id_url).sheets[0].to_frame()
    sample_info = wm.get_samples() # TODO: check whether to use this or not
    cohorts_per_batch = {} # will be dict of cohorts in each batch
    all_changed_cohorts = set()
    for i in range(len(samplesetnames)):
        # get appropriate subset of the samples for each batch
        batch_sample_info = sample_info[sample_info['batch'] == samplesetnames[i][:-4]]
        cohorts_in_batch = set()
        cohorts_with_pairs = [] # check: currently not used.
        # for each batch, make pairsets by cohort
        for val in cohorts['ID'].values:
            cohortsamples = batch_sample_info[batch_sample_info["cohorts"] == val].index.tolist()
            tumorsamplesincohort = batch_sample_info[batch_sample_info["cohorts"] == val][batch_sample_info['sample_type']=="Tumor"].index.tolist()
            if len(cohortsamples)>0:
                cohorts_in_batch.update([val])
        batch_name = samplesetnames[i]
        cohorts_per_batch[batch_name] = cohorts_in_batch
        all_changed_cohorts.update(cohorts_in_batch) # add all the new cohorts in this batch to the full list

    # list of the cohort pairsets affected by the new samples we added
    all_pairsets = wm.get_pair_sets().index.tolist()
    cohort_pairsets = set(all_changed_cohorts) - (set(all_changed_cohorts) - set(all_pairsets))

    return sample_info, all_pairsets, cohorts_per_batch, cohort_pairsets, all_changed_cohorts


def check_cohort_abbrev_availability(df, cohorts2id_url):
    _, no_cohort_match = get_list_of_cohort_abbrevs_in_df(df, cohorts2id_url)
    return no_cohort_match


def get_list_of_cohort_abbrevs_in_df(df, cohorts2id_url):
    # the df must contain the columns: `Collection` (this is found in the External ID sheet)

    # cohorts and cohort abbreviations
    cohorts = sheets.get(cohorts2id_url).sheets[0].to_frame()

    # match collection data and error out
    cohortlist = []
    no_cohort_match = set()
    for k, val in df['Collection'].iteritems():
        res = cohorts[cohorts['Name'] == val]
        if len(res) == 0:
            no_cohort_match.add(str(val))
            cohortlist.append('nan')
        else:
            cohortlist.append(res['ID'].values[0])

    assert len(no_cohort_match) == 0,\
        "We do not have a corresponding cohort abbreviation for these cohorts: [{0}]".format(", ".join(str(i) for i in no_cohort_match)) + \
        "\nStop running the pipeline and have the CCLF team fix the missing cohort information."
    return cohortlist, no_cohort_match

def create_preliminary_sample_and_metadata_tables(wto, wfrom, samplesetnames, external_sheets_url_list, cohorts2id_url, forcekeep=[]):
    # Inputs:
    # - wto: Dalmation workspace manager for the data processing workspace
    # - wfrom: Dalmation workspace manager for the data delivery workspace


    # we look at all the samples we already have in the CCLF QC processing workspace
    refsamples = wto.get_samples()
    refids = refsamples.index

    # get the External sheet data from google sheets
    gsheets = [sheets.get(url).sheets[0].to_frame() for url in external_sheets_url_list]

    # add a column with batch information (e.g. CCLF_TWIST1 vs CCLF_TWIST2)
    metadata = pd.concat(gsheets,sort=False, keys = samplesetnames)
    metadata = metadata.reset_index().rename(columns = {'level_0':'batch', "External ID":'external_id_validation'}).drop(['level_1'], axis = 'columns')


    print('Number of rows from the external ID gsheets: {}'.format(len(metadata)))

    # we use this gsheet package to get all the sheets into one dataframe
    # cohorts
    cohorts = sheets.get(cohorts2id_url).sheets[0].to_frame()

    # we look at all the samples we already have in Terra
    # we do some corrections just in case
    samples1 = wfrom.get_samples().replace(np.nan, '', regex=True)

    # creating sample_id (like in processing wm) for metadata and samples1
    newmetadata = metadata.dropna(0, subset=['Collaborator Sample ID','Sample Type','Exported DNA SM-ID'])
    # print("Dropped indices: "+str(set(metadata.index.tolist())-set(newmetadata.index.tolist())))
    print("Dropped samples: "+str(set(metadata['Exported DNA SM-ID'].tolist())-set(newmetadata['Exported DNA SM-ID'].tolist())))
    print("Note that we expect to drop one per batch, the technical sequencing control.")

    print('New length after dropping any samples with no value for the Collaborator Sample ID, Sample Type, or Exported DNA SM-ID: '+str(len(newmetadata)))
    metadata=newmetadata

    ttype = [i for i in metadata["Sample Type"]]
    metadata['sample_id'] = [str(val['Collaborator Sample ID'][:-1]) + '-' + str(val['Sample Type']) + '-' + str(val['Exported DNA SM-ID']) for i, val in metadata.iterrows()]

    samples1.index = [i.split('_')[2] for i, val in samples1.iterrows()]

    samples1['sample_id'] = [str(val["individual_alias"]) + '-' + str(val['sample_type']) + '-' + i for i, val in samples1.iterrows()]

    print("Number of samples already in Terra: ", samples1.index.isin(refids).tolist().count(False))

    metadata.index = metadata['Exported DNA SM-ID']
    # print("Number of samples in the input External ID metadata: ",len(metadata.index.tolist()))

    # filtering on what already exists in the processing wm (refids)
    newsamples = samples1[(~samples1.index.isin(refids)) | samples1.index.isin(forcekeep)]
    post_filter_len = len(metadata.index.tolist())
    print('Number of samples remaining after filtering out samples already in the processing wm: '+ str(post_filter_len))


    # keep the new samples that we have metadata for
    tokeep = set(metadata.index) & set(newsamples.index)
    if len(tokeep) != post_filter_len:
        print("\nThe data delivery workspace (which is used to create `newsamples`) doesn't contain all the samples specified in the External Sheet(s) (which are used to create `metadata`). \nThe",len(set(metadata.index) - set(newsamples.index)),"samples we couldn't find are:", set(metadata.index) - set(newsamples.index), "\nThis can occur if a line is blacklisted due to too few reads, for example.")



    # useful to merge the two df, sm-id is one of the only unique id here
    num_newsamples = len(newsamples[~newsamples.index.isin(tokeep)])
    if num_newsamples > 0:
        print('We could not add these ' + str(num_newsamples) + ' samples from the data delivery wm as we don\'t have metadata for them: ' + '\n'
              + str(newsamples[~newsamples.index.isin(tokeep)].index))
    newsamples = newsamples[newsamples.index.isin(tokeep)]
    newmetadata = metadata[metadata.index.isin(tokeep)].sort_index().drop_duplicates("Exported DNA SM-ID")
    print("Final dimensions of the newsamples table: ", newsamples.shape)
    print("Final dimensions of the newmetadata table: ", newmetadata.shape)

    return newsamples, newmetadata


def check_required_metadata_columns(df, required_metadata_cols, cohorts2id_url, drop=False):
    # check collection (cohort) data and raise warning if missing
    check_cohort_abbrev_availability(df, cohorts2id_url)

    # If samples in the df are lacking any of the required metadata columns, flag these samples.
    samplesMissingMetadata = df.iloc[[j for j,i in enumerate(df[required_metadata_cols].isna().values.sum(1)) if i !=0]].index.tolist()
    nMissing = len(samplesMissingMetadata)

    # Return if no samples are missing metadata
    if nMissing == 0:
        print("No samples are missing any required metadata. Continue running the pipeline.")
        return df

    print("Since they do not have full data, we will be dropping " + str(nMissing) + " samples: \n" +  str(samplesMissingMetadata))

    # examine which columns in the External Sheet had missing information
    print('\nNumber of NAs for each required column:')
    print(df[required_metadata_cols].isna().sum())

    assert drop, '\nStop running the pipeline. Ask the CCLF team to fill out the missing values in the External Sheet.'

    if drop:
        # only keep samples that have all the appropriate information
        df = df.iloc[[j for j,i in enumerate(df[required_metadata_cols].isna().values.sum(1)) if i == 0]]
        print('We have now dropped the samples for which we didn\'t have all the required metadata columns.')

    return df


def create_sample_df_for_terra(wm, df, cohorts2id_url, picard_aggregation_type_validation='PCR', source='CCLF'):
    # Create sample df to upload to Terra
    # Inputs: df containing filtered set of samples plus metadata, cohorts and abbreviations gsheet url
    # - wm: Dalmation workspace manager for the data processing workspace
    # Output: sample_info table

    # cohorts and cohort abbreviations
    cohorts = sheets.get(cohorts2id_url).sheets[0].to_frame()
    cohortlist, _ = get_list_of_cohort_abbrevs_in_df(df, cohorts2id_url)

    # we look at all the samples we already have in the data processing workspace
    refsamples = wm.get_samples()

    print('creating new sample information df')
    # from this filtered set of samples (df) we create a dataframe which will get uploaded to terra
    sample_info = df[['crai_or_bai_path', 'cram_or_bam_path']].copy()
    sample_info['batch'] = df['batch'].astype(str)
    sample_info['individual_id'] = df['Collaborator Participant ID'].astype(str)
    sample_info['reference_id'] = df['Exported DNA SM-ID'].astype(str)
    sample_info['patient_id'] = df['Participant ID'].astype(str)
    sample_info['participant'] = df['Collaborator Participant ID'].astype(str)
    sample_info['aggregation_product_name_validation'] = df['bait_set'].astype(str)
    # here we add this number as the reference id might be present many times already for different samples
    # in the processing workspace

    # start building external_id_validation column
    sample_info['external_id_validation'] = 'nan'
    # for each SM-ID:
    for i in range(len(sample_info['reference_id'])):
        # get the external id for the sample
        ext_id_for_sample = df[df.index == sample_info['reference_id'][i]]['external_id_validation'].values[0]
        # replace any "/" that exist with "_"; otherwise get errors because looks like new directory when try to build file paths
        ext_id_for_sample = ext_id_for_sample.replace('/', '_') #[ext_id_for_sample.replace('/', '_') for ext_id in ext_id_for_sample]

        # tack on a number to distinguish external IDs that we have run more than once
        # using str.contains because we want to ignore the tacked on numbers we've added to the ext_id (e.g. _1, _2)
        # num of samples with this ext_id already in the workspace
        num_in_workspace = refsamples[refsamples.external_id_validation.str.contains(ext_id_for_sample)].shape[0]
        # num of samples with this ext_id that we've already seen in the data we're adding
        try:
            num_already_seen_here = sample_info[sample_info.external_id_validation.str.contains(ext_id_for_sample)].shape[0]
            num_to_add = num_in_workspace + num_already_seen_here + 1
        except:
            num_to_add = num_in_workspace + 1
        sample_info['external_id_validation'][i] = ext_id_for_sample + '_' + str(num_to_add)

    sample_info['bsp_sample_id_validation'] = df.index.str.strip()
    sample_info['stock_sample_id_validation'] = df['Stock DNA SM-ID'].str.strip()
    sample_info['sample_type'] = df['Sample Type'].str.strip()
    # TODO: I don't know what the "picard_aggregation_type_validation" is used for. Right now, the value is completely meaningless.
    sample_info['picard_aggregation_type_validation'] = [picard_aggregation_type_validation] * sample_info.shape[0]
    sample_info['tumor_subtype'] = df['Tumor Type'].str.strip()
    sample_info['squid_sample_id_validation'] = sample_info['external_id_validation']
    sample_info['source_subtype_validation'] = df['Original Material Type'].str.strip()
    sample_info['processed_subtype_validation'] = df['Material Type'].str.strip()
    sample_info['primary_disease'] = df['Primary Disease'].str.strip()
    sample_info['media'] = df['Media on Tube'].str.strip()
    sample_info['Collection'] = df['Collection'].str.strip()
    sample_info['cohorts'] = cohortlist
    sample_info['tissue_site'] = df['Tissue Site'].astype(str)
    sample_info['source'] = [source] * sample_info.shape[0]
    sample_info['sample_id'] = df.index.astype(str)

    sample_info = sample_info.set_index('sample_id')

    return sample_info


######## Functions for creating the pairs, pairsets, samplesets

def create_pairs_table(sample_info, wm):
    # Function to create the pairs table to upload to Terra.
    # Inputs:
    # - sample_info: the sample info table for the new samples that will be uploaded to Terra
    # - wm: Dalmation workspace manager for the data processing workspace
    refsamples = wm.get_samples()

    normals = [r["participant"] for i, r in sample_info.iterrows() if r['sample_type'] == "Normal"]
    normalsid = [i for i, r in sample_info.iterrows() if r['sample_type'] == "Normal"]
    tumors = [r["participant"] for i, r in sample_info.iterrows() if r['sample_type'] == "Tumor"]
    tumorsid = [i for i, r in sample_info.iterrows() if r['sample_type'] == "Tumor"]
    prevtumors = [val["participant"] for k, val in refsamples.iterrows() if val.sample_type == "Tumor"]
    prevnormals = [val["participant"] for k, val in refsamples.iterrows() if val.sample_type == "Normal"]

    print("creating new pairs...")
    # do we have new tumors/normals for our previous ones
    newpairs = {'pair_id': [], 'case_sample': [], 'control_sample': [], 'participant': [], 'match_type':[]}

    toreprocess_normals = set(tumors) & set(prevnormals)
    for val in toreprocess_normals:
        if val != 'nan':
            for tumor_id in sample_info[sample_info['participant'] == val][sample_info[
                    'sample_type'] == 'Tumor'].index.tolist():
                normal_id = refsamples[refsamples['participant'] == val][refsamples[
                  'sample_type'] == 'Normal'].index.tolist()[0]
                newpairs['pair_id'].append(tumor_id + '_' + normal_id)
                newpairs['case_sample'].append(tumor_id)
                newpairs['control_sample'].append(normal_id)
                newpairs['participant'].append(val)
                newpairs['match_type'].append("Tumor_Normal")

    paired = set(tumors) & set(normals)
    for val in set(tumors) - toreprocess_normals:
        if val != 'nan':
            for tumor_id in sample_info[sample_info['participant'] == val][sample_info[
                    'sample_type'] == 'Tumor'].index.tolist():
                normal_id = sample_info[(sample_info['participant'] == val) & (sample_info[
                  'sample_type'] == 'Normal')].index.tolist()[0] if val in paired else 'NA'
                newpairs['pair_id'].append(tumor_id + "_" + normal_id)
                newpairs['case_sample'].append(tumor_id)
                newpairs['control_sample'].append(normal_id)
                newpairs['participant'].append(val)
                newpairs['match_type'].append("Tumor_Normal" if val in paired else 'Tumor_NA')

    newpairs = pd.DataFrame(newpairs).set_index('pair_id')
    print("Done")
    return newpairs


def create_dict_of_pairs_per_sampleset(pairs, sample_info, samplesetnames, save=True, filename='dict_pairs_per_batch.npy'):
    # Pair => case_sample => look into samples => retrieve the batch, assign to the key
    dict_pairs_per_batch = {}
    for samplesetname in samplesetnames:
        dict_pairs_per_batch[samplesetname] = []

    for pair_id, pair in pairs.iterrows():
        case_sample = pair['case_sample']
        # Retrieve the batch which sample belongs to
        pair_s_batch = sample_info.loc[case_sample]['batch']
        dict_pairs_per_batch[pair_s_batch].append(pair_id)

    if save:
        print('Saving dictionary of pairs per sampleset to %s' % filename)
        np.save(filename, dict_pairs_per_batch)

    return dict_pairs_per_batch


def create_dict_of_samples_per_sampleset(sample_info, samplesetnames, save=False, filename='dict_samples_per_batch.npy'):
    # Samples => retrieve the batch, assign to the key
    dict_samples_per_batch = {}
    for samplesetname in samplesetnames:
        dict_samples_per_batch[samplesetname] = []

    # Retrieve the batch which sample belongs to
    for sample_id, sample in sample_info.iterrows():
        batch = sample_info.loc[sample_id]['batch']
        dict_samples_per_batch[batch].append(sample_id)

    if save:
        print('Saving dictionary of samples per sampleset to %s' % filename)
        np.save(filename, dict_samples_per_batch)

    return dict_samples_per_batch


# TODO: WIP
# def create_dict_of_cohorts_per_sampleset(sample_info, samplesetnames, save=False, filename='dict_cohorts_per_batch.npy'):
#     # Samples => retrieve the batch, assign to the key
#     dict_samples_per_batch = {}
#     for samplesetname in samplesetnames:
#         dict_samples_per_batch[samplesetname] = []
#
#     # Retrieve the batch which sample belongs to
#     for sample_id, sample in sample_info.iterrows():
#         batch = sample_info.loc[sample_id]['batch']
#         dict_samples_per_batch[batch].append(sample_id)
#
#     if save:
#         print('Saving dictionary of samples per sampleset to %s' % filename)
#         np.save(filename, dict_samples_per_batch)
#
#     return dict_samples_per_batch

def create_aggregate_samplesets(normalsid, proc_workspace, wto):
    print("Creating and uploading sample sets for all samples in workspace, and all normals in workspace...")
    all_samples = wto.get_samples().index.tolist()
    all_samples.remove('NA')

    # remove duplicates
    normalsid = list(set(normalsid))

    # create sample set for all normals in workspace (will use all combined normals in PoN for mutation calling)
    terra.addToSampleSet(proc_workspace, samplesetid="All_normals_TWIST", samples=normalsid)
    terra.addToSampleSet(proc_workspace, samplesetid="All_samples_TWIST", samples=all_samples)
    return

def create_sample_sets_per_batch(sample_info, samplesetnames, samplesetnames_all, samplesetnames_tumors, samplesetnames_normals, proc_workspace):
    print("Creating and uploading sample sets for each batch...")
    # want to create a sample set for each batch
    for i, current_batch in enumerate(samplesetnames):
        # get appropriate subset of the samples
        batch_sample_info = sample_info[sample_info['batch'] == current_batch]
        # define batch-specific tumors and normals
        batch_normals = [r["participant"] for _, r in batch_sample_info.iterrows() if r['sample_type'] == "Normal"]
        batch_normalsid = [k for k, _ in batch_sample_info.iterrows() if _['sample_type'] == "Normal"]
        batch_tumors = [r["participant"] for _, r in batch_sample_info.iterrows() if r['sample_type'] == "Tumor"]
        batch_tumorsid = [k for k,_ in batch_sample_info.iterrows() if _['sample_type'] == "Tumor"]
        # create 3 batch-level sample sets: all, tumors, normals
        terra.addToSampleSet(proc_workspace, samplesetid=samplesetnames_all[i], samples=batch_sample_info.index.tolist())
        terra.addToSampleSet(proc_workspace, samplesetid=samplesetnames_tumors[i], samples=batch_tumorsid)
        terra.addToSampleSet(proc_workspace, samplesetid=samplesetnames_normals[i], samples=batch_normalsid)
    return


def create_pair_sets_per_batch(samplesetnames, samplesetnames_pairs, proc_workspace, dict_pairs_per_batch):
    print("Uploading a pair set for each batch...")
    for i, current_batch in enumerate(samplesetnames):
        terra.addToPairSet(proc_workspace, samplesetnames_pairs[i], dict_pairs_per_batch[current_batch])
    return



def create_samplesets_and_pairsets_per_cohort(sample_info, samplesetnames, proc_workspace, cohorts2id_url, newpairs):
    print("Creating and uploading pairsets and samplesets by cohort (iterating over each batch)...")
    cohorts = sheets.get(cohorts2id_url).sheets[0].to_frame()
    cohorts_per_batch = {}
    for i, current_batch in enumerate(samplesetnames):

        # get appropriate subset of the samples for each batch
        batch_sample_info = sample_info[sample_info['batch'] == samplesetnames[i]]
        cohorts_in_batch = []
        # for each batch, make pairsets and samplesets by cohort
        for val in cohorts['ID'].values:
            cohortsamples = batch_sample_info[batch_sample_info["cohorts"] == val].index.tolist()
            tumorsamplesincohort = batch_sample_info[batch_sample_info["cohorts"] == val][batch_sample_info['sample_type']=="Tumor"].index.tolist()
            pairsamples = newpairs[newpairs['case_sample'].isin(tumorsamplesincohort)].index.tolist()
            if len(cohortsamples)>0:
                cohorts_in_batch.append(val)
                terra.addToSampleSet(proc_workspace, val, cohortsamples)

            if len(pairsamples)>0:
                terra.addToPairSet(proc_workspace,val, pairsamples)

        batch_name = samplesetnames[i]
        cohorts_per_batch.update(batch_name = cohorts_in_batch)
    return


# old aggregate function to generate all the cohort and batch specific sammple sets and pair sets
def create_pairsets_and_samplesets(wto, samplesetnames, proc_workspace, sample_info, dict_pairs_per_batch, samplesetnames_all, samplesetnames_tumors, samplesetnames_normals):
    # TODO: Uploading to Terra takes time. Right now, I iterate over cohorts and upload during the iteration. This means that if 2 batches have samples belonging to the same cohort, I upload to that cohort's sample set in Terra 2 times. This is inefficient. We could greatly speed this up (if we run multiple batches at once) by only uploading the cohort-level pairsets and samplesets at the end. That would mean tracking for each cohort all of the new samples (across each batch), then updating Terra at the end.
    # create a pair set for each batch.
    cohorts_per_batch = {}
    for i, current_batch in enumerate(samplesetnames):
        # upload a pair set for each batch
        terra.addToPairSet(proc_workspace, samplesetnames_pairs[i], dict_pairs_per_batch[current_batch])

        # get appropriate subset of the samples for each batch
        batch_sample_info = sample_info[sample_info['batch'] == samplesetnames[i]]
        cohorts_in_batch = []
        cohorts_with_pairs = [] # check: currently not used.
        # for each batch, make pairsets by cohort
        for val in cohorts['ID'].values:
            cohortsamples = batch_sample_info[batch_sample_info["cohorts"] == val].index.tolist()
            tumorsamplesincohort = batch_sample_info[batch_sample_info["cohorts"] == val][batch_sample_info['sample_type']=="Tumor"].index.tolist()
            pairsamples = newpairs[newpairs['case_sample'].isin(tumorsamplesincohort)].index.tolist()
            if len(cohortsamples)>0:
                cohorts_in_batch.append(val)
                terra.addToSampleSet(proc_workspace, val, cohortsamples)

            if len(pairsamples)>0:
                cohorts_with_pairs.append(val)
                terra.addToPairSet(proc_workspace,val, pairsamples)

        batch_name = samplesetnames[i]
        cohorts_per_batch.update(batch_name = cohorts_in_batch)

    print("creating sample sets for each batch...")
    # want to create a sample set for each batch
    for i, current_batch in enumerate(samplesetnames):
        # get appropriate subset of the samples
        batch_sample_info = sample_info[sample_info['batch'] == current_batch]
        # define batch-specific tumors and normals
        batch_normals = [r["participant"] for _, r in batch_sample_info.iterrows() if r['sample_type'] == "Normal"]
        batch_normalsid = [k for k, _ in batch_sample_info.iterrows() if _['sample_type'] == "Normal"]
        batch_tumors = [r["participant"] for _, r in batch_sample_info.iterrows() if r['sample_type'] == "Tumor"]
        batch_tumorsid = [k for k,_ in batch_sample_info.iterrows() if _['sample_type'] == "Tumor"]
        # create 3 batch-level sample sets: all, tumors, normals
        terra.addToSampleSet(proc_workspace, samplesetid=samplesetnames_all[i], samples=batch_sample_info.index.tolist())
        terra.addToSampleSet(proc_workspace, samplesetid=samplesetnames_tumors[i], samples=batch_tumorsid)
        terra.addToSampleSet(proc_workspace, samplesetid=samplesetnames_normals[i], samples=batch_normalsid)

    print("creating sample sets for all samples in workspace, and all normals in workspace...")
    # create sample set for all normals in workspace (will use all combined normals in PoN for mutation calling)
    normalsid.extend([k for k, _ in refsamples.iterrows() if _.sample_type == "Normal"])
    terra.addToSampleSet(proc_workspace, samplesetid="All_normals_TWIST", samples=normalsid)

    # create sample sets for all samples in workspace
    all_samples = wto.get_samples().index.tolist()
    all_samples.remove('NA')
    terra.addToSampleSet(proc_workspace, samplesetid="All_samples_TWIST", samples=all_samples)
