import os
import re
import pandas as pd
import numpy as np
from firecloud import api as firecloud_api
import datetime
import glob

"""
 _____ ___ ____  _____ ____ _     ___  _   _ ____    ____  _   _  ____    _    ____  
|  ___|_ _|  _ \| ____/ ___| |   / _ \| | | |  _ \  / ___|| | | |/ ___|  / \  |  _ \ 
| |_   | || |_) |  _|| |   | |  | | | | | | | | | | \___ \| | | | |  _  / _ \ | |_) |
|  _|  | ||  _ <| |__| |___| |__| |_| | |_| | |_| |  ___) | |_| | |_| |/ ___ \|  _ < 
|_|   |___|_| \_\_____\____|_____\___/ \___/|____/  |____/ \___/ \____/_/   \_\_| \_\
"""

"""
INSTRUCTIONS:
FILES TO UPDATE PRIOR TO RUNNING:
cohort_files/bsp_latest_all_samples_ --> REPLACE WITH LATEST BSP SPREADSHEET
remote_files/ --> ALL THE FILES HERE
walkupseqfiles/ --> ALL FILES HERE AND COMBINE

"""
########################################################
# Sample functions
########################################################
def get_samples(paths_to_batches_info, google_bucket_id, sublist=None):
    """
    Compile samples from multiple batches
    Args: Self-explanatory
        - paths_to_samples_info: .xlsx file containing paths to files containing sample_info
        - sublist: list of tsca_ids to only compile data from certain batches. If None, compile data from all batches.
    Returns: 
        - df with samples from all batches
    """
    paths_to_samples = pd.read_excel(paths_to_batches_info, index_col=0)
    df_list = []

    for tsca_id, paths in paths_to_samples.iterrows():
        if sublist is not None and tsca_id not in sublist:
            continue
        # Make data Firecloud-compatible
        batch_data = prepare_batch_samples_for_metadata_export(paths.path_to_samples_info, tsca_id, google_bucket_id)
        df_list.append(batch_data)

    all_samples = pd.concat(df_list, axis=0)
    
    # Add cohort codes to data
    cohort_formatted_names = pd.read_table('cohort_files/cohort_names_dictionary.txt', header=None, names=['Collection', 'cohort_code'])
    all_samples = pd.merge(all_samples, cohort_formatted_names, on='Collection', how='left')
    return all_samples

def add_cohort_to_old_batches(all_samples):
    """
    Older batches didn't come with the Collection attribute, so it must be added
    Args:
        - all_samples: all samples, without cohort data
    """
    # Retrieve list of samples with corresponding cohort from bsp.broadinstitute.org
    samples_with_cohort = pd.read_excel('cohort_files/bsp_latest_all_samples_TSCA22.xls')
    # Add column for join
    samples_with_cohort['bsp_sample_id_validation'] = samples_with_cohort['Sample ID']

    # FC doesn't accept cohort names with non-alphanumeric characters, so use cohort codes instead
    # Load dictionary of {long cohort name : short cohort code}
    cohort_formatted_names = pd.read_table('cohort_files/cohort_names_dictionary.txt', header=None, names=['Collection', 'cohort_code'])
    # Add cohort codes to samples_with_cohort
    samples_with_cohort = pd.merge(samples_with_cohort, cohort_formatted_names, on='Collection', how='inner')

    # Add cohort data to all samples
    data = pd.merge(all_samples, samples_with_cohort[['bsp_sample_id_validation', 'cohort_code', 'Collection']], \
                         on='bsp_sample_id_validation', \
                         how='left')

    # Merge two `Collection` columns created
    data.loc[pd.isnull(data['cohort_code_x']), 'cohort_code'] = data.loc[pd.isnull(data['cohort_code_x']), 'cohort_code_y']
    data.loc[pd.isnull(data['cohort_code_y']), 'cohort_code'] = data.loc[pd.isnull(data['cohort_code_y']), 'cohort_code_x']
    data.loc[pd.isnull(data['Collection_x']), 'Collection'] = data.loc[pd.isnull(data['Collection_x']), 'Collection_y']
    data.loc[pd.isnull(data['Collection_y']), 'Collection'] = data.loc[pd.isnull(data['Collection_y']), 'Collection_x']
    data = data.drop(['cohort_code_x', 'cohort_code_y'], axis=1)
    data = data.drop(['Collection_x', 'Collection_y'], axis=1)
    return data

def prepare_batch_samples_for_metadata_export(path_to_samples_info, tsca_id, google_bucket_id):
    """Prepare the file to export samples metadata to firecloud
    Args:
        path_id: path to file ending in {}.import_samples.txt
        tsca_id: TSCAXX
        google_bucket_id: id of google bucket ('gs://google_bucket_id')
    Returns:
        pd.DF of data ready for export
    """
    # export raw data
    data = pd.read_table(path_to_samples_info)
    # Rename columns to match firecloud requirements
    data = data.rename(columns={'sample_id':'entity:sample_id', 'individual_id':'participant_id'})
    # Locations of BAM files in google bucket
    path_in_bucket_full = "gs://%s/seq_data/%s" % (google_bucket_id, tsca_id)
    # Extract bam filename
    data['bam_filename'] = data.apply(lambda row: row['clean_bam_file_capture'].split('/')[-1], axis=1)
    # Create bai filename (change extension on .bam file)
    data['bai_filename'] = data.apply(lambda row: "%s%s" %(row['bam_filename'][:-3], 'bai'), axis=1)
    # Change BAM path from xchip to Google cloud
    data['clean_bam_file_capture'] = \
        data.apply( lambda row: "%s/%s/%s" \
                   %(path_in_bucket_full, row['external_id_validation'], row['bam_filename']), axis=1)
    # Add location of .bai file 
    data['clean_bai_file_capture'] = \
        data.apply( lambda row: "%s/%s/%s" \
                   %(path_in_bucket_full, row['external_id_validation'], row['bai_filename']), axis=1)
    # Add TSCA ID
    data['tsca_id'] = tsca_id

    # Reorganize columns (entity:sample_id at the beginning)
    columns = ['entity:sample_id'] + [col for col in data if col != 'entity:sample_id']
    data = data[columns]
    return data

def save_and_upload_samples(data, namespace, workspace, tsca_id):
    """Create FC uploading file and upload patients to FC
    Args:
        - data: participants df
    Writes: 
        - {tsca_id}/fc_upload_patients_tsca_{tsca_id}.txt
    """
    os.system('mkdir -p %s'%tsca_id)
    filename = '%s/fc_upload_samples_tsca_%s.txt' % (tsca_id, tsca_id)
    data.to_csv(filename, sep='\t', index=False)
    res = upload_entities_from_tsv(namespace, workspace, filename)
    return res

########################################################
# Pairs functions
########################################################
def create_pairs_list(all_samples):
    """Creates DF with pairs for firecloud metadata export.
    Args:
        - all_samples: all samples.
    """
    dfs = []
    # Find match normals for tumor samples only
    tumor_samples = all_samples[all_samples.sample_type=="Tumor"]
    i = 0
    for index, row in tumor_samples.iterrows():
        # Find all samples from same individual (same individual_id, different sample_id)
        patient_samples = all_samples[ (all_samples['participant_id'] == row['participant_id']) \
                                          & (all_samples['entity:sample_id'] != row['entity:sample_id']) ]

        # NOTE: If more than one match tumor tissue or match normal found, select first one found.
        # The match normal is used to compute allelic fractions in Mutect2, so for now we ignore the conditions it workspaces grown in.

        ######## Match normal: Add match normal
        match_normal = patient_samples[ patient_samples['sample_type'] == "Normal"]
        #   > No match normal found
        if match_normal.empty: 
            control_sample_id = "NA"
            control_sample_tsca_id = "NA"
        #   > Match normal found
        elif match_normal.shape[0] > 0:
            match_normal = match_normal.iloc[0]
            control_sample_id = match_normal['entity:sample_id']
            control_sample_tsca_id = match_normal['tsca_id']
        
        # Create DF with Tumor/Normal pair set
        pair_id = "%s_%s_TN" % (row['entity:sample_id'], control_sample_id)
        df_dict = {'entity:pair_id': pair_id, 'case_sample_id': row['entity:sample_id'], \
                    'control_sample_id': control_sample_id, 'participant_id': row['participant_id'], 'match_type': 'tumor_normal', \
                    'case_sample_tsca_id': row['tsca_id'], 'control_sample_tsca_id': control_sample_tsca_id}
        dfs.append(pd.DataFrame(df_dict, index=[i], columns=df_dict.keys()))
        i+=1
        
        ######## Tumor tissue: Add primary tumor tissue
        match_primary_tumor = patient_samples[ ( patient_samples['external_id_validation'].str.contains('primary|prim|tissue|tiss|Primary|Tissue') ) & \
        										(patient_samples['sample_type'] == "Tumor")]
        #    > No primary tumor tissue found
        if match_primary_tumor.empty:
            control_sample_id = "NA"
            control_sample_tsca_id = "NA"
        #    > Sample itself is a primary tumor tissue
        elif any(substring in row['external_id_validation'] for substring in ['primary', 'prim', 'tissue', 'tiss', 'Primary', 'Tissue']):
            control_sample_id = "NA"
            control_sample_tsca_id = "NA"
        #    > Tumor tissue found
        elif match_primary_tumor.shape[0] > 0:
            match_primary_tumor = match_primary_tumor.iloc[0]
            control_sample_id = match_primary_tumor['entity:sample_id']
            control_sample_tsca_id = match_primary_tumor['tsca_id']
        
        # Create DF with Tumor/Primary pair set
        pair_id = "%s_%s_TP" % (row['entity:sample_id'], control_sample_id)
        df_dict = {'entity:pair_id': pair_id, 'case_sample_id': row['entity:sample_id'], \
                    'control_sample_id': control_sample_id, 'participant_id': row['participant_id'], 'match_type': 'tumor_primary', \
                    'case_sample_tsca_id': row['tsca_id'], 'control_sample_tsca_id': control_sample_tsca_id}
        dfs.append(pd.DataFrame(df_dict, index=[i], columns=df_dict.keys()))
        i+=1
   
    return pd.concat(dfs, axis=0)

def save_and_upload_pairs(namespace, workspace, pairs, blacklist=[]):
    """Updates pairs to firecloud. 
    NOTE: All pairs need to be updated with every new batch,
    as it may contain match normals or primary matches for previous batches.
    Args:
        - pairs: df of all pairs, as created by create_pairs_list
    Returns: 
        - res: json response from http request
    Creates: 
        - ./Pairs/fc_upload_pairs.txt file
    """
    pairs = pairs[ ~pairs['case_sample_id'].isin(blacklist)]
    os.system('mkdir -p Pairs')
    filename = './Pairs/fc_upload_pairs.txt'
    pairs.to_csv(filename, '\t', index=False)
    res = upload_entities_from_tsv(namespace, workspace, filename)
    return res

########################################################
# Pair set functions
########################################################
def prepare_batch_pairsets_for_metadata_export(all_samples, pairs, tsca_id):
    """Creates pair sets for FC export, both tumor-normal and tumor-primary pairs
    Args:
        - all_samples: all samples
        - pairs: pairs df as created by create_pairs_list
    Returns:
        - Tuple of dfs: (tumor-normal, tumor-primary)
    """
    tn_pairs = pairs[(pairs['match_type'] == "tumor_normal") & (pairs['case_sample_tsca_id']==tsca_id)]
    tp_pairs = pairs[(pairs['match_type'] == "tumor_primary") & (pairs['case_sample_tsca_id']==tsca_id)]

    tn_pairsets = pd.merge(tn_pairs, all_samples[['entity:sample_id', 'tsca_id']], \
                            left_on='case_sample_id', right_on='entity:sample_id', \
                            how='inner')[['tsca_id', 'entity:pair_id']] \
                            .rename(columns={'tsca_id': 'membership:pair_set_id', 'entity:pair_id': 'pair_id'})

    tp_pairsets = pd.merge(tp_pairs, all_samples[['entity:sample_id', 'tsca_id']], \
                            left_on='case_sample_id', right_on='entity:sample_id', \
                            how='inner')[['tsca_id', 'entity:pair_id']] \
                            .rename(columns={'tsca_id': 'membership:pair_set_id', 'entity:pair_id': 'pair_id'})

    # Append _TN/_TP to the end of the tumor-normal/tumor-primary pair set ids
    tn_pairsets['membership:pair_set_id'] = tn_pairsets['membership:pair_set_id'].apply(lambda x: "%s_TN"%x)
    tp_pairsets['membership:pair_set_id'] = tp_pairsets['membership:pair_set_id'].apply(lambda x: "%s_TP"%x)
    
    return (tn_pairsets, tp_pairsets)

def prepare_cumulative_pairsets_for_metadata_export(pairs, name):
	"""Creates pair sets for FC export, both tumor-normal and tumor-primary pairs. 
	Pair sets are cumulative, not batch-based.
	Args:
		- pairs: pairs df as created by create_pairs_list
		- name: name of cumulative pairset, usually TSCA# of current tsca_id
	Returns:
	    - Tuple of dfs: (tumor-normal, tumor-primary) of cumulative pairs
	"""
	tn_pairsets = pairs.loc[pairs['match_type'] == "tumor_normal", ["entity:pair_id"]].rename(columns={"entity:pair_id":"pair"})
	tp_pairsets = pairs.loc[pairs['match_type'] == "tumor_primary", ["entity:pair_id"]].rename(columns={"entity:pair_id":"pair"})
	tn_pairsets['membership:pair_set_id'] = "Cum_TN_%s_all"%name
	tp_pairsets['membership:pair_set_id'] = "Cum_TP_%s_all"%name
		
	return tn_pairsets[['membership:pair_set_id', 'pair']], tp_pairsets[['membership:pair_set_id', 'pair']]
	
def upload_pairsets(namespace, workspace, pairsets, pairset_type):
    os.system('mkdir -p Pairs')
    filename = './Pairs/fc_upload_pairsets_%s.txt'%pairset_type
    pairsets.to_csv(filename, '\t', index=False)
    res = upload_entities_from_tsv(namespace, workspace, filename)
    return res

###############################################
# Participant functions
###############################################
def prepare_participants_for_metadata_export(path_to_samples_info, tsca_id):
    """Create participant entities DF for Firecloud. 
    Participants need to exist before you can upload their respective samples
    """    
    raw = pd.read_table(path_to_samples_info)
    print( "%d Participants in this batch" % raw['individual_id'].unique().shape[0] )
    # Data to upload
    data = pd.DataFrame(raw.individual_id.drop_duplicates()).rename(columns={'individual_id':'entity:participant_id'})
    return data

def save_and_upload_participants(data, namespace, workspace, tsca_id):
    """Create FC uploading file and upload patients to FC
    Args:
        - data: participants df
    Writes: 
        - {tsca_id}/fc_upload_patients_tsca_{tsca_id}.txt
    """
    os.system('mkdir -p %s'%tsca_id)
    filename = './%s/fc_upload_patients_%s.txt' % (tsca_id, tsca_id)
    data.to_csv(filename, '\t', index=False)
    res = upload_entities_from_tsv(namespace, workspace, filename)
    return res

########################################################
# Sample set functions
########################################################
def prepare_batch_sample_set_for_metadata_export(path, tsca_id):
    """Create dfs and write files with sample_set metadata for Firecloud.
    Three sample sets are created: all samples, normal samples, and tumor samples in batch
    A sample for a given batch 
    Args:
        - path: path to file ending in {}.import_samples.txt
        - tsca_id: batch tsca id
    Returns:
        - (all_samples, tumor_samples, normal_samples)
    Writes:
        - files for all_samples, tumor_samples, normal_samples
    """
    raw = pd.read_table(path)
    print( "%d Samples in this batch" % raw.shape[0] )

    # Create dfs to upload
    all_samples = pd.concat([pd.DataFrame(index=raw.index, columns=['membership:sample_set_id'], data=tsca_id), \
                      raw[ ['sample_id', 'sample_type'] ]], axis=1)


    tumors  = all_samples.loc[ all_samples['sample_type'] == "Tumor", ['membership:sample_set_id', 'sample_id'] ]
    tumors.loc[: , 'membership:sample_set_id'] = "%s_T"%tsca_id
    
    normals = all_samples.loc[ all_samples['sample_type'] == "Normal", ['membership:sample_set_id', 'sample_id'] ]
    normals.loc[: , 'membership:sample_set_id'] = "%s_N"%tsca_id

    all_samples = all_samples.drop('sample_type', axis=1)
    return (all_samples, tumors, normals)

def filter_existing_samples(df, sample_id_colname, remote_samples):
    """If sample set is uploaded more than once, samples will be duplicated, which is undesirable behavior.
    Therefore, need to download samples existing remotely (in FC), and remove them from sample set to be uploaded.
    Args:
        - df: dataframe to filter
        - sample_id_colname: name of column containing sample_id
    """
    remote_sample_ids = remote_samples['entity:sample_id'].tolist()
    df_clean = df[~df[sample_id_colname].isin(remote_sample_ids)]
    return df_clean

def save_and_upload_batch_sample_sets(batch_samples, batch_tumors, batch_normals, tsca_id, namespace, workspace):
    """Create FC uploading file and upload patients to FC
    """
    # Save to file
    os.system('mkdir -p %s'%tsca_id)
    batch_samples_filename = './%s/fc_upload_sample_set_tsca_%s.txt' % (tsca_id, tsca_id)
    batch_tumors_filename = './%s/fc_upload_sample_set_tsca_%s_tumors.txt' % (tsca_id, tsca_id)
    batch_normals_filename = './%s/fc_upload_sample_set_tsca_%s_normals.txt' % (tsca_id, tsca_id)
    
    batch_samples.to_csv(batch_samples_filename , sep="\t", index=False )
    batch_tumors.to_csv(batch_tumors_filename , sep="\t", index=False )
    batch_normals.to_csv(batch_normals_filename , sep="\t", index=False )

    r1 = upload_entities_from_tsv(namespace, workspace, batch_samples_filename)
    r2 = upload_entities_from_tsv(namespace, workspace, batch_tumors_filename)
    r3 = upload_entities_from_tsv(namespace, workspace, batch_normals_filename)
    return (r1, r2, r3)

########################################################
# Cohort functions
########################################################
def prepare_cohorts_for_metadata_export(all_samples, blacklist=[]):
    """Creates sample sets corresponding to cohorts for Firecloud export.
    Args:
        - all_samples: with cohort
    Returns: 
        - metadata with samples and cohorts they belong to for FC upload
    """
    # Prepare for FC export format
    data = all_samples.rename(columns={'cohort_code': 'membership:sample_set_id', 'entity:sample_id': 'sample_id'})
    data_clean = data[['membership:sample_set_id', 'sample_id']]

    # Remove blacklist
    data_clean = data_clean[ ~data_clean['sample_id'].isin(blacklist)]

    return data_clean

def save_and_upload_cohorts(data, latest_tsca_id, namespace, workspace):
    """Save and upload cohort metadata to FC
    """
    # Write FC import file 
    filename = 'cohort_files/fc_upload_sample_set_cohorts_%s.txt'%latest_tsca_id
    data.to_csv(filename, index=False, sep="\t")
    res = upload_entities_from_tsv(namespace, workspace, filename)
    return res

def prepare_cohort_pairsets_for_metadata_exports(latest_tsca_id, pairs, all_samples, blacklist=[]):
    """Create DF with cohort pairsets, used to create cohort reports for SNVs.
    Args:
        - pairs: pair list
        - all_samples: all samples with cohort
    Returns: 
        - DF with [cohort_code, pair_id]
    """
    # Create list of pairs
    clean_pairs_list = pairs[ ~pairs['case_sample_id'].isin(blacklist)]

    # Keep only tumor-normal pairs, as only these pairs are used for cohort reports
    clean_pairs_list = clean_pairs_list[clean_pairs_list["match_type"]=="tumor_normal"]
    
    # Add cohorts to pairs
    pairs_with_cohort = pd.merge(clean_pairs_list, all_samples[['entity:sample_id', 'cohort_code']], \
             left_on='case_sample_id', right_on='entity:sample_id')
    # Prepare DF for FC export
    pairs_with_cohort_clean = pairs_with_cohort[['cohort_code', 'entity:pair_id']] \
        .rename(columns={'entity:pair_id': 'pair_id', 'cohort_code': 'membership:pair_set_id'})

    return pairs_with_cohort_clean

def save_and_upload_cohort_pairsets(namespace, workspace, pairsets):
    os.system('mkdir -p Pairs')
    filename = './Pairs/fc_upload_cohort_pairsets.txt'
    pairsets.to_csv(filename, '\t', index=False)
    res = upload_entities_from_tsv(namespace, workspace, filename)
    return res

def save_and_upload_cohort_all_samples(all_samples, name, namespace, workspace, blacklist=[]):
    """Create and upload cohort all samples
    Args: 
        - Self-explanatory
    """
    df = all_samples[['entity:sample_id']].rename(columns={'entity:sample_id': 'sample_id'})
    df['membership:sample_set_id'] = name

    # Re-arrange columns
    cols = ['membership:sample_set_id', 'sample_id']
    df = df[cols]

    # Blacklist
    df = df[ ~df['sample_id'].isin(blacklist) ]
    df.to_csv('all_samples/fc_upload_%s.txt'%name, index=None, sep="\t")
    res = upload_entities_from_tsv(namespace, workspace, 'all_samples/fc_upload_%s.txt'%name)
    return res

def save_and_upload_cohort_all_tumors(all_samples, name, namespace, workspace, blacklist=[]):
    """Create and upload cohort of all tumor samples across all batches
    Args: 
        - Self-explanatory
        - name: cohort name (usually Cum_Tumors_{LATEST_TSCA_ID}_all)
        - paths_to_samples_info: .xlsx file containing paths to files containing sample_info
    """
    tumor_samples = all_samples[all_samples.sample_type == "Tumor"]

    # Prepare column names
    df = tumor_samples[['entity:sample_id']].rename(columns={'entity:sample_id': 'sample_id'})
    df['membership:sample_set_id'] = name

    # Re-arrange columns
    cols = ['membership:sample_set_id', 'sample_id']
    df = df[cols]

    # Blacklist
    df = df[ ~df['sample_id'].isin(blacklist) ]
    df.to_csv('tumor_samples/fc_upload_%s.txt'%name, index=None, sep="\t")
    res = upload_entities_from_tsv(namespace, workspace, 'tumor_samples/fc_upload_%s.txt'%name)
    return res

########################################################
# PoN Functions
########################################################
def create_panel_of_normals_advanced(tsca_id, all_samples, num_normals_per_cohort_involved = 3, batches_to_pick_from = []):
    """Create a panel of normals for batch with @tsca_id.
    Add @num_normals_per_cohort_involved normal samples from all the cohorts involved.
    You can restrict the batches you pool normals from to the list of batches in @batches_to_pick_from.
    """
    # Get all samples
    batch_samples = all_samples[all_samples['tsca_id']==tsca_id]
    # Batch normals    
    batch_normals = batch_samples[batch_samples['sample_type']=="Normal"]
    # Number of normals in batch
    num_normals_from_batch = batch_normals.shape[0]
    # Cohorts of samples in batch
    cohorts_involved = batch_samples['cohort_code'].unique()

    # Only select normals from the restricted batches
    restricted_normals = all_samples[(all_samples['tsca_id'].isin(batches_to_pick_from)) & (all_samples['sample_type']=="Normal")]
    # Merge all normals from cohorts involved
    cohorts_normals_lst = []
    for cohort_involved in cohorts_involved:
        cohort_normals = restricted_normals[(restricted_normals['cohort_code']==cohort_involved)] \
                            .iloc[:num_normals_per_cohort_involved]
        cohorts_normals_lst.append(cohort_normals)

    cohorts_normals = pd.concat(cohorts_normals_lst)

    # Final PoN: batch normals + normals cohorts involved
    final_pon = pd.concat([batch_normals, cohorts_normals])
    num_normals = final_pon.shape[0]
    final_pon_name = "PoN_%s_%s_batch_normals_%s_normals_per_cohort_%s_total"\
                        %(tsca_id, num_normals_from_batch, num_normals_per_cohort_involved, num_normals)

    # Prepare for FC format
    final_pon['membership:sample_set_id'] = final_pon_name
    final_pon['sample_id'] = final_pon['entity:sample_id']
    final_pon = final_pon[['membership:sample_set_id', 'sample_id']]
    return final_pon, final_pon_name

def create_panel_of_normals_from_small_batch(tsca_id, all_samples, num_normals_from_batch = -1, num_normals_per_cohort_involved = 3, num_normals = 25):
    """Create panel of normals with samples from a given small batch.
    Small batches are 24-48 samples (instead of the traditional 96).
    The main difference is that there are only 2-4 normals per batch, so we must include normals from all the cohorts
    involved in this batch.
    Args: 
        - all_samples: with cohort
        - N: (int) number of normals from the batch to include in the PoN
        - tsca_id: tsca to build PoN for
        - num_normals_from_batch: number of normal samples from batch (-1 for all)
        - num_normals_per_cohort_involved: number of normal samples from every cohort involved in batch
        - num_normals: total number of normals
    """
    # Get all samples
    batch_samples = all_samples[all_samples['tsca_id']==tsca_id]
    # Batch normals    
    batch_normals = batch_samples[batch_samples['sample_type']=="Normal"].iloc[:num_normals_from_batch]
    # Cohorts of samples in batch
    cohorts_involved = batch_samples['cohort_code'].unique()
    
    # Merge all normals from cohorts involved
    cohorts_normals_lst = []
    for cohort_involved in cohorts_involved:
        cohort_normals = all_samples[(all_samples['cohort_code']==cohort_involved) & (all_samples['sample_type']=="Normal")] \
                            .iloc[:num_normals_per_cohort_involved]
        cohorts_normals_lst.append(cohort_normals)

    cohorts_normals = pd.concat(cohorts_normals_lst)

    # Number of normals necessary to complete a total of @num_normals in PoN
    num_missing_normals = num_normals - cohorts_normals.shape[0] - batch_normals.shape[0]
    # If missing normals, select at random from the rest of the samples
    if num_missing_normals > 0:
        random_normals = all_samples[all_samples['sample_type']=="Normal"].sample(n=num_missing_normals)

    # Final PoN: batch normals + normals cohorts involved + random normals to complete if necessary
    final_pon = pd.concat([batch_normals, cohorts_normals, random_normals])
    # If num_normals_from_batch is -1, using all normals in batch
    if num_normals_from_batch == -1:
        num_normals_from_batch = "all"
    final_pon_name = "PoN_%s_%s_batch_normals_%s_normals_per_cohort_%s_total"\
                        %(tsca_id, num_normals_from_batch, num_normals_per_cohort_involved, num_normals)

    # Prepare for FC format
    final_pon['membership:sample_set_id'] = final_pon_name
    final_pon['sample_id'] = final_pon['entity:sample_id']
    final_pon = final_pon[['membership:sample_set_id', 'sample_id']]
    return final_pon, final_pon_name

def create_panel_of_normals(paths, N, name):
    """Create panel of normals sample set for Firecloud from multiple TSCA batches.
    Randomly select N samples from samples present in files listed in paths.
    Args:
        paths: (list) paths to file ending in {}.import_samples.txt
        N: (int) number of samples in panel of normals
        name: (string) name of Panel of Normals
    """
    dfs = [ pd.read_table(paths[0]) ]
    for i, path in enumerate(paths[1:]):
        df_to_concat = pd.read_table(path)
        dfs.append(df_to_concat)
    df = pd.concat(dfs, axis=0)
    # Shuffle samples to pick from all batches
    df = df.sample(frac=1).reset_index(drop=True)
    normals = df[df.sample_type=="Normal"][:N]['sample_id']
    if N==-1: print ("Creating panel of %d normals" %normals.shape[0])
    else: print ("Creating panel of %d normals" %N)
    
    # Compile data
    data = pd.concat([pd.DataFrame(index=normals.index, columns=['membership:sample_set_id'], data=name), \
                        normals], axis=1)

    return data
        
def upload_pon(pon_df, pon_name, namespace, workspace):
    """Upload PoN to FC
    Args:
        - pon_df: dataframe with normal samples in PoN
        - pon_name: name of PoN
    """
    os.system('mkdir -p PoNs')
    filename = './PoNs/fc_upload_PoN_%s.txt' % (pon_name)
    pon_df.to_csv(filename, '\t', index=False)
    res = upload_entities_from_tsv(namespace, workspace, 'PoNs/fc_upload_PoN_%s.txt'%pon_name)
    return res

################################################
# Helper Functions
###############################################
def upload_entities_from_tsv(namespace, workspace, entities_tsv_file):
    """Upload entities from tsv file
    Args: 
        Self-explanatory
        entities_tsv_file: path to tsv file
    Returns: 
        HTTP Response
    """
    res = firecloud_api.upload_entities_tsv(namespace, workspace, entities_tsv=entities_tsv_file)
    return res

def delete_pair(namespace, workspace, pair_id):
    """Delete pair from workspace/namespace
    """
    body = [{"entityType": "pair", "entityName": pair_id}]
    res = firecloud_api.delete_entities(namespace, workspace, body)
    return res

def delete_pair_set(namespace, workspace, pair_set_id):
    """Delete pair set from workspace/namespace
    Args: 
        Self-explanatory
    Returns: 
        HTTP Response
    """
    body = [{"entityType": "pair_set", "entityName": pair_set_id}]
    res = firecloud_api.delete_entities(namespace, workspace, body)
    return res

def delete_sample(namespace, workspace, sample_id):
    """Delete sample from workspace/namespace
    Args: 
        Self-explanatory
    Returns: 
        HTTP Response
    """
    body = [{"entityType": "sample", "entityName": sample_id}]
    res = firecloud_api.delete_entities(namespace, workspace, body)
    return res

def delete_sample_set(namespace, workspace, sample_set_id):
    """Delete sample set from workspace/namespace
    Args: 
        Self-explanatory
    Returns: 
        HTTP Response
    """
    body = [{"entityType": "sample_set", "entityName": sample_set_id}]
    res = firecloud_api.delete_entities(namespace, workspace, body)
    return res

def delete_participant(namespace, workspace, participant_id):
    """Delete participant from workspace/namespace
    Args: 
        Self-explanatory
    Returns: 
        HTTP Response
    """
    body = [{"entityType": "participant", "entityName": participant_id}]
    res = firecloud_api.delete_entities(namespace, workspace, body)
    return res

def download_remote_samples(namespace, workspace):
    """Download remote samples from Firecloud
    Writes:
        - remote_samples.txt: samples in FC
    """
    res = firecloud_api.get_entities_tsv(namespace, workspace, "sample")
    with open('remote_files/remote_samples.txt', 'w') as outfile:
        outfile.write(res.text)
    return

def merge_walkupseq_files(latest_tsca_id):
    """Merge all walkupseq files in the walkupseq directory
    These files have metadata (primary disease, media type) to add to sample metadata
    """
    paths = glob.glob('walkupseq_files/*sample_info*')

    dfs = []
    for f in paths:
        tmp = pd.read_table(f, encoding='latin1')
        dfs.append(tmp)

    df = pd.concat(dfs, axis=0)
    df.to_csv('walkupseq_files/walkupseq_all_combined_%s.txt'%latest_tsca_id, sep="\t", index=None)
    return df

def remove_blacklisted_samples(all_samples):
    """Remove blacklisted samples
    """
    # Blacklisted samples
    blacklisted = pd.read_table("samples_blacklist.txt", header=None, names=["entity:sample_id"])
    blacklisted_sample_ids = blacklisted["entity:sample_id"].values.tolist()
    all_samples = all_samples[~all_samples["entity:sample_id"].isin(blacklisted_sample_ids)]
    return all_samples

################################################
### MAIN
###############################################

################################# TODO: #####################################
# UPDATE REMOTE FILES WITH TSCA JUST RUN AFTER THE FIRST TIME BOTH PROGRAMS HAVE SUCCEEDED
# (upload_data_new_batch and update_cohorts)
#############################################################################
# Based on results from commands (r4*, r5, r6)
# For example, if running TSCA 23, tsca_id = "TSCA23" and latest_tsca_id = "TSCA22".
# Run upload_data_new_batch

def upload_data_new_batch(tsca_id, latest_tsca_id, paths_to_batches_info, namespace, workspace, google_bucket_id):
    """
    Upload data for new batch.
    This function only needs to be called once for a new batch.
    Only run once, otherwise you will create duplicate samples in the sample sets
    Args:
        - tsca_id: id of TSCA for which data is to be uploaded
        - latest_tsca_id: id of latest TSCA available
        - paths_to_batches_info: file with paths to .import_samples.txt files (paths_to_batches_info.xlsx)
        - Rest are self-explanatory
    """
    paths_to_batches_info_df = pd.read_excel(paths_to_batches_info, index_col=0)
    path_to_samples_info  = paths_to_batches_info_df.loc[tsca_id, 'path_to_samples_info']

    # DF of remote [sample < > sample set ]
    remote_sample_sets = pd.read_table('remote_files/sample_set_membership_%s.tsv'%latest_tsca_id)
    # DF of remote [pair < > pair set]
    remote_pair_sets = pd.read_table('remote_files/pair_set_membership_%s.tsv'%latest_tsca_id)

    all_samples = get_samples(paths_to_batches_info, google_bucket_id)
    # Add cohorts for older batches
    all_samples = add_cohort_to_old_batches(all_samples)

    ##### Remove blacklisted samples ##
    # Blacklisted samples
    blacklisted = pd.read_table("samples_blacklist.txt", header=None, names=["entity:sample_id"])
    blacklisted_sample_ids = blacklisted["entity:sample_id"].values.tolist()
    all_samples = all_samples[~all_samples["entity:sample_id"].isin(blacklisted_sample_ids)]

    ########## Participants ##########
    print("Uploading participants...")  
    participants = prepare_participants_for_metadata_export(path_to_samples_info, tsca_id)
    r1 = save_and_upload_participants(participants, namespace, workspace, tsca_id)
    ##################################

    ##########  Samples  ############
    print("Uploading samples...")
    batch_samples = prepare_batch_samples_for_metadata_export(path_to_samples_info, tsca_id, google_bucket_id)
    r2 = save_and_upload_samples(batch_samples, namespace, workspace, tsca_id)
    #################################

    ##########   Pairs   #############
    print("Uploading pairs...")
    pairs = create_pairs_list(all_samples)
    r3 = save_and_upload_pairs(namespace, workspace, pairs)
    #################################

    ##########  Sample Sets  #########
    print("Uploading sample sets...")
    batch_sample_set, batch_tumor_set, batch_normal_set = prepare_batch_sample_set_for_metadata_export(path_to_samples_info, tsca_id)
    # Remove the samples that have already been uploaded 
    uploaded_sample_ids = remote_sample_sets['sample'].tolist()
    batch_sample_set_clean = batch_sample_set[~batch_sample_set['sample_id'].isin(uploaded_sample_ids)]
    batch_tumor_set_clean = batch_tumor_set[~batch_tumor_set['sample_id'].isin(uploaded_sample_ids)]
    batch_normal_set_clean = batch_normal_set[~batch_normal_set['sample_id'].isin(uploaded_sample_ids)]
    r4a, r4b, r4c = save_and_upload_batch_sample_sets(batch_sample_set_clean, batch_tumor_set_clean, batch_normal_set_clean, tsca_id, namespace, workspace)
    #################################

    ##########  PoNs  ###############
    print("Uploading PoNs...")
    
    # Number of latest tsca id
    latest_tsca_id_int = int(re.findall('\d+', latest_tsca_id )[0])
    # Array with list of all previous TSCA ids
    previous_tsca_ids = ["TSCA%s"%i for i in np.arange(14, latest_tsca_id_int+1)]
    previous_tsca_ids.insert(0, "TSCA1213")

    pon, name = create_panel_of_normals_advanced(tsca_id, all_samples,\
                    num_normals_per_cohort_involved = 3, \
                    batches_to_pick_from = previous_tsca_ids)

    # Only upload PoN if it hasn't been uploaded already
    if not name in remote_sample_sets['membership:sample_set_id'].unique().tolist():
        r5 = upload_pon(pon, name, namespace, workspace)    
    else: 
        print("PoN already exists...")
        r5 = {}
    #################################
    
    ##########  Pair Set  ###########
    print("Uploading pair sets...")
    # Upload cumulative pair sets
    tn_cum_pairsets, tp_cum_pairsets = prepare_cumulative_pairsets_for_metadata_export(pairs, tsca_id)
    r6 = upload_pairsets(namespace, workspace, tn_cum_pairsets, "TN")
    r7 = upload_pairsets(namespace, workspace, tp_cum_pairsets, "TP")

    # Batch pair sets
    tn_pairsets, tp_pairsets = prepare_batch_pairsets_for_metadata_export(all_samples, pairs, tsca_id)
    uploaded_pair_ids = remote_pair_sets['pair'].tolist()
    tn_pairsets_clean = tn_pairsets[~tn_pairsets['pair_id'].isin(uploaded_pair_ids)]
    tp_pairsets_clean = tp_pairsets[~tp_pairsets['pair_id'].isin(uploaded_pair_ids)]

    r8 = upload_pairsets(namespace, workspace, tn_pairsets_clean, "TN")
    r9 = upload_pairsets(namespace, workspace, tp_pairsets_clean, "TP")
    #################################

    return (r1, r2, r3, r4a, r4b, r4c, r5, r6, r7, r8, r9)

def update_cohorts(tsca_id, latest_tsca_id, paths_to_batches_info, namespace, workspace, google_bucket_id):
    """
    Update cohorts (sample sets that span multiple batches)
    This function needs to be called once for a new batch. Before updating a cohort sample set, it removes samples that 
    already belong to that cohort remotely, because if we don't they will be duplicated.
    Args: 
        - latest_tsca_id: id of latest TSCA available
        - paths_to_batches_info: file with paths to .import_samples.txt files (paths_to_batches_info.xlsx)
        - Rest are self-explanatory
    """
    # Pre-requisites
    all_samples = get_samples(paths_to_batches_info, google_bucket_id)
    # Add cohorts for older batches
    all_samples = add_cohort_to_old_batches(all_samples)
    # Remove blacklisted samples
    all_samples = remove_blacklisted_samples(all_samples)

    pairs = create_pairs_list(all_samples)

    # DF of remote samples
    remote_samples = pd.read_table('remote_files/remote_samples_%s.txt'%latest_tsca_id)
    # DF of remote [sample < > sample set ]
    remote_sample_sets = pd.read_table('remote_files/sample_set_membership_%s.tsv'%latest_tsca_id)
    # DF of remote [pair < > pair set]
    remote_pair_sets = pd.read_table('remote_files/pair_set_membership_%s.tsv'%latest_tsca_id)

    #### UPDATE COHORT SAMPLE SETS
    # DF of [samples < > sample set] to be updated
    cohorts = prepare_cohorts_for_metadata_export(all_samples)
    # Remove the samples that already belong to the cohort in FC 
    sample_ids_in_cohort = remote_sample_sets['sample'].tolist()
    cohorts_clean = cohorts[~cohorts['sample_id'].isin(sample_ids_in_cohort)]
    r1 = save_and_upload_cohorts(cohorts_clean, latest_tsca_id, namespace, workspace)

    #### UPDATE COHORT PAIR SETS
    #  Retrieve cohort pairsets
    cohort_pairsets = prepare_cohort_pairsets_for_metadata_exports(latest_tsca_id, pairs, all_samples, blacklist=[])
    # Remove the pairs that already belong to the cohort in FC
    pair_ids_in_cohort = remote_pair_sets['pair'].tolist()
    cohort_pairsets_clean = cohort_pairsets[~cohort_pairsets['pair_id'].isin(pair_ids_in_cohort)]
    r2 = save_and_upload_cohort_pairsets(namespace, workspace, cohort_pairsets_clean)

    # Remove samples that already exist in FC
    # remote_sample_ids = remote_samples['entity:sample_id'].tolist()
    # all_samples_clean = all_samples[~all_samples['entity:sample_id'].isin(remote_sample_ids)]
    r3 = save_and_upload_cohort_all_samples(all_samples, "Cum_%s_all"%tsca_id, namespace, workspace, blacklist=[])
    r4 = save_and_upload_cohort_all_tumors(all_samples, "Cum_Tumors_%s_all"%tsca_id, namespace, workspace, blacklist=[])

    ### Create cumulative PoN (all batches)
    pon_name = 'Cum_PoN_%s_all'%tsca_id
    paths_to_batches_info_df = pd.read_excel(paths_to_batches_info, index_col=0)
    cumulative_pon = create_panel_of_normals(paths_to_batches_info_df.path_to_samples_info.tolist(), -1, pon_name)
    r5 = upload_pon(cumulative_pon, pon_name, namespace, workspace)
    
    return (r1, r2, r3, r4, r5)

