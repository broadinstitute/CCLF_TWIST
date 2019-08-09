from glob import glob
import os
import re
import pandas as pd
import numpy as np
import dalmatian

##################################################################
### Prepare data (BAM files) for export to Firecloud
### M. Cadosch
### 12/17
##################################################################

##################################################################
# EDIT: Params
##################################################################
# Should be in email
sn_id = "SN0169772"
site = "HT33MBCX2"
tsca_id = "TSCA45"
date = "201904"
external_sheet_path = "/xchip/clf/seq_data/walkup_seq_sample_info/TSCA45 External Sample ID.xlsx"
plate_prefix = "1"
##################################################################

TSCA_version = "TSCA Rapid Cancer Detection Panel v2"
picard_aggregation_type_validation = "PCR"

target_dir_path = "/xchip/clf/seq_data/process_for_fc/%s__%s_%s"%(tsca_id.lower(), date, sn_id)

##################################################################
### Helper functions
##################################################################
def format_position(position):
    """Format position, e.g. A01 -> A1
    """
    return re.sub("([A-Z])0?([0-9]+)", r"\1\2", position)

def get_new_sample_external_id(new_sample_base_external_id, duplicates_external_ids):
    """
    Given the base external id of the new sample @new_sample_base_external_id, 
    and the list of current samples with the same base external id @duplicates_external_ids, 
    returns the formatted new sample external id, with the max index, to ensure no duplicates.
    Params:
        - new_sample_base_external_id: external id of duplicated sample to add
        - duplicates_external_ids: external ids of existing duplicate samples
    Returns:
        - new_sample_external_id: formatted sample id (if A already exists, return A_2, A_3, etc.)
    """
    # Index of duplicate external ids
    duplicates_indices = [ re.sub(r"(%s_?)(.*)"%new_sample_base_external_id, r"\2", x) for x in duplicates_external_ids]
    # Empty string indicates index 0
    duplicates_indices = ['0' if x=='' else x for x in duplicates_indices]
    # Indices as ints
    duplicates_indices = list(map(int, duplicates_indices))
    # Sort list
    duplicates_indices.sort()
    # Get max index
    new_sample_index = duplicates_indices[-1] + 1
    new_sample_external_id = "%s_%s"%(new_sample_base_external_id, new_sample_index)
    return new_sample_external_id

def get_current_samples():
    """Get current samples from FC
    """
    namespace = "nci-mimoun-bi-org"
    workspace = "CCLF_TSCA_2_0"
    wm = dalmatian.WorkspaceManager(namespace, workspace)
    # Current samples
    curr_samples = wm.get_samples()
    return curr_samples

##################################################################
### Read external sheet with metadata and save as .txt file
##################################################################
def read_external_sheet_and_save_as_txt_file(tsca_id, date, external_sheet_path):
    print("Reading external sheet with metadata and saving as .txt file")
    metadata_raw = pd.read_excel(external_sheet_path)
    metadata_raw.to_csv("/xchip/clf/seq_data/walkup_seq_sample_info/%s_walkupseq_%s_sample_info.txt"%(tsca_id.lower(), date), index=None, sep="\t")
    metadata_raw.to_csv("/xchip/clf/seq_data/process_for_fc/metadata_exports/walkupseq_files/%s_walkupseq_%s_sample_info.txt"%(tsca_id.lower(), date), index=None, sep="\t")
    metadata = pd.read_table("/xchip/clf/seq_data/walkup_seq_sample_info/%s_walkupseq_%s_sample_info.txt"%(tsca_id.lower(), date))
    return metadata

##################################################################
### Rename samples with duplicate external_ids
##################################################################
def rename_duplicate_external_ids(new_samples, target_dir_path, tsca_id, date, sn_id):
    """Rename duplicate external ids, i.e. if sample.external_id A already exists, rename it to A_2, A_3, and so on.
    Params:
        - new_samples: pd.DF with samples to be added to FC. Must have External ID column
    Returns: 
        - new_samples: pd.DF with 'external_id_validation_no_duplicates'
    """
    # Get current samples from FC
    curr_samples = get_current_samples()

    # Placeholder column
    new_samples["external_id_validation_no_duplicates"] = new_samples["External ID"]
    duplicate_sample_ids = []
    # Iterate over new samples, and find any samples in current set with the same external id
    for idx, sample in new_samples.iterrows():
        # # Reached end of df
        # if pd.isnull(sample["External ID"]):
        #     break
        duplicates = curr_samples[curr_samples["external_id_validation"].str.contains(sample["External ID"])]
        # Found existing sample with same external id
        if duplicates.shape[0] == 1 and duplicates["external_id_validation"].values[0] == sample["External ID"]:
            print("Found duplicate for %s, renaming to %s"%(sample["External ID"],  sample["External ID"] + "_2"))
            new_samples.loc[idx, "external_id_validation_no_duplicates"] = sample["External ID"] + "_2"
            duplicate_sample_ids.append(sample["External ID"] + "_2")
        # Found existing sample(s) with same base external id
        if duplicates.shape[0] > 1 and sample["External ID"] in duplicates["external_id_validation"].values:
            new_sample_external_id = get_new_sample_external_id(sample["External ID"], duplicates["external_id_validation"].values)
            print("Found duplicate for %s, renaming to %s"%(sample["External ID"], new_sample_external_id))
            new_samples.loc[idx, "external_id_validation_no_duplicates"] = new_sample_external_id
            duplicate_sample_ids.append(new_sample_external_id)

    duplicate_sample_ids_df = pd.DataFrame({"duplicate_sample_ids": duplicate_sample_ids})
    return new_samples, duplicate_sample_ids_df

##################################################################
### Copy all unidentified files to new directory where they will be renamed
##################################################################
def copy_unidentified_files_to_new_directory(sn_id, plate_prefix, target_dir_path):
    print("Creating directory for new batch...")
    res = os.system("mkdir -p %s"%target_dir_path)
    print(res)

    print("Copying all unidentified files to new directory where they will be renamed...")
    unidentified_files = [fn for fn in glob("/xchip/clf/seq_data/get.broadinstitute.org/pkgs/%s/%s_*"%(sn_id, plate_prefix)) if not os.path.basename(fn).endswith("fastq.gz")]
    for uf in unidentified_files:
        print("Copying %s"%uf)
        cmd = "cp -r %s %s/"%(uf, target_dir_path)
        res = os.system(cmd)
        print(res)

    return

#################################################################
# Rename files according to external id
################################################################
def rename_files_according_to_external_id(target_dir_path, plate_prefix, metadata):
    print("Starting renaming process...")
    # For each sample, find position name in plate and its corresponding external_id
    for idx, row in metadata.iterrows():
        position = row['Position']
        position = format_position(row['Position'])
        external_id = row['external_id_validation_no_duplicates'] # External ID
        external_id_fmt = external_id.replace("/", "")  # Format external_id: all upper case and no slashes

        print("Renaming files for sample in position [%s] with external_id [%s]"%(position, external_id_fmt))
        
        # Retrieve all the files in the current directory starting with @position
        path = "%s/%s_%s_*" % (target_dir_path, plate_prefix, position)
        sample_files = glob(path)
        if len(sample_files) == 0:
            print("No files found for this sample. Moving on to next sample...\n")
            continue
        
        # Create directory for this sample
        print("Creating directory [%s] for this sample..."%external_id_fmt)
        os.system("mkdir -p %s/%s"%(target_dir_path, external_id_fmt))
        
        # For each sample file (name starting with @position), rename position to corresponding external id
        print("Renaming all files with name starting with [%s] to name starting with [%s]"%(position, external_id))
        for filename in sample_files:
            new_name = re.sub(r"(%s_)([A-Z0-9]+)(_.*)"%plate_prefix, r"\1%s\3"%external_id_fmt, filename)
            # print("\t -> Renaming %s to %s"%(filename, new_name))
            os.rename(filename, new_name)
            # print("\t -> Moving %s to its sample directory..."%new_name)
            os.system("mv %s %s/%s"%(new_name, target_dir_path, external_id_fmt))
            # print("\n")
        print("---"*50, "\n\n")

    return

##################################################################
### Clean unnecessary files
##################################################################
def clean_unnecessary_files(plate_prefix, target_dir_path):
    print("Cleaning unnecessary files...")
    res = os.system("rm -r %s/%s_Solexa-*"%(target_dir_path, plate_prefix))
    res = os.system("rm -r %s/info_logs*"%(target_dir_path))
    res  = os.system("rm -r %s/*library*"%(target_dir_path))
    res  = os.system("rm -r %s/*html"%(target_dir_path))
    return

##################################################################
### Create sample info file
##################################################################
def create_sample_info_file(metadata, target_dir_path, duplicates):
    sample_info = pd.DataFrame()
    sample_info['individual_id'] = metadata['Collaborator Participant ID']
    sample_info['aggregation_product_name_validation'] = TSCA_version
    sample_info['external_id_validation'] = metadata['external_id_validation_no_duplicates'] # External ID
    sample_info['external_id_validation'] = sample_info['external_id_validation'].apply(lambda x: x.replace("/", ""))
    sample_info['bsp_sample_id_validation'] = metadata['Exported DNA SM-ID']
    sample_info['stock_sample_id_validation'] = metadata['Stock DNA SM-ID']
    sample_info['sample_type'] = metadata['Sample Type']
    sample_info['picard_aggregation_type_validation'] = picard_aggregation_type_validation
    sample_info['sample_id'] = sample_info['individual_id'] \
        .str.cat(sample_info['sample_type'], sep="-") \
        .str.cat(sample_info['bsp_sample_id_validation'], sep="-")
    sample_info['tumor_subtype'] = metadata['Tumor Type']
    sample_info['squid_sample_id_validation'] = sample_info['external_id_validation']
    sample_info['source_subtype_validation'] = metadata['Original Material Type']
    sample_info['processed_subtype_validation'] = metadata['Material Type']
    sample_info['primary_disease'] = metadata['Primary Disease']
    sample_info['media'] = metadata['Media on Tube']
    sample_info['Collection'] = metadata['Collection']
    sample_info['tissue_site'] = metadata['Tissue Site']

    ### Add BAM clean_bam_file_capture
    for idx, row in sample_info.iterrows():
        external_id = sample_info.loc[idx, 'external_id_validation']
        external_id_fmt = external_id.replace("/", "")
        bam_path_finds = glob("%s/%s/*.aligned.bam"%(target_dir_path, external_id_fmt))
        # bam_path_finds = glob("%s/%s/*.aligned.duplicates_marked.bam"%(target_dir_path, external_id_fmt))
        # Ensure bam path was found
        if len(bam_path_finds) == 0:
            print("No bam path found for sample %s"%external_id_fmt)
            bam_path = np.nan
        else:
            bam_path = bam_path_finds[0]
        sample_info.loc[idx, 'clean_bam_file_capture'] = bam_path


    #### Remove samples with no BAM file
    samples_no_bam = sample_info[pd.isnull(sample_info['clean_bam_file_capture'])]
    print("There are %s samples with no BAM file found"%(samples_no_bam.shape[0]))
    samples_no_bam.to_csv("%s/%s_%s_%s.missing_bam.txt"%(target_dir_path, tsca_id.lower(), date, sn_id), sep="\t", index=None)

    #### Save samples metadata
    sample_info_clean = sample_info[pd.notnull(sample_info['clean_bam_file_capture'])]
    sample_info_clean.to_csv("%s/%s_%s_%s.import_samples.txt"%(target_dir_path, tsca_id.lower(), date, sn_id), sep="\t", index=None)

    duplicates.to_csv("%s/%s_%s_%s.duplicate_sample_ids.txt"%(target_dir_path, tsca_id.lower(), date, sn_id), sep="\t", index=None)
    return


### IMPORTANT: 
# Run for each plate prefix (1_ and 2_), except the create_sample_info_file call
metadata = read_external_sheet_and_save_as_txt_file(tsca_id, date, external_sheet_path)
metadata, duplicates = rename_duplicate_external_ids(metadata, target_dir_path, tsca_id, date, sn_id)
print(metadata.columns)
copy_unidentified_files_to_new_directory(sn_id, plate_prefix, target_dir_path)
rename_files_according_to_external_id(target_dir_path, plate_prefix, metadata)
clean_unnecessary_files(plate_prefix, target_dir_path)
create_sample_info_file(metadata, target_dir_path, duplicates)
