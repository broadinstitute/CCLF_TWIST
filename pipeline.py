# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
# from IPython import get_ipython

# %%
# https://github.com/jkobject/JKBIO
# JKBio repo commit: 912087536d3cf6a7f1cbb00f9b131bc645780ee9 (9120875)
# from __future__ import print_function
import os.path
import dalmatian as dm
import pandas as pd
import sys

pathtoJK_parent = "../"
sys.path.append(pathtoJK_parent)
from JKBio import terra
import CCLF_processing as cclf
from IPython.core.debugger import set_trace
from src.helper import *

from JKBio.utils import *
import numpy as np
from gsheets import Sheets

# get_ipython().run_line_magic('load_ext', 'autoreload')
# get_ipython().run_line_magic('autoreload', '2')

"""
To generate the required files to read in from Google sheets:
Log into the Google Developers Console with the Google account whose spreadsheets you want to access.
Create (or select) a project and enable the Drive API and Sheets API (under Google Apps APIs).

https://console.developers.google.com/

Go to the Credentials for your project and create New credentials > OAuth client ID > of type Other.
In the list of your OAuth 2.0 client IDs click Download JSON for the Client ID you just created.
Save the file as client_secrets.json in your home directory (user directory).
Another file, named storage.json in this example, will be created after successful authorization
to cache OAuth data.

On you first usage of gsheets with this file (holding the client secrets),
your webbrowser will be opened, asking you to log in with your Google account to authorize
this client read access to all its Google Drive files and Google Sheets.
"""
sheets = Sheets.from_files('~/.client_secret.json', '~/.storage.json')
replace = {'T': 'Tumor', 'N': 'Normal', 'm': 'Unknown', 'L': 'Unknown'}

# %% [markdown]
# # CCLF TWIST Pipeline
#
# *go to the [readme](./README.md) to see more about execution*
#
#
# This pipeline has the following major steps:
# 1. Pull in information about the TWIST batch(es) from Google sheet(s).
# 2. Create a TSV of the new sample information
# 3. Create a TSV of the new sample set information (e.g. cohorts)
# 4. Upload the sample information and sample set TSVs to the Terra workspace
# 5. Run Terra workflows to get copy number (CNV) and mutation (SNV) information, and to create copy number heat maps by batch and by cohort.
#
# %% [markdown]
# # Initialization
# Pull in information about the TWIST batch(es) from Google sheet(s).
#
# **Note:** Each time, the `samplesetnames` and the `gsheeturllist` need to be updated.

# %%
# create sample set names for each batch in *chronological* order (e.g. CCLF_TWIST1 before CCLF_TWIST2)
# if you only have one batch to run, still make it a list e.g. ["CCLF_TWIST1"]
# this ensures that the pipeline will run as designed
samplesetnames = ['CCLF_TWIST36']


# list of the external sheets produced for each batch you want to run through the pipeline
# TO EDIT:
gsheeturllist = ['https://docs.google.com/spreadsheets/d/1X0VvxsefdWdSDWPNYDd0V76ZOJpdyJQ43B78EciRtvs/edit#gid=0']

# generate the sample set names we will use in Terra
samplesetnames_normals = [s + '_normals' for s in samplesetnames]
samplesetnames_tumors = [s + '_tumors' for s in samplesetnames]
samplesetnames_pairs = [s + '_pairs' for s in samplesetnames]
samplesetnames_all = [s + '_all' for s in samplesetnames]

# workspace where we are pulling in the data from
data_workspace="terra-broad-cancer-prod/Cancer_Cell_Line_Factory_CCLF_PanCancer_PanelSeq"
# workspace where we are running the workflows
proc_workspace="nci-mimoun-bi-org/PANCAN_TWIST copy"

# TODO: these are hard-coded in the helper.py file. However, I don't know how these columns are being used (if at all) and I think these may be meaningless.
# source="CCLF"
# picard_aggregation_type_validation="PCR"

# mapping abbreviations to full names/descriptions
cohorts2id="https://docs.google.com/spreadsheets/d/1R97pgzoX0YClGDr5nmQYQwimnKXxDBGnGzg7YPlhZJU"

# %% [markdown]
# # Connection errors? Reload the chunk below
# For example,
# `ConnectionError: ('Connection aborted.', OSError("(54, 'ECONNRESET')"))`

# %%
wfrom = dm.WorkspaceManager(data_workspace)
wto = dm.WorkspaceManager(proc_workspace)

# %% [markdown]
# # Getting the samples
#
# - we load the samples from data workspace and load the metadata files
# - we remove data that has already been processed
# - we create the final ids

# %%
# create preliminary versions of the sample and metadata tables
newsamples, newmetadata = create_preliminary_sample_and_metadata_tables(wto, wfrom, samplesetnames, external_sheets_url_list=gsheeturllist, cohorts2id_url=cohorts2id)

# %% [markdown]
# # Creating the sample information dataframe
# Create a dataframe of the new sample information
#
# **Note:** It can be difficult to recreate the sample_info variable below after you have already uploaded TSVs to Terra since this pipeline specifically looks for samples that do not already exist in the workspace. When running the pipeline on a new batch of data, **I recommend writing the final sample_info to a file.**
#
# **Note 2:** We replace all "/" in the External IDs with "_". This prevents errors when filepaths are created using the external IDs in Terra.
# %% [markdown]
# ## Required metadata columns
# We do not include samples that are missing information in any of the following columns in the external sheet:
# - Collaborator Participant ID
# - Exported DNA SM-ID
# - Stock DNA SM-ID
# - Participant ID <- This is the patient identifier
# - Sample Type
# - Original Material Type
# - Material Type
# - Primary Disease <- Only the technical controls won't have this information.
# - Collection
# - Tissue Site
#
# Without this list of metadata, the samples will not be added to Terra.

# %%
# merge the data from the External Sheet(s) and the data from the data source (e.g. Broad genomics delivery)
df = pd.concat([newmetadata, newsamples], axis=1, sort=True)

# specify required metadata columns
tolook = ['Collaborator Participant ID','Exported DNA SM-ID', 'Stock DNA SM-ID', 'Participant ID', 'Sample Type','Tissue Site', 'Original Material Type', 'Material Type','Primary Disease', 'Collection']


# %%
# If any samples are missing some of the required metadata, stop now and ask the CCLF team to fill out the missing values in the External Sheet.
check_required_metadata_columns(df, tolook, cohorts2id_url=cohorts2id, drop=False)


# %%
# only keep samples that have all the appropriate metadata information
df = check_required_metadata_columns(df, required_metadata_cols=tolook, cohorts2id_url=cohorts2id, drop=True)


# %%
# generate sample df to upload to Terra
sample_info = create_sample_df_for_terra(wto, df, cohorts2id_url=cohorts2id)

# sanity check: this should be what you plan on uploading to Terra
print(sample_info.shape)
print(sample_info.head())


# %%
# Run this chunk to save the sample_info TSV to a file. I highly recommend this when running a pipeline on a new batch.
# This way, if anything goes wrong in the workspace, you can fall back to this.

## check: create directory "data/sample_infos" if does not exist
filepath = 'data/sample_infos/%s_sample_info.tsv' % '_'.join(samplesetnames)
sample_info.to_csv(filepath, sep='\t')


# %%
# read in the file you just saved
filepath = 'data/sample_infos/%s_sample_info.tsv' % '_'.join(samplesetnames)
sample_info = pd.read_csv(filepath, sep = '\t', na_filter = False)
sample_info = sample_info.set_index('sample_id')
print(sample_info.shape)
sample_info.head()

# %% [markdown]
# # Creating the pairs
# Create a TSV of the new pairs information.

# %%
newpairs = create_pairs_table(sample_info, wto)

# %% [markdown]
# # Uploading samples and pairs to Terra

# %%
print("uploading new samples...")
wto.upload_samples(sample_info)
if not "NA" in wto.get_samples().index.tolist():
    wto.upload_samples(pd.DataFrame({'sample_id':['NA'], 'participant_id':['NA']}).set_index('sample_id'))
    wto.upload_samples(sample_info)


# %%
print("uploading pairs...")
wto.upload_entities('pair', newpairs)

# %% [markdown]
# # Create pair sets and sample sets
#
# In the following cell, we create:
# - a pair set for each batch
# - sample sets for each batch
# - sample sets for each cohort
#
# And then we upload these entities to the Terra workspace.

# %%
# Create pairs per batch dictionary
dict_pairs_per_batch = create_dict_of_pairs_per_sampleset(newpairs, sample_info, samplesetnames, save=True)

# Load from saved file
dict_pairs_per_batch = np.load('dict_pairs_per_batch.npy',allow_pickle='TRUE').item()


# %%
#  functionalized version
create_sample_sets_per_batch(sample_info, samplesetnames, samplesetnames_all, samplesetnames_tumors, samplesetnames_normals, proc_workspace)

create_pair_sets_per_batch(samplesetnames, samplesetnames_pairs, proc_workspace, dict_pairs_per_batch)

create_samplesets_and_pairsets_per_cohort(sample_info, samplesetnames, proc_workspace, cohorts2id, newpairs)

# get a list of all normals in the processing workspace
# TODO: I think by this point in the pipeline all of the samples will have been uploaded to Terra. Thus, I think we may only need the following one-liner to get a list of all normals:
# all_normals = [i for i, r in wto.get_samples().iterrows() if r['sample_type'] == "Normal"]

normalsid = [i for i, r in sample_info.iterrows() if r['sample_type'] == "Normal"]
print(len(normalsid))

refsamples = wto.get_samples()
normalsid.extend([k for k, _ in refsamples.iterrows() if _.sample_type == "Normal"])
print(len(normalsid))

create_aggregate_samplesets(normalsid, proc_workspace, wto)

# %% [markdown]
# ***
# ***
# # Push to Git repo (CCLF_TWIST) now!
# This way, other people can easily take over the process of running the pipelines and feel confident that they have the most up-to-date version of this Jupyter notebook.
# ***
# ***
# %% [markdown]
# # Running Terra Worlflows
# Run Terra workflows to get copy number (CNV) and mutation (SNV) information, and to create copy number heat maps by batch and by cohort.
#
# The order of running the workflows is as follows:
# - RenameBAM_TWIST
# - CalculateTargetCoverage_PANCAN,
#     + DepthOfCov_PANCAN
# - CreatePanelOfNormalsGATK_PANCAN, (edit the output config "normals_pon attribute"))
#     + DepthOfCovQC_PANCAN
# - CallSomaticCNV_PANCAN (edit the input config to match the output from CreatePanelOfNormalsGATK_PANCAN)
# - PlotSomaticCNVMaps_PANCAN: we plot CN heat maps for each batch and also for each cohort
# - MutationCalling_Normals_TWIST
# - FilterGermlineVariants_NormalSample_TWIST
# (edit the "PoN_name" config for CreatePoNSNV_Mutect1 and CreatePoNSNV_Mutect2)
# - CreatePoNSNV_Mutect1,
#     + CreatePoNSNV_Mutect2
# - SNV_PostProcessing_Normals,
#     + MutationCalling_Tumors_TWIST (edit the input config to match pon_mutect1, pon_mutect2)
# - FilterGermlineEvents_TumorSample
# - SNVPostProcessing_TWIST,
#     + FNG_Compile_Pileup_Cnt
# - FNG_Compile_db_slow_download
# - **do manual step here on local machine: need to merge two fingerprinting tables**
# - FNG_Query_db
#
# More information about the pipeline exist here: https://cclf.gitbook.io/tsca/
#
# **Note 1:** If for som reason, one of the terra submission function gives no output and it does not seem to submit anything to terra, it might be that you have been logged out of terra you will have to reload the workspace manager and package.
#
# **Note 2:** If you get the preflight error "expression and etype must BOTH be None or a string value", check the workflow configuration using wto.get_config("NAME_OF_WORKFLOW"). This error usually occurs when you pass in expression and etype information, but the etype is already set as the "rootEntity" aka the default for the workflow. You can fix this by either changing the workflow configuration in Terra, or by not passing in the etype or expression. If you want to see why this error occurs, look at the preflight function in lapdog.py (https://github.com/broadinstitute/lapdog/blob/master/lapdog/lapdog.py).

# %%
print("Creating Terra submissions: remember you can only cancel \n or interact with terra submissions from the Terra website. \n https://app.terra.bio/#workspaces/"+proc_workspace.replace(" ", "%20")+"/job_history")

RenameBAM_TWIST = terra.createManySubmissions(proc_workspace, "RenameBAM_TWIST", samplesetnames_all,
                                              entity='sample_set', expression='this.samples')

print("waiting for 'Rename'")
terra.waitForSubmission(proc_workspace, RenameBAM_TWIST)


# %%
CalculateTargetCoverage_PANCAN = terra.createManySubmissions(proc_workspace, "CalculateTargetCoverage_PANCAN", samplesetnames_all,
                                              entity='sample_set', expression='this.samples')
DepthOfCov_PANCAN = terra.createManySubmissions(proc_workspace, "DepthOfCov_PANCAN", samplesetnames_all,
                                              entity='sample_set', expression='this.samples')

print("waiting for 'CalculateTargetCoverage' & 'DepthOfCov_PANCAN'")
combined_list = CalculateTargetCoverage_PANCAN + DepthOfCov_PANCAN
terra.waitForSubmission(proc_workspace, combined_list)


# %%
## Updates the config for each batch id
CreatePanelOfNormalsGATK_PANCAN = []
DepthOfCovQC_PANCAN = []
for ind, batch_id in enumerate(samplesetnames):
    # get current config for workflow that creates the PON for CNV calling
    createPON_config = wto.get_config('CreatePanelOfNormalsGATK_PANCAN')
    # edit the config
    createPON_config['outputs']['CreatePanelOfNormals.combined_normals'] = 'workspace.combined_normals_' + batch_id
    createPON_config['outputs']['CreatePanelOfNormals.normals_pon'] = 'workspace.pon_normals_' + batch_id
    createPON_config['outputs']
    # update the config in Terra
    wto.update_config(createPON_config)

    # create batch-specific PON to be used for CNVs
    CreatePanelOfNormalsGATK_PANCAN.append(wto.create_submission("CreatePanelOfNormalsGATK_PANCAN", samplesetnames_normals[ind]))
    DepthOfCovQC_PANCAN.append(wto.create_submission("DepthOfCovQC_PANCAN", samplesetnames_all[ind], etype='sample_set', expression='this.samples'))


# %%
print("waiting for 'DepthOfCovQC_PANCAN' & 'CNV_CreatePoNForCNV'")
combined_list = DepthOfCovQC_PANCAN + CreatePanelOfNormalsGATK_PANCAN
terra.waitForSubmission(proc_workspace, combined_list)


# %%
CallSomaticCNV_PANCAN = []
for ind, batch_id in enumerate(samplesetnames):
    # get current config
    CNV_config = wto.get_config('CallSomaticCNV_PANCAN')
    CNV_config['inputs']['CallSomaticCNV.normals_pon']

    # edit the config
    CNV_config['inputs']['CallSomaticCNV.normals_pon'] = 'workspace.pon_normals_' + batch_id
    CNV_config['inputs']

    # update the config in Terra
    wto.update_config(CNV_config)
    CallSomaticCNV_PANCAN.append(wto.create_submission("CallSomaticCNV_PANCAN", samplesetnames_all[ind], etype='sample_set', expression='this.samples', use_callcache = True))


# %%
print("waiting for 'CallSomaticCNV_PANCAN'")
terra.waitForSubmission(proc_workspace, CallSomaticCNV_PANCAN)


# %%
# in case you want to rerun the notebook and generate the sample key variables
sample_info, all_pairsets, cohorts_per_batch, cohort_pairsets, all_changed_cohorts =     regenerate_variables(wto, samplesetnames_all, cohorts2id)


# %%
all_changed_cohorts


# %%
# create CNV map for each batch
terra.createManySubmissions(proc_workspace, "PlotSomaticCNVMaps_PANCAN", samplesetnames_all, use_callcache = False)
# create CNV map for each cohort
terra.createManySubmissions(proc_workspace, "PlotSomaticCNVMaps_PANCAN", list(all_changed_cohorts), use_callcache = False)

print("submitted final jobs for CNV pipeline")
print("you don't need to wait before moving onto the next cell")


# %%
MutationCalling_Normals_TWIST = terra.createManySubmissions(proc_workspace, "MutationCalling_Normals_TWIST", samplesetnames_normals,
                                              entity='sample_set', expression='this.samples')
print("waiting for 'MutationCalling_Normals_TWIST'")
terra.waitForSubmission(proc_workspace, MutationCalling_Normals_TWIST)


# %%
FilterGermlineVariants_NormalSample_TWIST = terra.createManySubmissions(proc_workspace, "FilterGermlineVariants_NormalSample_TWIST", samplesetnames_normals,
                                              entity='sample_set', expression='this.samples', use_callcache=True)
print("waiting for 'SNV_FilterGermline'")
terra.waitForSubmission(proc_workspace, FilterGermlineVariants_NormalSample_TWIST)


# %%
# get current config
mutect1_config = wto.get_config('CreatePoNSNV_Mutect1')
mutect2_config = wto.get_config('CreatePoN_SNV_MuTect2')

# edit the config
mutect1_config['inputs']['CreatePanelOfNormals.PoN_name'] = '"Cum_PoN_' + samplesetnames[-1] + '_all_vcf_mutect1"'
mutect2_config['inputs']['CreatePanelOfNormals.PoN_name'] = '"Cum_PoN_' + samplesetnames[-1] + '_all_vcf_mutect2"'
mutect1_config['outputs']['CreatePanelOfNormals.normals_pon_vcf'] = 'workspace.Cum_PoN_' + samplesetnames[-1] + '_all_vcf_mutect1'
mutect2_config['outputs']['CreatePanelOfNormals.createPanelOfNormals.normals_pon_vcf'] = 'workspace.Cum_PoN_' + samplesetnames[-1] + '_all_vcf_mutect2'

# update the config in Terra
wto.update_config(mutect1_config)
wto.update_config(mutect2_config)

# create PON for SNV from all the normals we have in the workspace so far
CreatePoNSNV_Mutect1 = wto.create_submission('CreatePoNSNV_Mutect1', "All_normals_TWIST")
CreatePoN_SNV_MuTect2 = wto.create_submission('CreatePoN_SNV_MuTect2', "All_normals_TWIST")


# %%
print("waiting for 'CreatePoN_SNV_MuTect2' & 'CreatePoNSNV_Mutect1'")
terra.waitForSubmission(proc_workspace, [CreatePoNSNV_Mutect1, CreatePoN_SNV_MuTect2])

# %% [markdown]
# ## Note: It may be okay if some samples fail the MutationCalling_Tumors_TWIST workflow. Samples will fail if no mutations made it through Mutect1 and Mutect2's filters.
# The MutationCalling_Tumors_TWIST pipeline has been updated to use GATK4, and there are many more pre-filters for Mutect2 that greatly reduce the computation time required. As part of this change, however, we discovered that the next step (FilterMutectCalls) will fail if the vcf it gets from Mutect2 is empty. This can happen if all the variants are filtered out. Thus, long story short, if the sample fails at the FilterMutectCalls step and the log file shows that there were no variants left after Mutect1 and Mutect2, then this failure is not something to worry about.
#
# The details can be found at https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2 by searching for "Read filters". In addition to the prefilters described in that section, Mutect2 also prefilters sites that are in the matched normal and the PoN.

# %%
SNV_PostProcessing_Normals = []
MutationCalling_Tumors_TWIST = []
for ind, batch_id in enumerate(samplesetnames):

    # get config
    mutcall_tumor = wto.get_config('MutationCalling_Tumors_TWIST')

    # edit the config
    mutcall_tumor['inputs']['MutationCalling_Tumor.pon_mutect1'] = 'workspace.Cum_PoN_' + samplesetnames[-1] + '_all_vcf_mutect1'
    mutcall_tumor['inputs']['MutationCalling_Tumor.pon_mutect2'] = 'workspace.Cum_PoN_' + samplesetnames[-1] + '_all_vcf_mutect2'

    # check config
    print(mutcall_tumor['inputs']['MutationCalling_Tumor.pon_mutect1'])
    print(mutcall_tumor['inputs']['MutationCalling_Tumor.pon_mutect2'])

    # update the config in Terra
    wto.update_config(mutcall_tumor)

    # create submission
    SNV_PostProcessing_Normals += [wto.create_submission("SNV_PostProcessing_Normals", samplesetnames_normals[ind])]

    MutationCalling_Tumors_TWIST += [wto.create_submission("MutationCalling_Tumors_TWIST", samplesetnames_pairs[ind], etype='pair_set', expression='this.pairs')]


# %%
print("waiting for 'SNV_PostProcessing_Normals' & 'MutationCalling_Tumors_TWIST'")
combined_list = SNV_PostProcessing_Normals + MutationCalling_Tumors_TWIST
terra.waitForSubmission(proc_workspace, combined_list)


# %%
## note: you might see that some of the cohorts fail on this workflow. That can be expected: the workflow needs cohorts with at least 2 acceptable cell lines to run (if only 1, then the workflow will fail)
FilterGermlineEvents_TumorSample = terra.createManySubmissions(proc_workspace, 'FilterGermlineEvents_TumorSample', samplesetnames_pairs, 'pair_set', expression='this.pairs')
print("waiting for 'FilterGermlineEvents_TumorSample'")
terra.waitForSubmission(proc_workspace, FilterGermlineEvents_TumorSample)


# %%
# create aggregate SNV tsvs for each batch
terra.createManySubmissions(proc_workspace, "SNVPostProcessing_TWIST", samplesetnames_pairs)
# create aggregate SNV tsvs for each cohort
terra.createManySubmissions(proc_workspace, "SNVPostProcessing_TWIST", list(cohort_pairsets))
print("Submitted final jobs for SNV pipeline")

# %% [markdown]
# ## Fingerprinting (FNG)

# %%
# determine the FNG pileup counts for each sample
FNG_Compile_Pileup_Cnt = terra.createManySubmissions(proc_workspace, "FNG_Compile_Pileup_Cnt", samplesetnames_all, entity='sample_set', expression='this.samples')
print("waiting for 'FNG_Compile_Pileup_Cnt'")
terra.waitForSubmission(proc_workspace, FNG_Compile_Pileup_Cnt)


# %%
# create the proper "super" set of samples to run through the FNG compiler
# if processing multiple batches at once, this will correspond to a sample "super" set containing the sample IDs from all the new batches. If processing a single batch, this will just be that batch's sample set.
samples_to_add = []
for set_name in samplesetnames_all:
#     samples_to_add += wto.get_sample_attributes_in_set(set_name).index.tolist()
    samples_to_add += wto.get_sample_sets().loc[set_name, 'samples']
samples_to_add
print(samples_to_add)

fng_sampleset_id = "_".join(samplesetnames)
print(fng_sampleset_id)

terra.addToSampleSet(workspace = proc_workspace, samplesetid = fng_sampleset_id, samples = samples_to_add)


# %%
## Update the output config to create the new FNG database
# get current config for the FNG compiling workflow
fngCompile_config = wto.get_config('FNG_Compile_db_slow_download')
# edit the config
# TODO: would be nice to be able to change the name of the outputted FNG database. Right now, all named the same.
# fngCompile_config['inputs']['FNG_Compile_db.compile_db.output_file_name'] = '"fingerprinting_db_through_' + samplesetnames[-1] + '.txt"'
fngCompile_config['outputs']['FNG_Compile_db.compile_db.fingerprinting_db'] = 'workspace.fingerprinting_db_through_' + samplesetnames[-1]
fngCompile_config['outputs']['FNG_Compile_db.compile_db.fingerprinting_db_current'] = 'workspace.fingerprinting_db'
fngCompile_config['outputs']

print(fngCompile_config)
# update the config in Terra
wto.update_config(fngCompile_config)


# %%
# create FNG db using Method Version 7
# TODO: update method to change the output file name to something more descriptive
# pass in a sample set containing all of the new samples you're processing
# do not use call cache; we need to see if the github repo has been updated and thus must clone each time
FNG_Compile_db_slow_download = wto.create_submission("FNG_Compile_db_slow_download", fng_sampleset_id, use_callcache=False)
print("waiting for 'FNG_Compile_db'")
terra.waitForSubmission(proc_workspace, FNG_Compile_db_slow_download)

# %% [markdown]
# ### "FNG_Compile_db_slow_download" command is problematic currently
# **NOTE**: The "FNG_Compile_db_slow_download" command is problematic currently in this workspace because this workspace only contains TWIST samples, but we want to be able to look at fingerprinting data from both TSCA and TWIST. Currently, Gwen merges the previous fingerprinting_db.txt file with the newly created fingerprinting_db.txt. We'll have to repeat this merging procedure unless we edit the workflow. Gwen has started this process, but hasn't finished the edits (just need to build the proper docker container). So for now (1/15/20), still need to do the merging locally (unfortunately).
#
# See the R file: "src/FNG_TWIST_and_TSCA_merge" in the CCLF_TWIST GitHub repo.
# %% [markdown]
# ### TODO: this is where the merging work in R needs to be performed on local.
# See the R file: `CCLF_TWIST/src/FNG_TWIST_and_TSCA_merge.Rmd` in the CCLF_TWIST GitHub repo.
# -> The details of what you need to do are in this Rmd file.
# -> Once the local work is done, continue running through the chunks below.
#
# I also have the same merging script in `CCLF_TWIST/workflow/scripts/merge_fingerprinting_dfs.R`, and a WDL script to run this at `CCLF_TWIST/workflow/scripts/merge_fng_databases.wdl`. This has been uploaded to Terra, too, and I have made a workflow for this in the data processing workspace: "merge_FNGs". However, the Docker image is missing the Argparse package. I think once this is added to the Docker image, the workflow should run smoothly on Terra.
#

# %%
# 8/4/2021, Javad's efforts to implement the R code in the python cell
# This seems to be redundant though. For some reason, the binding of the
# data already has happened in the FNG_Compile_db_slow_download output
# This seems to at least be the case for TWIST35

import pandas as pd

workspace_attr = wto.get_workspace_metadata()['workspace']['attributes']
fp_old = pd.read_csv(workspace_attr['fingerprinting_db_through_CCLF_TWIST34'], sep='\t')
fp_new = pd.read_csv(workspace_attr['fingerprinting_db_through_' + samplesetnames[-1]], sep='\t')
assert not fp_new.duplicated().any()
assert not fp_old.duplicated().any()
assert set(fp_new.columns) == set(fp_old.columns)
assert (set(fp_old['batch']) - set(fp_new['batch'])) == set([])

# not sure why this is necessary. It produces identical results
fp_all = pd.concat([fp_new, fp_old])
fp_all.drop_duplicates(inplace=True)
fp_all.reset_index(drop=True, inplace=True)
assert fp_new.shape == fp_all.shape
assert fp_all.equals(fp_new)
fp_all_filename = '/tmp/fingerprinting_db_through_{}.txt'.format(samplesetnames[-1])
print('It seems like the merging step of the pipeline does not do anything and can be dropped.\nUntil this refactoring is implemented will save the Python-generated merged file to {}'.format(fp_all_filename))
fp_all.to_csv(fp_all_filename, sep='\t', index=False)



# %%
# for each batch, query the FNG database
FNG_Query_db = terra.createManySubmissions(proc_workspace, "FNG_Query_db", samplesetnames_all)
print("Submitted final FNG Job")
terra.waitForSubmission(proc_workspace, FNG_Query_db)

print('\n\n\n')
print('#'*40)
print('Done')
print('#'*40)

# %% [markdown]
# # You've finished running through the pipeline!
# You should have all the SNV, CNV, and FNG results ready in Terra. Update the Asana task.
# %% [markdown]
# <!-- # If got a new cohort label / abbreviation and need to update data that already exists in Terra: -->

# %%
# # Get data from Gsheet metadata
# metadata = pd.concat(gsheets,sort=False, keys = samplesetnames)
# metadata = metadata.reset_index().rename(columns = {'level_0':'batch', "External ID":'external_id_validation'}).drop(['level_1'], axis = 'columns')
# metadata.index = metadata['Exported DNA SM-ID']

# display(metadata.head())

# # Pull relevant sample_info from Terra
# sample_info = wto.get_samples()
# sample_info = sample_info[sample_info["batch"].isin(samplesetnames)]
# display(sample_info.head())
# # display(sample_info.loc[:,["cohorts", "Collection"]])

# # Merge new Metadata with stuff existing in Terra (in particular, we often want to update the Collections and cohorts columns)
# updated = pd.concat([sample_info.drop(columns=['Collection']), metadata["Collection"].reindex(sample_info.index)], axis=1)
# updated.head()
# updated.columns.tolist()
# updated = updated.reindex(columns=(['Collection', 'cohorts'] + list([a for a in updated.columns if a not in ['Collection', 'cohorts']]) ))
# updated

# updated = getCohortAbbreviations(updated)
# updated
# print("Final 'updated' df, looking at just the cohorts and Collection columns:")
# display(updated.loc[:,["cohorts", "Collection"]].head())
