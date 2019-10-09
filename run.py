# By Jérémie Kalfon @jkobject
# For the BroadInsitute of MIT and Harvard
# jkobject@gmail.com
# 07/2019

from __future__ import print_function
import os.path
import dalmatian as dm
import pandas as pd
import sys
import TerraFunction as terra
from Helper import *
import numpy as np
from gsheets import Sheets

# https://github.com/jkobject/JKBIO

"""
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


def update(samplesetname,
           date="2019",
           data_namespace="broad-genomics-delivery",
           data_workspace="Cancer_Cell_Line_Factory_CCLF_PanCancer_PanelSeq",
           proc_namespace="nci-mimoun-bi-org",
           proc_workspace="PANCAN_TWIS_Dev",
           source="CCLF",
           site="HT33MBCX2",
           tsca_id="TSCA45",
           TSCA_version="TSCA Rapid Cancer Detection Panel v2",
           picard_aggregation_type_validation="PCR",
           forcekeep=[],
           cohorts2id="https://docs.google.com/spreadsheets/d/1R97pgzoX0YClGDr5nmQYQwimnKXxDBGnGzg7YPlhZJU/edit#gid=872582930",
           gsheeturllist=["https://docs.google.com/spreadsheets/d/1LR8OFylVClxf0kmZpAdlVjrn3RBcfZKpNoDYtKdnHB8",
                          "https://docs.google.com/spreadsheets/d/128dkFhL1A0GqTjmR7iMvBZE8j6ymO8krBL9WX-wUAn4",
                          "https://docs.google.com/spreadsheets/d/1S3DqBdVkd9dLP1PDYcdSWuD2Iy2gJpzuYBhvmP37UxU",
                          "https://docs.google.com/spreadsheets/d/1kVIeIw66AxWLhAZlqUnAY17S87Rtfhijf1o3x0hG3Jw"]):
  """
  get the non overlapping samples from a data workspace to a processing workspace

  Adds them in a manner consistent to the CCLF processing model for analysis and fingerprinting,
  creating the necessary pairs and sets.
  It will output nothing but will have updated the different Terra tsvs so that it contains:
  - the new samples
  - the new participants
  - pairs for each samples tumor_normals
  - sets for all samples, current batch and all_normals
  - pair sets for the batch and the different cohorts.

  args:
  - date: (opts) str if one wants to add a processing date to the current batch
  - samplesetname: str the name of the sampleset, i.e. batch
  - data_namespace: str
  - data_workspace: str
  - proc_namespace: str
  - proc_workspace: str
  - source: str
  - site: str (opts)
  - tsca_id: str (opts)
  - TSCA_version: str (opts)
  - picard_aggregation_type_validation: (opts)
  - forcekeep: (opts) list[str] different samples you would want to reupload
  - cohorts2id: (opts) str path to a googlesheet containing the match : cohorts_name / cohort id
  - gsheeturllist: (opts) list[str] url to google sheets where metadata for samples might be

  """
  wfrom = dm.WorkspaceManager(data_namespace, data_workspace)
  wto = dm.WorkspaceManager(proc_namespace, proc_workspace)
  # we look at all the samples we already have
  # we look at all the samples we already have
  refsamples = wto.get_samples()
  refids = refsamples.index
  cohorts = sheets.get(cohorts2id).sheets[0].to_frame()
  # we use this gsheet package to get all the sheets into one dataframe
  metadata = pd.concat([sheets.get(url).sheets[0].to_frame() for url in gsheeturllist])

  # we do some corrections just in case
  samples1 = wfrom.get_samples().replace(np.nan, '', regex=True)

  # creating sample_id (like in processing workspace) for metadata and samples1
  metadata = metadata.dropna(0, subset=['Collaborator Sample ID'])
  ttype = [i for i in metadata["Sample Type"]]
  metadata['sample_id'] = [val['Collaborator Participant ID'] + '-' + val['Sample Type'] + '-' + val['Exported DNA SM-ID'] for i, val in metadata.iterrows()]

  sample_id = [val["individual_alias"] + '-' + val['sample_type'] + '-' + i.split('_')[2] for i, val in samples1.iterrows()]
  samples1.index = sample_id

  # filtering on what already exists in the processing workspace (refids)
  newsamples = samples1[(~samples1.index.isin(refids)) | samples1.index.isin(forcekeep)]
  tokeep = set(metadata.index) & set(newsamples.index)

  # usefull to merge the two df, sm-id is one of the only unique id here
  if len(newsamples[~newsamples.index.isin(tokeep)]) > 0:
    print('we could not add these as we dont have metadata for them:' + str(newsamples[~newsamples.index.isin(tokeep)]))
  newsamples = newsamples[newsamples.index.isin(tokeep)]
  newmetadata = metadata[metadata.index.isin(tokeep)]

  print('creating new df')
  df = pd.concat([newmetadata, newsamples], axis=1, sort=True)
  # from this new set we create a dataframe which will get uploaded to terra
  sample_info = df[['crai_or_bai_path', 'cram_or_bam_path']]
  sample_info['individual_id'] = df['Collaborator Participant ID']
  sample_info['reference_id'] = df['Exported DNA SM-ID']
  sample_info['participant'] = df['Collaborator Participant ID']
  sample_info['aggregation_product_name_validation'] = [TSCA_version] * sample_info.shape[0]
  # here we add this number as the reference id might be present many times already for different samples
  # in the processing workspace
  sample_info['external_id_validation'] = [i + '_' + str(refsamples[refsamples['external_id_validation'] == i].shape[1]) if refsamples[refsamples['external_id_validation'] == i].shape[0] > 0 else i for i in sample_info['reference_id']]
  sample_info['bsp_sample_id_validation'] = df.index
  sample_info['stock_sample_id_validation'] = df['Stock DNA SM-ID']
  sample_info['sample_type'] = df['Sample Type']
  sample_info['picard_aggregation_type_validation'] = [picard_aggregation_type_validation] * sample_info.shape[0]
  sample_info['tumor_subtype'] = df['Tumor Type']
  sample_info['squid_sample_id_validation'] = sample_info['external_id_validation']
  sample_info['source_subtype_validation'] = df['Original Material Type']
  sample_info['processed_subtype_validation'] = df['Material Type']
  sample_info['primary_disease'] = df['Primary Disease']
  sample_info['media'] = df['Media on Tube']
  sample_info['Collection'] = df['Collection']
  # match collection data and error out
  cohortlist = []
  for k, val in sample_info['Collection'].iteritems():
    res = cohorts[cohorts['Name'] == val]
    if len(res) == 0:
      raise "we do not have a correponsding cohort for this collection"
    cohortlist.append(res['ID'].values[0])
  sample_info['cohorts'] = cohortlist

  sample_info['tissue_site'] = df['Tissue Site']
  sample_info['source'] = [source] * sample_info.shape[0]
  sample_info['sample_id'] = df.index

  sample_info = sample_info.set_index('sample_id')

  # creating the sample_sets
  normals = [r["participant"] for i, r in sample_info.iterrows() if r['sample_type'] == "Normal"]
  normalsid = [i for i, r in sample_info.iterrows() if r['sample_type'] == "Normal"]
  tumors = [r["participant"] for i, r in sample_info.iterrows() if r['sample_type'] == "Tumor"]
  tumorsid = [i for i, r in sample_info.iterrows() if r['sample_type'] == "Tumor"]
  prevtumors = [val["participant"] for k, val in refsamples.iterrows() if val.sample_type == "Tumor"]
  prevnormals = [val["participant"] for k, val in refsamples.iterrows() if val.sample_type == "Normal"]

  print("creating new pairs")
  # do we have new tumors/normals for our previous ones
  newpairs = {'pair_id': [], 'case_sample': [], 'control_sample': [], 'participant': []}

  toreprocess_normals = set(tumors) & set(prevnormals)
  for val in toreprocess_normals:
    for tumor_id in sample_info[sample_info['participant'] == val][sample_info[
            'sample_type'] == 'Tumor'].index.tolist():
      normal_id = refsamples[refsamples['participant'] == val][refsamples[
          'sample_type'] == 'Normal'].index.tolist()[0]
      newpairs['pair_id'].append(tumor_id + '_' + normal_id)
      newpairs['case_sample'].append(tumor_id)
      newpairs['control_sample'].append(normal_id)
      newpairs['participant'].append(val)

  paired = set(tumors) & set(normals)
  for val in set(tumors) - toreprocess_normals:
    for tumor_id in sample_info[sample_info['participant'] == val][sample_info[
            'sample_type'] == 'Tumor'].index.tolist():
      normal_id = sample_info[(sample_info['participant'] == val) & (sample_info[
          'sample_type'] == 'Normal')].index.tolist()[0] if val in paired else 'NA'
      newpairs['pair_id'].append(tumor_id + "_" + normal_id)
      newpairs['case_sample'].append(tumor_id)
      newpairs['control_sample'].append(normal_id)
      newpairs['participant'].append(val)

  newpairs = pd.DataFrame(newpairs).set_index('pair_id')

  print("all the entities need to exist! Else it will raise an error and block further uploads to Terra")
  print("uploading new samples")
  wto.upload_samples(sample_info)
  if not "NA" in wto.get_samples().index.tolist():
    wto.upload_samples(pd.DataFrame({'sample_id': ['NA'], 'participant_id': ['NA']}).set_index('sample_id'))

  print("creating pairs and pairsets")
  wto.upload_entities('pair', newpairs)
  wto.update_pair_set(samplesetname + '_pairs', newpairs.index)
  cohorts_in_batch = []
  cohorts_with_pairs = []
  for val in cohorts['ID'].values:
    cohortsamples = sample_info[sample_info["cohorts"] == val].index.tolist()
    tumorsamplesincohort = sample_info[sample_info["cohorts"] == val][sample_info['sample_type'] == "Tumor"].index.tolist()
    pairsamples = newpairs[newpairs['case_sample'].isin(tumorsamplesincohort)].index.tolist()
    if len(cohortsamples) > 0:
      cohorts_in_batch.append(val)
      try:
        terra.addToSampleSet(wto, val, cohortsamples)
      except KeyError:  # we may not have this set yet
        wto.update_sample_set(val, cohortsamples)
    if len(pairsamples) > 0:
      cohorts_with_pairs.append(val)
      try:
        terra.addToPairSet(wto, val, pairsamples)
      except KeyError:  # we may not have this set yet
        wto.update_pair_set(val, pairsamples)
  print("creating a sample set")
  wto.update_sample_set(sample_set_id=samplesetname + "_all", sample_ids=sample_info.index.tolist())
  wto.update_sample_set(sample_set_id=samplesetname + "_tumors", sample_ids=tumorsid)
  wto.update_sample_set(sample_set_id=samplesetname + "_normals", sample_ids=normalsid)
  normalsid.extend([k for k, val in refsamples.iterrows() if val.sample_type == "Normal"])
  # Same as cum pon but better
  wto.update_sample_set(sample_set_id="All_normals", sample_ids=normalsid)
  all_samples = wto.get_samples()..index.tolist()
  all_samples.remove('NA')
  wto.update_sample_set(sample_set_id="All_samples", sample_ids=all_samples)


def submit(samplesetname):
  print("Creating Terra submissions: remember you can only cancel \
      or interact with terra submissions from the Terra website. \
      https://app.terra.bio/#workspaces/" + proc_namespace.replace(" ", "%20") + "/" + proc_workspace.replace(" ", "%20") + "/job_history")

  RenameBAM_TWIST = wto.create_submission("RenameBAM_TWIST", samplesetname + "_all", 'sample_set', expression='this.samples')
  print("waiting for 'Rename'")
  terra.waitForSubmission(wto, [RenameBAM_TWIST])

  CalculateTargetCoverage_PANCAN = wto.create_submission('CalculateTargetCoverage_PANCAN', samplesetname + "_all", 'sample_set', expression='this.samples')
  DepthOfCov_PANCAN = wto.create_submission('DepthOfCov_PANCAN', samplesetname + "_all", 'sample_set', expression='this.samples')
  print("waiting for 'CalculateTargetCoverage_PANCAN' & 'DepthOfCov_PANCAN'")
  terra.waitForSubmission(wto, [CalculateTargetCoverage_PANCAN, DepthOfCov_PANCAN])

  CreatePanelOfNormalsGATK_PANCAN = wto.create_submission('CreatePanelOfNormalsGATK_PANCAN', 'All_normals')
  DepthOfCovQC_PANCAN = wto.create_submission('DepthOfCovQC_PANCAN', samplesetname + "_all", 'sample_set', expression='this.samples')
  print("waiting for 'DepthOfCovQC_PANCAN' & 'CNV_CreatePoNForCNV'")
  terra.waitForSubmission(wto, [DepthOfCovQC_PANCAN, CreatePanelOfNormalsGATK_PANCAN])

  CallSomaticCNV_PANCAN = wto.create_submission('CallSomaticCNV_PANCAN', samplesetname + "_all", 'sample_set', expression='this.samples')
  print("waiting for 'CallSomaticCNV_PANCAN'")
  terra.waitForSubmission(wto, [CallSomaticCNV_PANCAN])

  MutationCalling_Normals_TWIST = wto.create_submission("MutationCalling_Normals_TWIST", samplesetname + "_normals", 'sample_set', expression='this.samples')
  print("waiting for 'MutationCalling_Normals_TWIST'")
  terra.waitForSubmission(wto, [MutationCalling_Normals_TWIST])

  FilterGermlineVariants_NormalSample_TWIST = wto.create_submission('FilterGermlineVariants_NormalSample_TWIST', samplesetname + "_normals", 'sample_set', expression='this.samples')
  print("waiting for 'SNV_FilterGermline'")
  terra.waitForSubmission(wto, [FilterGermlineVariants_NormalSample_TWIST])

  CreatePoNSNV_Mutect1 = wto.create_submission('CreatePoNSNV_Mutect1', "All_normals")
  CreatePoN_SNV_MuTect2 = wto.create_submission('CreatePoN_SNV_MuTect2', "All_normals")
  print("waiting for 'CreatePoN_SNV_MuTect2' & 'CreatePoNSNV_Mutect1'")
  terra.waitForSubmission(wto, [CreatePoNSNV_Mutect1, CreatePoN_SNV_MuTect2])

  PlotSomaticCNVMaps_PANCAN = wto.create_submission('PlotSomaticCNVMaps_PANCAN', samplesetname + "_all")
  for val in cohorts_in_batch:
    wto.create_submission("PlotSomaticCNVMaps_PANCAN", val)
  print("submitted final jobs for CNV pipeline")

  SNV_PostProcessing_Normals = wto.create_submission('SNV_PostProcessing_Normals', samplesetname + "_normals")
  MutationCalling_Tumors_TWIST = wto.create_submission('MutationCalling_Tumors_TWIST', samplesetname + '_pairs', 'pair_set', expression='this.pairs')
  print("waiting for 'SNV_PostProcessing' & 'MutationCalling_Tumors_TWIST'")
  terra.waitForSubmission(wto, [SNV_PostProcessing_Normals, MutationCalling_Tumors_TWIST])

  FilterGermlineEvents_TumorSample = wto.create_submission('FilterGermlineEvents_TumorSample', samplesetname + '_pairs', 'pair_set', expression='this.pairs')
  print("waiting for 'FilterGermlineEvents_TumorSample'")
  terra.waitForSubmission(wto, FilterGermlineEvents_TumorSample)

  SNVPostProcessing_TWIST = wto.create_submission('SNVPostProcessing_TWIST', samplesetname + '_pairs', "pair_set")
  print("Submitted final jobs for SNV pipeline")

  FNG_Compile_Pileup_Cnt = wto.create_submission("FNG_Compile_Pileup_Cnt", samplesetname + "_all", 'sample_set', expression='this.samples')
  print("waiting for 'FNG_Compile_Pileup_Cnt'")
  terra.waitForSubmission(wto, [FNG_Compile_Pileup_Cnt])

  FNG_Compile_db_slow_download = wto.create_submission("FNG_Compile_db_slow_download", "All_samples")
  print("waiting for 'FNG_Compile_db'")
  terra.waitForSubmission(wto, [FNG_Compile_db_slow_download])

  FNG_Query_db = wto.create_submission("FNG_Query_db", samplesetname + "_all")
  print("Submitted final FNG Job")

  print('Done')
