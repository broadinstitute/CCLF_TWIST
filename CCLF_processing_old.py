import dalmatian as dm
import pandas as pd
import Helper
import os
import numpy as np
import TerraFunction as terra
import ipdb


def getCNHeatMap(workspace1="CCLF_Targeted", namespace1="nci-mimoun-bi-org",
								 mergeCNworkflow = 'PlotSomaticCNVMaps_PANCAN',
								 output_filename = 'copy_number_heatmap.png',
								 list_sample_ids = None,
								 samplesetname = None,

								 outputloc= None
								 ):
	#####
	### NEED TO COMPLETE / FIX ###
	#####
	# produce a merged copy number heatmap for all the samples in list_sample_ids
	wm1 = dm.WorkspaceManager(namespace1, workspace1)
	# wm1.update_sample_set(samplesetname, list_sample_ids) # (check: what if set already exists?) create new sample set containing the samples in the list
	terra.addToSampleSet(namespace1+'/'+workspace1, samplesetname, list_sample_ids)
	pathto_cnvpng = "cnv_calls_tn_img"

	# run the PlotSomaticCNVMaps workflow to produce a merged, labeled copy number heat map
	submission_id = wm1.create_submission(mergeCNworkflow, samplesetname)
	# submission_id = wm1.create_submission(mergeCNworkflow, samplesetname,'sample_set',expression='this.samples')
	print("Creating merged CN heat map...")
	terra.waitForSubmission(namespace1+'/'+workspace1, submission_id)

	# copy the merged CN heat map plot to the output directory
	sample_set = wml.get_sample_sets()
	cnfile = sample_set["sample_set_id"==samplesetname][pathto_cnvpng]
	os.system('gsutil cp ' + cnfile + ' ' + outputloc + output_filename)
	print("Copied the merged CN heat map to the output location.")




def getReport(workspace1="CCLF_TSCA_2_0_3_HCMI", namespace1="nci-mimoun-bi-org", # CCLF_Targeted instead of CCLF_TSCA_2_0_3_HCMI or 2_0_2
							pathto_cnvpng='segmented_copy_ratio_img',
							pathto_snv='filtered_variants',
							pathto_seg='cnv_calls',
							workspacewes='CCLF_WES', namespacewes='nci-mimoun-bi-org',
							pathto_cnvpng_wes='segmented_copy_ratio_img',
							pathto_snv_wes='mafLite',
							pathto_seg_wes='tumor_seg',
							is_from_pairs=True,
							datadir='gs://cclf_results/targeted/',
							tempdir='temp/cclfmerge/',
							specificlist=None,
							specificlist_disease=None
							):
	"""
	This function gathers all TSV and PNG files related to mutations and CNVs for the specified participants,
	annotates the images and adds metadata column to the TSVs as appropriate, and outputs the results to a single
	Google bucket. Currently, the Google bucket is organized into subfolders by disease type
	(i.e. all Neuroblastoma participants will be in the Neuroblastoma folder, and each participant will have its own folder).

	# Args
	- specificlist = filename/filepath of file containing a list of participant IDs to acquire data for


	# Args specific to WES:
	- workspacewes
	- namespacewes
	- pathto_cnvpng_wes
	- pathto_snv_wes
	- pathto_seg_wes

	# Args specific to targeted sequencing:
	- workspace1
	- namespace1
	- pathto_cnvpng
	- pathto_snv
	- pathto_seg

	# Overview:
	For each participant, we look for all information we can find in each workspace / for each datatyle (WES, TSCA, TWIST).
	1.

	# TODO: instead, pass in dict of workspace:dataset? e.g.
	{'nci-mimoun-bi-org/CCLF_TSCA_2_0_3_HCMI':'TSCA', 'nci-mimoun-bi-org/CCLF_WES':'WES'}
	That way, could loop through the keys in the dictionary.

	Or we could also make a dict of dicts so that we could have
	the variables all be called the same thing! (e.g. get rid of pathto_seg vs pathto_seg_wes):
	{'nci-mimoun-bi-org/CCLF_TSCA_2_0_3_HCMI':
		{dataset:'TSCA', 'pathto_cnvpng':'segmented_copy_ratio_img', 'pathto_snv':'filtered_variants'},
	'nci-mimoun-bi-org/CCLF_WES':
		{dataset:'WES', 'pathto_cnvpng':'segmented_copy_ratio_img', 'pathto_snv':'mafLite'}
	}
	"""
	# Create and update source info dictionary (dict of dicts)
	source_info = {
	'TSCA':
		{'workspace':'nci-mimoun-bi-org/CCLF_TSCA_2_0_3_HCMI',
		'pathto_cnvpng':'segmented_copy_ratio_img',
		'pathto_snv':'filtered_variants',
		'pathto_seg':'cnv_calls',
		'external_id_col':'external_id_validation',
		'is_targeted':True},
	'TWIST':
		{'workspace':'nci-mimoun-bi-org/PANCAN_TWIST copy',
		'pathto_cnvpng':'segmented_copy_ratio_img',
		'pathto_snv':'filtered_variants',
		'pathto_seg':'cnv_calls',
		'external_id_col':'external_id_validation',
		'is_targeted':True},
	'WES':
		{'workspace':'nci-mimoun-bi-org/CCLF_WES',
		'pathto_cnvpng':'segmented_copy_ratio_img',
		'pathto_snv':'mafLite', # mafLite in pairs TSV, not in sample TSV
		'pathto_seg':'cnv_calls', # NOT tumor_seg
		'external_id_col':'external_id_capture',
		'is_targeted':False}
	}

	print('you need to be on macOS for now')
	# Add key data to source_info dict of dicts
	for key in source_info.keys():
		# Add workspace
		source_info[key]['wm'] = dm.WorkspaceManager(source_info[key]['workspace'])
		# Add samples
		source_info[key]['samples'] = source_info[key]['wm'].get_samples()
		# Add pairs
		if is_from_pairs:
			source_info[key]['pairs'] = source_info[key]['wm'].get_pairs()
		# Add participants
		source_info[key]['participants'] = source_info[key]['samples'].participant.tolist()

	if type(specificlist) is str:
		# we consider it is a filename
		specificlist = pd.read_csv(specificlist).tolist()
	elif specificlist is None:
		print("taking all samples from all the data workspaces combined")
		specificlist = []
		# TODO: check if want get_sample sets, or more likely just the list of participants or samples?
		for i in source_info.keys():
			specificlist.append(source_info[key]['wm'].get_sample_sets().index)
		specificlist = list(set(specificlist)) # TODO: better way to remove duplicates from the list?

	###################
	# get data for each participant ID (e.g. PEDS182) by looking through each dataset
	###################
	# TODO: go through and identify portions that are specific to targeted only (e.g. depth of cov QC filtering, I believe)
	for val in specificlist:
		# Extract the primary disease information and include in the output Google storage directory
		primary_disease = pd.read_csv(specificlist_disease, index_col = "participant_id").loc[val,"primary_disease"]
		outputloc = datadir + primary_disease.replace(' ', '_') + '/' + val + '/'

		found = False
		images = [] # store all the CN images for a participant

		all_ext_ids = pd.DataFrame(columns = ['external_id', 'sample_id', "dataset", "is_targeted", "condition", "is_normal", "has_mut", "has_cnv"])
		all_failed_ext_ids = pd.DataFrame(columns = ['external_id', 'sample_id', "dataset", "is_targeted", "condition", "is_normal"])

		for dset in source_info.keys():
			# Define some variables here to make the following code more readable
			samples = source_info[dset]['samples']
			pairs = source_info[dset]['pairs']
			participants = source_info[dset]['participants']
			pathto_seg = source_info[dset]['pathto_seg']
			pathto_snv = source_info[dset]['pathto_snv']
			pathto_cnvpng = source_info[dset]['pathto_cnvpng']
			external_id_col = source_info[dset]['external_id_col']
			is_targeted = source_info[dset]['is_targeted']
			# dataset = source_info[dset]['dataset']

			# Skip  to next participant if this one isn't in the dset's  participants
			if val not in participants:
				continue  # TODO: continue or break? I think it might be a break because I want to skip back to the next value. Maybe I should change the order of for loops for participant vs dset?

			print('Getting data for {} from {}...'.format(str(val), str(dset)) )
			# TODO: FIXME: these need to be inside a for loop or two. Need to fix this.
			source_info[dset]['found'] = False # found_TSCA
			source_info[dset]['ext_id_df'] = pd.DataFrame(columns=["external_id", "dataset", "is_targeted"])
			source_info[dset]['ext_ids'] = [] # ext_ids
			source_info[dset]['sample_ids'] = [] # sample_ids
			# ext_ids = source_info[dset]['ext_ids']
			# sample_ids = source_info[dset]['sample_ids']
			# ext_id_df = source_info[dset]['ext_id_df']
			# ext_id_df.append({'external_id':__, "dataset":__, "is_targeted":__})


			mutfile = pd.DataFrame()
			cnfile = pd.DataFrame()

			sample_subset = samples[samples.participant == val]
			# source_info[dset]['sample_ids'] = [] # sample_ids
			# sample_ids = source_info[dset]['sample_ids']
			sample_ids = [] # will collect all the sample_ids from the targeted data; sample set to make CN heat map for
			normals_for_participant = []
			# ext_ids = []
			for k, condition in sample_subset.iterrows():
				if condition['sample_type'] == "Tumor":
					external_id = condition[external_id_col]
					# TODO: do we have a depth of cov QC for WES samples? Is it built in to a pipeline?
					# Add external IDs that fail the depth of coverage QC to a list & skip further analysis
					if is_targeted and condition['depth_of_cov_qc_result'] == 'fail':
						all_failed_ext_ids = all_failed_ext_ids.append({'external_id':external_id, 'sample_id':k, "dataset":dset, "is_targeted":is_targeted, "condition":condition["media"], "is_normal":False}, ignore_index=True)
						continue
					# Track metadata (sample ID, external ID, whether we found this sample in the dataset)
					found = True
					source_info[dset]['found'] = True

					# Condition is NA for WES, and will be updated below if targeted sequencing
					cond_name = np.nan
					if is_targeted:
						cond_name = condition['media']
						all_ext_ids = all_ext_ids.append(pd.Series({'external_id':external_id, 'sample_id':k, "dataset":dset, "is_targeted":is_targeted, "condition":condition["media"], "is_normal":False, "has_mut":False, "has_cnv":False}, name = external_id), ignore_index=False)
					else:
						all_ext_ids = all_ext_ids.append(pd.Series({'external_id':external_id, 'sample_id':k, "dataset":dset, "is_targeted":is_targeted, "condition":np.nan, "is_normal":False, "has_mut":False, "has_cnv":False}, name = external_id), ignore_index=False)

					# Copy the seg file and CN horizontal plot locally
					print("Getting seg file and CNV png for ", k,"...")
					seg_path = condition[pathto_seg]
					cnvpng_path = condition[pathto_cnvpng]
					if type(seg_path) is str: # TODO: if float, means that it was emtpy / nan
						# all_ext_ids[all_ext_ids['external_id'] == external_id]['has_cnv'] = True
						all_ext_ids.loc[external_id,'has_cnv'] = True
						os.system('gsutil cp ' + seg_path + ' ' + tempdir + 'copy_number.tsv')
					if type(cnvpng_path) is str:
						os.system('gsutil cp ' + cnvpng_path + ' ' + tempdir + external_id + '_copy_number_map.png')
						# Annotate local image with sample information and copy to output Google bucket (outputloc)
						# Text added to image is the external ID followed by the dataset (TWIST, TSCA, or WES) in parentheses
						imagedir = tempdir + external_id + '_copy_number_map.png'
						text = condition[external_id_col].replace("_", " ") + " (" + dset + ")"
						Helper.addTextToImage(imagedir, text, outputpath=imagedir)
						images.append(imagedir)
						# TODO: start here. This should print for WES samples too, right?
						print("Added imagedir to images")
						print("imagedir:", imagedir)
						print("images:", images)
						os.system('gsutil cp ' + imagedir + ' ' + outputloc + dset + "_" + external_id + '_copy_number_map.png')

					if is_from_pairs:
						pair_subset = pairs[pairs["case_sample"] == k]

						print("Getting tumor SNV info, plus matched normal SNV if it exists for case sample", k)

						# The snv file w/filtered variants is stored in different locations (sample vs pairs TSV, various columns) depending on the dset and sample type
						# However, here we are specifically working with tumor samples, which are always stored in the pairs TSVs
						snvs = pair_subset[pathto_snv]
						for snv in snvs:
							if snv is not np.nan:
								# all_ext_ids[all_ext_ids['external_id'] == external_id]['has_mut'] = True
								all_ext_ids.loc[external_id,'has_mut'] = True
								os.system('gsutil cp ' + snv + ' ' + tempdir + 'mutations.tsv')
								# add first matched normal to list if exists (aka add the normal's sample_id)
								# TODO: what to do if have multiple matched normals?
								try:
									# TODO: I not certain I'm tracking this sample ID properly. Does this sample ID make it into the overall TSV?
									sample_ids += [next(id for id in pair_subset["control_sample"] if id not in ["NA",np.nan])]
									if sample_ids[-1] not in set(all_ext_ids.loc[:,"sample_ids"]):
										raise Exception("The normal sample", sample_ids[-1], " is missing from the all_ext_ids TSV.")
									print('Getting pair data for %s...' % str(val))
								# only NA or nans; no real matched normal
								except:
									continue
								break

						# If a matched normal exists for the current participant and external ID,
						# copy the normal CN plot locally, annotate the image, then copy to output location
						# TODO: currently assuming only 1 distinct matched normal exists.
						print("If a matched normal exists for tumor sample", k,", copy the normal CN plot locally, annotate the image, then copy to output location")
						matched_normal_ids = pair_subset["control_sample"]
						print("matched_normal_ids:",matched_normal_ids)
						for normal_id in matched_normal_ids:
							if normal_id not in [np.nan, "null", "NA", ""] and normal_id not in normals_for_participant:
								normals_for_participant += [normal_id]
								print("Adding a matched normal,", normal_id,", for sample",k)
								# Add the matched normal sample to the external ID TSV (if it passes QC),
								# including the external ID, dataset (TWIST, TSCA, WES), and label it as a matched normal
								normal_sample = samples[samples.index == normal_id]
								normal_external_id = normal_sample[external_id_col][normal_id]
								normal_cond_name = np.nan
								if is_targeted:
									normal_cond_name = normal_sample['media'][normal_id]

								print("Adding matched normal to summary tables")
								if is_targeted and normal_sample.loc[normal_id,'depth_of_cov_qc_result'] == 'fail':
									all_failed_ext_ids = all_failed_ext_ids.append({'external_id':normal_external_id, 'sample_id':normal_id, "dataset":dset, "is_targeted":is_targeted, "condition":normal_sample["media"][normal_id], "is_normal":True}, ignore_index=True)
								elif is_targeted:
									# TODO: for normal samples, use None or "NA" for has_mut?
									all_ext_ids = all_ext_ids.append(pd.Series({'external_id':normal_external_id, 'sample_id':normal_id, "dataset":dset, "is_targeted":is_targeted, "condition":normal_sample["media"][normal_id], "is_normal":True, "has_mut":None, "has_cnv":False}, name = normal_id), ignore_index=False)
								else:
									all_ext_ids = all_ext_ids.append(pd.Series({'external_id':normal_external_id, 'sample_id':normal_id, "dataset":dset, "is_targeted":is_targeted, "condition":np.nan, "is_normal":True, "has_mut":None, "has_cnv":False}, name = normal_id), ignore_index=False)

								print("Getting seg file and CNV png for the matched normal", normal_id,"...")
								normal_seg_path = normal_sample[pathto_seg][normal_id]
								normal_cnvpng_path = normal_sample[pathto_cnvpng][normal_id]
								if type(seg_path) is str: # TODO: if float, means that it was emtpy / nan
									all_ext_ids.loc[external_id,'has_cnv'] = True
									# all_ext_ids[all_ext_ids['external_id'] == external_id]['has_cnv'] = True
									os.system('gsutil cp ' + normal_seg_path + ' ' + tempdir + 'normal_copy_number.tsv')

								# If available, annotate matched normal CN image with the external ID, dataset (TWIST, TSCA, WES), and label it as a matched normal
								if normal_sample[pathto_cnvpng][normal_id] not in ["NA", np.nan]:
									print("Annotating matched normal CN image with the external ID, dataset (TWIST, TSCA, WES), labeling it as a matched normal, and uploading")
									# Copy image locally
									imagedir = tempdir + normal_id + '_copy_number_map.png'
									os.system('gsutil cp ' + normal_sample[pathto_cnvpng][normal_id] + ' ' + imagedir)

									if os.path.exists(imagedir):
										text = normal_external_id.replace("_", " ") + " (" + dset + ", matched normal)"
										Helper.addTextToImage(imagedir, text, outputpath=imagedir)
										images.append(imagedir)
										print("Added path for normal sample to images")
										print("path:", str(imagedir))
										print("images:", images)
										# Copy to output location
										os.system('gsutil cp ' + imagedir + ' ' + outputloc + dset + "_" + normal_external_id + '_copy_number_map.png')
									break
					# If not is_from_pairs
					else:
						# TODO: currently, pathto_snv_unpaired not implemented. I'm not sure which column this would be in Terra.
						os.system('gsutil cp ' + condition[pathto_snv_unpaired] + ' ' + tempdir + 'mutations.tsv')

					# Add culture condition, external_id, and dataset (TWIST, TSCA, WES) columns to the tables
					# TODO: determine whether to combine targeted and WES into one table, or if have different colnames whether we should split it into two (targeted vs not)
					mut = pd.read_csv(tempdir + 'mutations.tsv', sep='\t', index_col = None)
					mut['condition'] = cond_name
					mut['external_id'] = external_id
					mut['data source'] = dset
					mutfile = mutfile.append(mut)

					cn = pd.read_csv(tempdir + 'copy_number.tsv', sep='\t', index_col = None)
					cn['condition'] = cond_name
					cn['external_id'] = external_id
					cn['data_source'] = dset
					cn['is_normal'] = False
					cnfile = cnfile.append(cn)

					normal_cn = pd.read_csv(tempdir + 'normal_copy_number.tsv', sep='\t', index_col = None)
					normal_cn['condition'] = normal_cond_name
					normal_cn['external_id'] = normal_external_id
					normal_cn['data_source'] = dset
					normal_cn['is_normal'] = True
					cnfile = cnfile.append(normal_cn)
			if found:
				# Write the CN and mutation files with added metadata locally, then copy to the output location
				cnfile.to_csv(tempdir + 'cn.tsv', sep='\t', index = False)
				mutfile.to_csv(tempdir + 'mut.tsv', sep='\t', index = False)
				os.system('gsutil cp ' + tempdir + 'cn.tsv ' + outputloc + dset + '_copy_number.tsv')
				os.system('gsutil cp ' + tempdir + 'mut.tsv ' + outputloc + dset + '_mutation.tsv')


			# At this level, aka for each dataset, need to update the DF for the external IDs and failed external ext_ids
			print("Looking at the results after looking through ", dset, " data:")
			print("all_ext_ids:", all_ext_ids)
			print("all_failed_ext_ids:", all_failed_ext_ids)

		# TODO: fix from here down.
		# Upload dataframe with metadata for each external IDs, including the source dataset
		print("Uploading dataframe with metadata for each external IDs, including the source dataset")
		all_ext_ids['disease'] = primary_disease
		all_ext_ids['participant'] = val
		all_ext_ids.to_csv(tempdir + 'all_external_ids.tsv', sep='\t', index=False)
		os.system('gsutil cp ' + tempdir + 'all_external_ids.tsv ' + outputloc + 'all_external_ids.tsv')

		# TODO: currently, just checking for targeted because don't have the same QC for WES.
		# Upload dataframe with metadata for each failed external ID, including the source dataset
		print("Uploading dataframe with metadata for each failed external ID, including the source dataset")
		all_failed_ext_ids['disease'] = primary_disease
		all_failed_ext_ids['participant'] = val
		all_failed_ext_ids.to_csv(tempdir + 'all_failed_external_ids.tsv', sep='\t', index=False)
		os.system('gsutil cp ' + tempdir + 'all_failed_external_ids.tsv ' + outputloc + 'all_failed_external_ids.tsv')

		# # Upload all images stored in images to the google bucket
		# for img_path in images: <- or I could just do this as I add to the list 'images'. I think I already do this :)
		#	 os.system('gsutil cp ' img_path + ' ' + outputloc + 'something_here.tsv')

		if not found:
			print("We did not find any targeted probe data or WES data for " + val)
		else:
			# merge all the horizontal CNV plots
			print("Uploading the merged CNV plot for", val)
			print("images list:", images)
			Helper.mergeImages(sorted(images), tempdir + 'merged.png')
			os.system('gsutil cp ' + tempdir + 'merged.png ' + outputloc + 'merged_copy_number_map.png')



# def get_snv_path(dset, external_id, sample_type, source_info): # maybe also participant id?
# 	'''For targeted data, samples[“filtered_variants”] if working with a normal sample.
# 	Otherwise, pairs[“filtered_variants”] if tumor (regardless if has matched normal)
# 	For WES data, I don’t see any column in the samples TSV.
# 	Looks like just pairs[“mafLite”] if it has been processed.'''
# 	if dset == "WES":
# 		# SNV info is in pairs dset
# 		# pairs[pathto_snv]
#
# 	# Targeted data currently stores the SNV information in two different locations
# 	else:
# 		if sample_type == "Normal":
# 			# samples[pathto_snv]
# 		else:
# 			# pairs[pathto_snv]
