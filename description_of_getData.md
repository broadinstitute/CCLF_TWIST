What do we need to get for each participant?
1. external IDs, sample IDs, media condition, is_targeted, disease type
2. metadata for each external ID:
  a. tumor vs normals
  b. find matched normals if they exist
3. seg file path
4. cnv img path
5. snv / maf file path
6. for each file path, whether the link exists, and whether it is broken
7. pass/fail depth of coverage QC (although I don't know who actually uses this metric! We certainly should not use it to decide which samples to use mutation data from. GATK explicitly discourages this, and says that its models are sophisticated enough to take )


What do we need to create for each participant?
(some of this will come from the jupyter notebook)
1. annotated CN images
2. merged CN image
3. CN heatmap
4. mutation summary table
5. a way to incorporate Moony's decisions (verified tumor/normal; need to extract data from the external sheets)

Ideally, create dictionary for each participant, then check the google storage paths systematically and in a parallelized mode.
{participant:
  {
    extID1: {seg_file_path: dfnadf,
             cnv_img_path: adsfaf,
             snv_path: asdfart},
    extID2: {seg_file_path: edfafnadf,
            cnv_img_path: eaadsfaf,
            snv_path: edfart},
    ...
  }
}

Create separate dict for identifying broken paths, mapping path to status.

After get status, need to update the participant dictionary to include things like:
1. has_cnv, has_snv, has_tpm, has_fusions, etc. when appropriate

Function that takes sample tsv(s), participant, returns the second dict level.

for workspace in workspaces:
  result = getMetadata(participant, workspace_sample_tsv)
  add result to participant dictionary



-------
Will take a long time and lots of effort to functionalize. For now, might just identify all gsutil steps and see if I can combine/simplify to reduce time.
