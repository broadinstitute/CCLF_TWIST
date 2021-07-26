# CCLF QC pipeline
# ISSUE: some tasks need to run on a single batch, but others need to run on all the new data combined.

import "https://raw.githubusercontent.com/broadinstitute/ccle_processing/master/RNA_pipeline/samtofastq_v1-0_BETA_cfg.wdl" as samtofastq_v1
import "https://raw.githubusercontent.com/broadinstitute/ccle_processing/master/RNA_pipeline/star_v1-0_BETA_cfg.wdl" as star_v1
import "https://raw.githubusercontent.com/broadinstitute/ccle_processing/master/RNA_pipeline/rnaseqc2_v1-0_BETA_cfg.wdl" as rnaseqc2_v1
import "https://raw.githubusercontent.com/broadinstitute/ccle_processing/master/RNA_pipeline/rsem_v1-0_BETA_cfg.wdl" as rsem_v1
import "https://raw.githubusercontent.com/broadinstitute/ccle_processing/master/RNA_pipeline/hg38_STAR_fusion.wdl" as hg38_STAR_fusion


workflow QC_pipeline {

  # Inputs sample level
  String sample_id
	String external_id

  # Inputs sample set level
  String sample_set_id # TODO: how deal with per sample vs per batch analyses? Should be fine, actually: eg this.samples.target_coverage vs this.target_coverage
  # needed for PlotSomaticCNVMaps_PANCAN
  Array[String] sample_types
  Array[String] sample_ids
  Array[String] participant_ids
  Array[String] external_validations

  # renameBAM_TWIST
  File bam_file
  File bai_file

  # CalculateTargetCoverage_PANCAN
  File gatk
  File interval_list
  Int padding

  # DepthofCoverage
  File ref_fasta
	File gene_list
	File ref_fasta_index
	File ref_fasta_dict

  # TODO: CreatePanelOfNormalsGATK_PANCAN (need to edit the config)
  File gatk_jar




  call RenameBAM_TWIST.renameBAM as renameBAM {
    input:
      sample_id=sample_id,
      bam_file=bam_file,
      bai_file=bai_file
  }

  call CalculateTargetCoverage_PANCAN.CalculateTargetCoverage as CalculateTargetCoverage {
  	input:
    	sample_id=sample_id,
      gatk=gatk,
      input_bam=renameBAM.renamed_bam,
      interval_list=interval_list,
      padding=padding
    }

  call DepthOfCov_PANCAN.depthOfCov as depthOfCov {
  	input:
      inputBam=renameBAM.renamed_bam,
      inputBamIndex=renameBAM.renamed_bai,
      sampleName=sample_id,
      refFasta=ref_fasta,
  		geneList=gene_list,
  	#	intervalList=interval_list, # TODO: NOT the same intervalList. workspace.PANCAN_target_intervals used here vs workspace.pancaner_list used elsewhere
  		refFastaIndex=ref_fasta_index,
  		refFastaDict=ref_fasta_dict
	}


  # TODO: need to update config, and deal with sample set as root entity. Can I loop over a variable, like a list of sample sets?
  call CreatePanelOfNormalsGATK_PANCAN.createPanelOfNormals as createPanelOfNormals {
    input:
        PoN_name=sample_set_id,
        gatk_jar=gatk,
        target_coverage_files=target_coverage_files #CalculateTargetCoverage.target_coverage
    }

  # TODO: need to do on a per-batch level
  call DepthOfCovQC_PANCAN.depthOfCovQC as depthOfCovQC {
  	input: refFasta=ref_fasta,
  		geneList=gene_list,
  		#intervalList=interval_list, # TODO: NOT the same intervalList. workspace.PANCAN_target_intervals used here vs workspace.pancaner_list used elsewhere
  		refFastaIndex=ref_fasta_index,
  		refFastaDict=ref_fasta_dict
  	}

  # TODO: need to edit the input config. Batch dependent. Calling a whole workflow.
  # call CallSomaticCNV_PANCAN.CallSomaticCNV as CallSomaticCNV {
  call CallSomaticCNV {
    input:
      sample_id=sample_id,
      external_id=external_id,
      gatk_jar=gatk,
      reference_fasta_dict=ref_fasta_dict,
      tumor_coverage=this.target_coverage, # TODO: sample level
      normals_pon=# TODO: this is batch specific. Output from createPanelOfNormals, e.g. "workspace.pon_normals_CCLF_TWIST34"
  }

  # TODO: each batch and each cohort. Run at sample set level.
  call PlotSomaticCNVMaps_PANCAN.plotSomaticCNVMaps as plotSomaticCNVMaps {
    input:
      tsca_id=sample_set_id,
      sample_ids=sample_ids,
      participant_ids=participant_ids,
      external_validations=external_validations,
      sample_types=sample_types,
      tumor_tn_files=this.samples.tumor_tn,
      tumor_seg_files=this.samples.tumor_seg,
      depth_of_cog_qc_results=this.samples.depth_of_cov_qc_result,
      tumor_seg_files_for_plotting=this.samples.tumor_seg_for_plotting
  }

  call PlotSomaticCNVMaps_PANCAN. as  {
    input:
  }
  MutationCalling_Normals_TWIST
  FilterGermlineVariants_NormalSample_TWIST
  # TODO: need to edit the "PoN_name" config for CreatePoNSNV_Mutect1 and CreatePoNSNV_Mutect2
  CreatePoNSNV_Mutect1
  CreatePoNSNV_Mutect2
  SNV_PostProcessing_Normals
  MutationCalling_Tumors_TWIST # TODO: (edit the input config to match pon_mutect1, pon_mutect2)
  FilterGermlineEvents_TumorSample
  SNVPostProcessing_TWIST
  FNG_Compile_Pileup_Cnt
  FNG_Compile_db_slow_download
  # TODO: add merge_FNG WDL!
  FNG_Query_db


  output {
  renameBAM
  CalculateTargetCoverage
  depthOfCov
  createPanelOfNormals
  depthOfCov_panel

  CallSomaticCNV_PANCAN # TODO: need to edit the input config
  PlotSomaticCNVMaps_PANCAN # TODO: each batch and each cohort
  MutationCalling_Normals_TWIST
  FilterGermlineVariants_NormalSample_TWIST
  # TODO: need to edit the "PoN_name" config for CreatePoNSNV_Mutect1 and CreatePoNSNV_Mutect2
  CreatePoNSNV_Mutect1
  CreatePoNSNV_Mutect2
  SNV_PostProcessing_Normals
  MutationCalling_Tumors_TWIST # TODO: (edit the input config to match pon_mutect1, pon_mutect2)
  FilterGermlineEvents_TumorSample
  SNVPostProcessing_TWIST
  FNG_Compile_Pileup_Cnt
  FNG_Compile_db_slow_download
  # TODO: add merge_FNG WDL!
  FNG_Query_db

    #samtofastq
    #star
    File bam_file=star.bam_file
    File bam_index=star.bam_index
    File transcriptome_bam=star.transcriptome_bam
    File chimeric_junctions=star.chimeric_junctions
    File chimeric_bam_file=star.chimeric_bam_file
    File read_counts=star.read_counts
    File junctions=star.junctions
    File junctions_pass1=star.junctions_pass1
    Array[File] logs=star.logs
    #rnaseqc
    File gene_tpm=rnaseqc2.gene_tpm
    File gene_counts=rnaseqc2.gene_counts
    File exon_counts=rnaseqc2.exon_counts
    File metrics=rnaseqc2.metrics
    File insertsize_distr=rnaseqc2.insertsize_distr
    #rsem
    File genes=rsem.genes
    File isoforms=rsem.isoforms
    #StarFusion
    File fusion_predictions=StarFusion.fusion_predictions
    File fusion_predictions_abridged=StarFusion.fusion_predictions_abridged
  }
}
