workflow PlotSomaticCNVMaps {
    String tsca_id
    String disk_size
    String docker_img
    String num_cpu
    Array[File] tumor_tn_files
    Array[File] tumor_seg_files
    Array[File] tumor_seg_files_for_plotting
    Array[String] sample_types
    Array[String] sample_ids
    Array[String] participant_ids
    Array[String] external_validations
    Array[String] depth_of_cov_qc_results
    
    call plotSomaticCNVMaps {
        input: 
            tsca_id=tsca_id,
            docker_img=docker_img,
            tumor_tn_files=tumor_tn_files,
            tumor_seg_files=tumor_seg_files,
            sample_ids=sample_ids,
            participant_ids=participant_ids,
            external_validations=external_validations,
            depth_of_cov_qc_results=depth_of_cov_qc_results,
            sample_types=sample_types,
            tumor_seg_files_for_plotting=tumor_seg_files_for_plotting,
            disk_size=disk_size,
            num_cpu=num_cpu
    }

    output {
        File cnv_calls_tn_img_unsegmented 						= plotSomaticCNVMaps.cnv_calls_tn_img_unsegmented
        File cnv_calls_tn_raw_unsegmented 						= plotSomaticCNVMaps.cnv_calls_tn_raw_unsegmented
		File cnv_calls_tn_raw_unsegmented_sample_ids 			= plotSomaticCNVMaps.cnv_calls_tn_raw_unsegmented_sample_ids
		File cnv_calls_tn_img_unsegmented_aligned 				= plotSomaticCNVMaps.cnv_calls_tn_img_unsegmented_aligned
		File cnv_calls_tn_raw_unsegmented_sample_ids_aligned 	= plotSomaticCNVMaps.cnv_calls_tn_raw_unsegmented_sample_ids_aligned
    }
}


task plotSomaticCNVMaps {
    String tsca_id
    String disk_size
    String num_cpu
    String docker_img
    Array[File] tumor_tn_files
    Array[File] tumor_seg_files
    Array[File] tumor_seg_files_for_plotting
    Array[String] sample_types
    Array[String] participant_ids
    Array[String] sample_ids
    Array[String] external_validations
    Array[String] depth_of_cov_qc_results
    
    command <<<
        cd /home
        git clone https://github.com/broadinstitute/CCLF_TWIST.git

        #### Plot unsegmented CNV Calls
        python CCLF_TWIST/workflows/scripts/plot_somatic_cnv_calls.py  --tsca_id ${tsca_id} \
                                                --tn_files ${sep=" " tumor_tn_files} \
                                                --external_ids ${sep=" " external_validations} \
                                                --sample_ids ${sep=" " sample_ids} \
                                                --participant_ids ${sep=" " participant_ids} \
                                                --sample_types ${sep=" " sample_types} \
                                                --depth_of_cov_qcs ${sep=" " depth_of_cov_qc_results}

		#### EXPERIMENTAL: Plot CNV calls with aligned intervals
		python CCLF_TWIST/workflows/scripts/plot_somatic_cnv_calls_align_intervals.py  \
												--tsca_id ${tsca_id} \
		                                        --tn_files ${sep=" " tumor_tn_files} \
		                                        --external_ids ${sep=" " external_validations} \
		                                        --sample_ids ${sep=" " sample_ids} \
		                                        --depth_of_cov_qcs ${sep=" " depth_of_cov_qc_results}                                            
        
    >>>
    
    output {
        File cnv_calls_tn_img_unsegmented 						= "${tsca_id}.cnv_calls_unsegmented.png"
        File cnv_calls_tn_raw_unsegmented 						= "${tsca_id}.cnv_calls_unsegmented.txt"
        File cnv_calls_tn_raw_unsegmented_sample_ids 			= "${tsca_id}.cnv_calls_unsegmented.sample_ids.txt"

        File cnv_calls_tn_img_unsegmented_aligned 				= "${tsca_id}.cnv_calls_unsegmented.aligned.png"
        File cnv_calls_tn_raw_unsegmented_sample_ids_aligned 	= "${tsca_id}.cnv_calls_unsegmented.sample_ids.aligned.txt"
    }
    
    runtime {
        docker: "${docker_img}"
        memory: "${disk_size}"
        cpu: "${num_cpu}"
        disks: "local-disk 40 HDD"
    }
}