import "RenameBAM_TWIST.wdl" as RenameBAM_TWIST


workflow CCLF_TWIST_pipeline {

    # RenameBAM_TWIST
    String sample_id
    File bam_file
    File bai_file

    call RenameBAM_TWIST.RenameBAM as RenameBAM{
        input:
            sample_id=sample_id
            bam_file=bam_file
            bai_file=bai_file
    }



}