# Given an old FNG database and a new table with FNG data, combine into a single table
workflow merge_FNG_databases_wrkflw {
    call merge_FNG_databases
}

task merge_FNG_databases {
    File old_FNG_database
    File new_FNG_data
    String sample_set_id
    File merge_FNG_databases_script

    Int memory
    Int disk_space
    Int num_preempt


    command {
        Rscript ${merge_FNG_databases_script} \
            --sample_set_id ${sample_set_id} \
            --old_fng_db ${old_FNG_database} \
            --new_fng_db ${new_FNG_data}
    }

    output {
        File combined_FNG_databases_file="fingerprinting_db_through_${sample_set_id}.txt"
        File combined_FNG_databases_for_workspace_tracking="fingerprinting_db_through_${sample_set_id}.txt"
    }

    runtime {
        docker: "rocker/tidyverse:3.6.1"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Gwen Miller"
    }
}
