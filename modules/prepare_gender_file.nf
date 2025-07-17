#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl = 2

process PREPARE_GENDER_FILE {
    container "${params.container__bcftools}"
    label 'single_proc'

    input:        
        path gender_files
        
    output:
        path ("gender.txt")

    script:                
        """
            cat ${gender_files.join(' ')} > gender.txt
        """
}
