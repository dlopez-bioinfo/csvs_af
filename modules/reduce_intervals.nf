#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl = 2

process REDUCE_INTERVALS {
    container "${params.container__bedtools}"
    label 'single_proc'

    input:
        tuple val(bed_name), val(raw_intervals)
        
    output:
        tuple val(bed_name), path("regions.bed")

    script:
        """
            echo "${raw_intervals}" | tr ',' '\n' | awk -F"_" '{print \$1"\t"\$2"\t"\$3}' > aux
            bedtools merge -i aux > regions.bed
        """
}
