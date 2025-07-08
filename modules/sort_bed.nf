#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl = 2

process SORT_BED {
    container "${params.container__bedtools}"
    label 'single_proc'
    tag "$bed_name"

    input:
        tuple val(bed_name), path(bed_path)
        
    output:
        tuple val(bed_name), path("*_sorted.bed")

    script:
        """
        out="${bed_name}_sorted.bed"
        bedtools sort -i ${bed_path} | sed -r 's/^chr//g' > \${out}
        """
}
