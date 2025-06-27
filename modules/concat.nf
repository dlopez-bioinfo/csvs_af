#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl = 2

process CONCAT {
    container "${params.container__bcftools}"
    label 'high_proc'

    publishDir "${params.output_folder}/", mode: 'copy', overwrite: true

    input:
        path(vcf_list)
        
    output:
        tuple path("csvs.vcf.gz"), path("csvs.vcf.gz.csi")

    script:
        def out = "csvs.vcf.gz"

        """
        bcftools concat --threads ${task.cpus} *.vcf.gz | bcftools sort -o ${out} -O z
        bcftools index ${out}
        """
}
