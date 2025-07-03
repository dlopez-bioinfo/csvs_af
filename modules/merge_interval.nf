#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl = 2

process MERGE_INTERVAL {
    container "${params.container__bcftools}"
    label 'single_proc'
    tag "$interval"

    input:
        tuple val(interval), path(vcf_list), path(vcf_idx)
        path ref_genome
        path ref_genome_fai

        
    output:
        tuple path("csvs_*.vcf.gz"), path("csvs_*.vcf.gz.csi")

    script:
        def out = "csvs_${interval}.vcf.gz"
        def vcf_files = vcf_list.collect { it.getName() }
        def region = interval.replaceAll('_', '-').replaceFirst('-', ':')

        if (vcf_files.size() == 1) {
            """
            bcftools view ${vcf_files[0]} -o ${out} -O z -G -r ${region}
            bcftools index ${out}
            """
        } else {
            """
            bcftools merge ${vcf_files.join(' ')} -0 -r ${region} | \\
                bcftools norm -D -d both -f ${ref_genome} | \\
                bcftools view -o ${out} -O z -G
            bcftools index ${out}
            """
        }
}
