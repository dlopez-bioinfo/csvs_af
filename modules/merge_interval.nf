#!/usr/bin/env nextflow

// Using DSL-2
//nextflow.enable.dsl = 2

process MERGE_INTERVAL {
    container "${params.container__bcftools}"
    label 'single_proc'


    input:
        tuple val(raw_intervals), path(vcf_list), path(vcf_index)
        path ref_genome
        path ref_genome_fai

        
    output:
        tuple path("*_merged.vcf.gz"), path("*_merged.vcf.gz.csi")

    script:        
        def vcf_files = vcf_list.collect { it.getName() }
        def out = raw_intervals.md5() + "_merged.vcf.gz"

        if (vcf_files.size() == 1) {
            """
            echo "${raw_intervals}" | tr ',' '\n' | awk -F"_" '{print \$1"\t"\$2"\t"\$3}' > regions.bed   

            bcftools view ${vcf_files[0]} -o ${out} -O z -G -R regions.bed
            bcftools index ${out}
            """
        } else {
            """
            echo "${raw_intervals}" | tr ',' '\n' | awk -F"_" '{print \$1"\t"\$2"\t"\$3}' > regions.bed   

            bcftools merge ${vcf_files.join(' ')} -0 -R regions.bed | \\
                bcftools norm -D -d both -f ${ref_genome} | \\
                bcftools view -o ${out} -O z -G
            bcftools index ${out}
            """
        }
}
