#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl = 2

process CONCAT {
    container "${params.container__bcftools}"
    label 'high_proc'

    publishDir "${params.output_folder}/", mode: 'copy', overwrite: true

    input:
        path(vcf_list)
        val(sample_md5)
        
    output:
        tuple path("csvs_*.vcf.gz"), path("csvs_*.vcf.gz.csi")

    script:
        def out = "csvs_${sample_md5}.vcf.gz"

        """
        bcftools merge  --threads ${task.cpus} *_merged.vcf.gz | \\
          bcftools +fill-tags -- -t "AC,AN,AF" | \\ 
          bcftools sort -o ${out} -Oz
        bcftools index ${out} --threads ${task.cpus}
        """
}
