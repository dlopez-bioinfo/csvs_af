#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl = 2

process CONCAT {
    container "${params.container__bcftools}"
    label 'high_proc'

    publishDir "${params.output_folder}/", mode: 'copy', overwrite: true

    input:
        path(vcf_list)
        path(sample_bed_file)
        
    output:
        tuple path("csvs_*.vcf.gz"), path("csvs_*.vcf.gz.csi")

    script:
        """
        OUT="csvs_"\$(md5sum ${sample_bed_file} |cut -f 1 -d ' ')".vcf.gz"

        bcftools merge  --threads ${task.cpus} *_merged.vcf.gz | \\
          bcftools +fill-tags -- -t "AC,AN,AF" | \\
          bcftools sort -o \${OUT} -Oz
        bcftools index \${OUT} --threads ${task.cpus}
        """
}
