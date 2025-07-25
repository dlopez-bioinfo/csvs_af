#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl = 2

process CONCAT {
    container "${params.container__bcftools}"
    label 'med_proc'

    publishDir "${params.output_folder}/", mode: 'copy', overwrite: true

    input:        
        path(vcf_list)
        path(sample_bed_file)
        
    output:
        tuple path("csvs_*_*/csvs_*_all.vcf.gz"), path("csvs_*_*/csvs_*_all.vcf.gz.csi")

    script:        
        """
        m=\$(md5sum ${sample_bed_file} |cut -f 1 -d ' ')
        OUT_DIR="csvs_${params.ref_genome_id}_\${m::6}/"
        OUT="\${OUT_DIR}/csvs_\${m::6}_all.vcf.gz"

        mkdir -p "\${OUT_DIR}"

        bcftools concat *.vcf.gz | \\
            bcftools sort | \\
            bcftools view -o \${OUT} -Oz --threads ${task.cpus}
        bcftools index \${OUT} --threads ${task.cpus}
        """
}
