#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl = 2

process NORMALIZE_VCF {
    container "${params.container__bcftools}"
    label 'high_proc'
    tag "$id"

    input:
        tuple val(id), path(vcf), val(bed_name), val(bed_path)
        path ref_genome
        
    output:
        tuple val(bed_name), val(bed_path), val(id), path("*_norm.vcf.gz"), path("*_norm.vcf.gz.csi")

    script:        
        def half_cpus = (task.cpus / 2).toInteger()
        def out = "${id}_norm.vcf.gz"
        """
            bcftools norm --threads ${half_cpus} -D -d both -f ${ref_genome} ${vcf} | \\
                bcftools view --exclude-uncalled --threads ${half_cpus} -o ${out} -O z 
            bcftools index --threads ${half_cpus} ${out}
        """
}
