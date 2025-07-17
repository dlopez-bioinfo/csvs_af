#!/usr/bin/env nextflow

// Using DSL-2
//nextflow.enable.dsl = 2

process MERGE_INTERVAL {
    container "${params.container__bcftools}"
    label 'high_proc'


    input:
        tuple path(region_file), path(vcf_dir), path(gender_file)
        path ref_genome
        path ref_genome_fai

        
    output:
        tuple path("*_merged.vcf.gz"), path("*_merged.vcf.gz.csi")

    script:                        
        def half_cpus = (task.cpus / 2).toInteger()
        def quarter_cpus = (half_cpus / 2).toInteger()

        """
        m=\$(md5sum ${region_file} |cut -f 1 -d ' ')
        out="\${m}_merged.vcf.gz"

        if [[ "\$(ls ${vcf_dir}/*.vcf.gz|wc -l)" -eq 1 ]]
        then
            bcftools view ${vcf_dir}/*.vcf.gz -o \${out} -Oz -R ${region_file} --threads ${task.cpus}
            
        else
            bcftools merge ${vcf_dir}/*.vcf.gz -0 --threads ${quarter_cpus} | \\
                bcftools norm -d exact -d both -f ${ref_genome} -T ${region_file} --threads ${quarter_cpus} | \\
                bcftools +fixploidy -- -s ${gender_file} | \\
                bcftools +fill-tags -- -t "AC,AN,AF" | \\
                bcftools view -G -o \${out} -O z --threads ${half_cpus}
        fi

        bcftools index \${out} --threads ${task.cpus}
        """        
}
