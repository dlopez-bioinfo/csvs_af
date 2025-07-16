#!/usr/bin/env nextflow

// Using DSL-2
//nextflow.enable.dsl = 2

process MERGE_INTERVAL {
    container "${params.container__bcftools}"
    label 'high_proc'


    input:
        tuple path(region_file), path(vcf_list), path(vcf_index), path(gender_list)
        path ref_genome
        path ref_genome_fai

        
    output:
        tuple path("*_merged.vcf.gz"), path("*_merged.vcf.gz.csi"), path("*_genderinterval.txt")

    script:                
        def vcf_files = vcf_list.collect { it.getName() }  
        def half_cpus = (task.cpus / 2).toInteger()
        def quarter_cpus = (half_cpus / 2).toInteger()

        if (vcf_files.size() == 1) {
            """
            m=\$(md5sum ${region_file} |cut -f 1 -d ' ')
            out="\${m}_merged.vcf.gz"

            gender_file="\${m}_genderinterval.txt"
            cp ${gender_list} \${gender_file}

            bcftools view ${vcf_files[0]} -o \${out} -O z -R ${region_file} --threads ${task.cpus}
            bcftools index \${out} --threads ${task.cpus}
            """
        } else {
            """
            m=\$(md5sum ${region_file} |cut -f 1 -d ' ')
            out="\${m}_merged.vcf.gz"

            #join gender information
            gender_file="\${m}_genderinterval.txt"
            cat *_gender.txt > \${gender_file}   

            bcftools merge *_norm.vcf.gz -0 --threads ${quarter_cpus} | \\
                bcftools norm -d exact -d both -f ${ref_genome} -T ${region_file} --threads ${quarter_cpus} | \\
                bcftools +fixploidy -- -s \${gender_file} | \\
                bcftools +fill-tags -- -t "AC,AN,AF" | \\
                bcftools view -o \${out} -O z --threads ${half_cpus}
            bcftools index \${out} --threads ${task.cpus}
            """
        }
}
