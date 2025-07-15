#!/usr/bin/env nextflow

// Using DSL-2
//nextflow.enable.dsl = 2

process MERGE_INTERVAL {
    container "${params.container__bcftools}"
    label 'high_proc'


    input:
        tuple val(raw_intervals), path(vcf_list), path(vcf_index)
        path ref_genome
        path ref_genome_fai

        
    output:
        tuple path("*_merged.vcf.gz"), path("*_merged.vcf.gz.csi")

    script:        
        def vcf_files = vcf_list.collect { it.getName() }
        def out = raw_intervals.md5() + "_merged.vcf.gz"
        def half_cpus = (task.cpus / 2).toInteger()
        def quarter_cpus = (half_cpus / 2).toInteger()

        if (vcf_files.size() == 1) {
            """
            echo "${raw_intervals}" | tr ',' '\n' | awk -F"_" '{print \$1"\t"\$2"\t"\$3}' > regions.bed   

            bcftools view ${vcf_files[0]} -o ${out} -O z -G -R regions.bed --threads ${task.cpus}
            bcftools index ${out} --threads ${task.cpus}
            """
        } else {
            """
            echo "${raw_intervals}" | tr ',' '\n' | awk -F"_" '{print \$1"\t"\$2"\t"\$3}' > regions.bed
            chr=\$(head -n 1 regions.bed | cut -f 1)   

            bcftools merge *_norm.vcf.gz -0 -r \${chr} --threads ${quarter_cpus} | \\
                bcftools norm -d exact -d both -f ${ref_genome} -T regions.bed --threads ${quarter_cpus} | \\
                bcftools +fill-tags -- -t "AC,AN,AF" | \\
                bcftools view -o ${out} -O z -G --threads ${half_cpus}
            bcftools index ${out} --threads ${task.cpus}
            """
        }
}
