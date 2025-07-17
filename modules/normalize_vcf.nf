#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl = 2

process NORMALIZE_VCF {
    container "${params.container__bcftools}"
    label 'high_proc'
    tag "$vcf_id"

    input:
        tuple val(vcf_id), path(vcf), val(bed_name), val(bed_path), val(gender)
        path ref_genome
        path ref_genome_fai
        
    output:
        tuple val(bed_name), val(bed_path), val(vcf_id), path("*_norm.vcf.gz"), path("*_norm.vcf.gz.csi"), path("*_gender.txt")

    script:        
        def half_cpus = (task.cpus / 2).toInteger()
        def out = "${vcf_id}_norm.vcf.gz"
        """
            echo "chrM MT" > chr_list.txt            
            for i in {1..22} X Y
            do
                echo "chr\${i} \${i}" >> chr_list.txt
            done 

            id=\$(bcftools query -l ${vcf})

            echo "\${id} ${gender}" > ${vcf_id}_gender.txt

            bcftools annotate -x INFO,^FORMAT/GT --force --rename-chrs chr_list.txt ${vcf} | \\
                bcftools norm --threads ${half_cpus} -d exact -d both -f ${ref_genome} --check-ref ws --targets \$(echo {1..22} X Y MT|sed 's/ /,/g') | \\
                bcftools +setGT -- -n . -i'FILTER!="PASS"' -tq | \\
                bcftools +fixploidy -- -s ${vcf_id}_gender.txt | \\
                bcftools +fill-tags -- -t "AC,AN,AF" | \\
                bcftools view --threads ${half_cpus} -o ${out} -O z 

            bcftools index --threads ${half_cpus} ${out}
        """
}
