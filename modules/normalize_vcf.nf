#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl = 2

process NORMALIZE_VCF {
    container "${params.container__bcftools}"
    label 'high_proc'
    tag "$id"

    input:
        tuple val(id), path(vcf), val(bed_name), val(bed_path), val(gender)
        path ref_genome
        path ref_genome_fai
        
    output:
        tuple val(bed_name), val(bed_path), val(id), path("*_norm.vcf.gz"), path("*_norm.vcf.gz.csi")

    script:        
        def half_cpus = (task.cpus / 2).toInteger()
        def out = "${id}_norm.vcf.gz"
        """
            echo "chrM MT" > chr_list.txt            
            for i in {1..22} X Y
            do
                echo "chr\${i} \${i}" >> chr_list.txt
            done 

            id=\$(bcftools query -l ${vcf})

            echo "\${id} ${gender}" > gender.txt

            bcftools annotate -x INFO,^FORMAT/GT --force --rename-chrs chr_list.txt ${vcf} | \\
                bcftools norm --threads ${half_cpus} -d exact -d both -f ${ref_genome} --check-ref ws --targets \$(echo {1..22} X Y MT|sed 's/ /,/g') | \\
                bcftools +setGT -- -n . -i'FILTER!="PASS"' -tq | \\
                bcftools +fixploidy -- -s gender.txt | \\
                bcftools +fill-tags -- -t all | \\
                bcftools view --exclude-uncalled --threads ${half_cpus} -o ${out} -O z 

            bcftools index --threads ${half_cpus} ${out}
        """
}
