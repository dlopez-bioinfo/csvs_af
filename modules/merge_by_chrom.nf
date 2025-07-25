#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl = 2

process MERGE_BY_CHROM {
    container "${params.container__bcftools}"
    label 'med_proc'

    publishDir "${params.output_folder}/", mode: 'copy', overwrite: true

    input:
        each(chr)
        path(vcf_list)
        path(gender_file)
        path(sample_bed_file)
        
    output:
        tuple path("csvs_*_*/csvs_*.vcf.gz"), path("csvs_*_*/csvs_*.vcf.gz.csi")

    script:
        """
        m=\$(md5sum ${sample_bed_file} |cut -f 1 -d ' ')
        OUT_DIR="csvs_${params.ref_genome_id}_\${m::6}/"
        OUT="\${OUT_DIR}/csvs_\${m::6}_${chr}.vcf.gz"

        mkdir -p "\${OUT_DIR}"

        bcftools merge -r ${chr} *_merged.vcf.gz | \\
            bcftools +fixploidy -- -s ${gender_file} | \\
            bcftools +fill-tags -- -t "AC,AN,AF" | \\
            bcftools sort | \\
            bcftools view -o \${OUT} -Oz --threads ${task.cpus}
        bcftools index \${OUT} --threads ${task.cpus}
        """
}
