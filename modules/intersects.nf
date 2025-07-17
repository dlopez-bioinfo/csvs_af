#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl = 2

process INTERSECTS {
    container "${params.container__bedtools}"
    label 'single_proc'

    input:
        tuple val(bed_names_input), path(bed_paths_input)
        path sizes
        
    output:
        path "intervals.txt"

    script:
        def bed_names = bed_names_input.collect { it }.join(' ')
        def bed_paths = bed_paths_input.collect { it.getName() }.join(' ')
        """
            # create dummy genome region file
            awk -F"\t" '{print \$1"\t0\t"\$2}' sizes.hs37d5  > genome.bed

            bedtools multiinter \\
                -header \\
                -empty \\
                -g ${sizes} \\
                -names ${params.genome_sample_str} ${bed_names} \\
                -i genome.bed ${bed_paths} > intervals.txt
        """
}
