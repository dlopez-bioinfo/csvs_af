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
            bedtools multiinter \\
                -header \\
                -empty \\
                -g ${sizes} \\
                -names ${bed_names} \\
                -i ${bed_paths} | \\
            awk -F "\t" '\$5!~/none/{if(NR>1)\$5="none,"\$5} {print}' OFS="\t" | \\
            sed -r 's/none/${params.genome_sample_str}/g' > intervals.txt
        """
}
