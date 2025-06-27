#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl = 2

process MAKE_INTERVALS {
    container "${params.container__bedtools}"
    label 'single_proc'

    input:
        path bed_info_list        
        path sizes
        
    output:
        path "intervals.txt"

    script:
        def bed_names = bed_info_list.collect { it }.join(' ')
        def bed_paths = bed_info_list.collect { it.getName() }.join(' ')
        """
            bedtools multiinter \\
                -header \\
                -empty \\
                -g ${sizes} \\
                -names ${bed_names} \\
                -i ${bed_paths} > aux.txt

            awk -F "\t" '\$5!~/none/{if(NR>1)\$5="none,"\$5} {print}' OFS="\t" aux.txt > intervals.txt 

            sed -i -r 's/none/${params.genome_sample_str}/g' intervals.txt
        """
}
