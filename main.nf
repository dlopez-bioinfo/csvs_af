#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl = 2

/*
================================================================
    MODULES
================================================================
*/
//TODO: include modules
include { MAKE_INTERVALS } from './modules/make_intervals'
include { NORMALIZE_VCF } from './modules/normalize_vcf'
include { MERGE_INTERVAL } from './modules/merge_interval'
include { CONCAT } from './modules/concat'
include { SORT_BED } from './modules/sort_bed'
include { PROCESS_INTERVALS } from './modules/process_intervals'
include { REDUCE_INTERVALS } from './modules/reduce_intervals'

/*
================================================================
    SUBWORKFLOWS
================================================================
*/
//TODO: include subworkflows

/*
================================================================
MAIN WORKFLOW
================================================================
*/
workflow {
    //**********
    //*  INIT  *
    //**********

    // Print input parameters
    log.info('Input parameters:')
    params.each { key, value ->
        if (value) {
            println("${key}: ${value}")
        }
    }
    log.info('--------------------------\n')


    // set main sample-bed channel
    Channel.fromPath(params.sample_bed_file)
        .splitText()
        .map { line -> 
            def (bed_path, vcf_path, gender) = line.trim().tokenize('\t')
            def sample_id = vcf_path.tokenize('/').last().replaceFirst(/\.vcf\.gz$/, '')
            def bed_name = bed_path.trim() == params.genome_sample_str ? bed_path : bed_path.md5()
            tuple(sample_id, vcf_path, bed_name, bed_path, gender)
        }
        //.dump(tag: 'SAMPLE_BED')
        .set { ch_sample_bed }

 
    // set bed channel
        ch_sample_bed
        .filter { id, vcf, bed_name, bed_path, gender -> bed_path.trim() != params.genome_sample_str }  
        .map { id, vcf, bed_name, bed_path, gender -> tuple(bed_name, bed_path) }
        .unique()
        //.dump(tag: 'BED_LIST')        
        .set { ch_bed_list }

     // Sort bed and rename to unique names
    SORT_BED(
        ch_bed_list
        )

    SORT_BED.out        
        .map { id, path -> tuple(id, file(path)) }        
        .collect(flat:false)         
        .map { list_of_tuples ->
            // Ordenar por id (primer elemento del tuple)
            def sorted = list_of_tuples.sort { a, b -> a[0] <=> b[0] } // sort by bed_name to allow caching

            def names = sorted*.getAt(0)
            def paths = sorted*.getAt(1)
            return tuple(names, paths)
        }
        //.dump(tag: 'SORTED_BED_LIST2')        
        .set { ch_sorted_bed_list }

    // create intervals    
    MAKE_INTERVALS(
        ch_sorted_bed_list, //[bed_name, bed_path]
        params.ref_genome_sizes) | PROCESS_INTERVALS


    // normalize VCFs
    NORMALIZE_VCF(
        ch_sample_bed,
        params.ref_genome,
        params.ref_genome_fai)

    
    // Create interval files
    PROCESS_INTERVALS.out
        .flatten()
        .splitCsv(header: true, sep: '\t')
        .map { row -> return tuple(row.list, row.id)}  // [bed_name, interval_list]
        .set { ch_intervals }


    // Reduce intervals
    REDUCE_INTERVALS(
        ch_intervals
    )     

    // Merge intervals
    REDUCE_INTERVALS.out
        .combine(NORMALIZE_VCF.out, by: 0)  //[ bed_name, bed_path, id, vcf, vcf_index, gender_file]        
        .groupTuple(by: 1)
        .map { bed_name, region_file, bed_path, sample_id, vcf_list, vcf_idx, gender_file_list->
            return tuple(region_file, vcf_list, vcf_idx, gender_file_list)
        }
        //.dump(tag: 'INTERVAL_VCF_LIST') // [bed_name, interval_list, bed_path, sample_id, vcf_list, vcf_idx, gender_file_list]
        .set { ch_interval_vcfs }

    
    MERGE_INTERVAL(
        ch_interval_vcfs,
        params.ref_genome,
        params.ref_genome_fai)


    CONCAT(
        MERGE_INTERVAL.out.collect(),
        file(params.sample_bed_file)
        )
}

