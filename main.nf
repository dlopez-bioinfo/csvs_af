#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl = 2

/*
================================================================
    MODULES
================================================================
*/
//TODO: include modules
include { INTERSECTS } from './modules/intersects'
include { NORMALIZE_VCF } from './modules/normalize_vcf'
include { MERGE_INTERVAL } from './modules/merge_interval'
include { CONCAT } from './modules/concat'
include { SORT_BED } from './modules/sort_bed'
include { PROCESS_INTERSECTS } from './modules/process_intersects'
include { PREPARE_INTERSECTIONS } from './modules/prepare_intersections'
include { PREPARE_GENDER_FILE } from './modules/prepare_gender_file'

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
        .splitCsv(sep: "\t", strip: true)
        .map { bed_path, vcf_path, gender ->
            def sample_id = vcf_path.tokenize('/').last().replaceFirst(/\.vcf\.gz$/, '')
            def bed_hash = bed_path.trim() == params.genome_sample_str ? bed_path : bed_path.md5()[0..7]
            tuple(sample_id, vcf_path, bed_hash, bed_path, gender)
        }        
        .tap { ch_sample_bed }  // emit by-sample intermediate channel [sample_id, vcf_path, bed_hash, bed_path, gender]
        .collectFile { item ->
            [ "samples.txt", item.join(",") + '\n' ]
        }
        .set { ch_sample_bed_file }
    


    //************************
    //*   NORMALIZE VCFS   *//
    //************************
    NORMALIZE_VCF(
        ch_sample_bed,
        params.ref_genome,
        params.ref_genome_fai
    )

    PREPARE_GENDER_FILE(
        NORMALIZE_VCF.out.collect { it.getAt(5) },  // list
    )



    //******************************
    //*   CREATE INTERSECTIONS   *//
    //******************************
    // set unique bed channel
    ch_sample_bed
        .filter { sample_id, vcf_path, bed_hash, bed_path, gender -> bed_path.trim() != params.genome_sample_str }  
        .map { sample_id, vcf_path, bed_hash, bed_path, gender -> tuple(bed_hash, bed_path) }
        .unique()
        .set { ch_bed_list }

    // Sort BEDs and rename to unique names
    SORT_BED(
        ch_bed_list
    )

    // Collect sorted beds and create a list of tuples with bed_hash and bed_path
    SORT_BED.out
        .toSortedList { a, b -> b[0] <=> a[0] }.map { list_of_pairs ->
            def ids  = list_of_pairs*.getAt(0)
            def beds = list_of_pairs*.getAt(1)
            tuple(ids, beds)
        }
        .set { ch_sorted_bed_list } //[bed_hash_list, bed_path_list]

    //  Create intersections from sorted beds
    INTERSECTS(
        ch_sorted_bed_list, 
        params.ref_genome_sizes
    ) | PROCESS_INTERSECTS | splitCsv(header: false, sep: '\t') | set { ch_intersections }  // [beds_list, regions]

    

    //*****************
    //*   MERGING   *//
    //*****************    
    PREPARE_INTERSECTIONS(
        ch_intersections
            .combine(ch_sample_bed_file),
        NORMALIZE_VCF.out.collect { it.getAt(3) },  // list of VCFs
        NORMALIZE_VCF.out.collect { it.getAt(4) },  // list of VCF indices
    )

    MERGE_INTERVAL(
        PREPARE_INTERSECTIONS.out
            .combine(PREPARE_GENDER_FILE.out), // [regions.bed, vcf_dir, gender_file]
        params.ref_genome,
        params.ref_genome_fai
    )
    
    CONCAT(
        MERGE_INTERVAL.out.collect(),
        PREPARE_GENDER_FILE.out,
        file(params.sample_bed_file)
    )

}

