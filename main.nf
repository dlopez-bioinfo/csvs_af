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
            def (vcf_path, bed_path) = line.trim().tokenize('\t')
            def sample_id = vcf_path.tokenize('/').last().replaceFirst(/\.vcf\.gz$/, '')
            def bed_name = bed_path == 'genome' ? 'genome' : bed_path.tokenize('/')[-1]
            tuple(sample_id, vcf_path, bed_name, bed_path)
        }
        .dump(tag: 'SAMPLE_BED')
        .set { ch_sample_bed }


    // set bed channel
    ch_sample_bed
        .filter { id, vcf, bed_name, bed_path -> bed_name.trim() != params.genome_sample_str }  
        .map { id, vcf, bed_name, bed_path -> bed_path.trim() }
        .distinct()
        .collect()
        .dump(tag: 'BED_LIST')
        .set { ch_bed_list }

    // create intervals    
    MAKE_INTERVALS(
        ch_bed_list,
        params.ref_genome_sizes)

    // normalize VCFs
    NORMALIZE_VCF(
        ch_sample_bed,
        params.ref_genome)

    // Create list of VCFs per interval
    MAKE_INTERVALS.out
        .splitCsv(header: true, sep: '\t')
        .flatMap { row ->
            def region_id = "${row.chrom}_${row.start}_${row.end}"
            def bed_list = row.list.split(',').collect { it.trim() }
            bed_list.collect { bed -> tuple(bed, region_id) }
        }
        .dump(tag:'INTERVAL')
        .combine(NORMALIZE_VCF.out, by: 0)
        .groupTuple(by: 1)
        .map { bed_names, interval, bed_paths, ids, vcfs, idx -> tuple(interval, vcfs, idx)} 
        .dump(tag:'GROUPED')
        .set { ch_interval_vcfs }

    MERGE_INTERVAL(
        ch_interval_vcfs,        
        params.ref_genome
        )

    CONCAT(
        MERGE_INTERVAL.out.collect()
        )
}

