#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Default parameter values (will be overridden by main workflow values)
params.output_dir = params.output_dir ?: "results/"
params.target_bed = params.target_bed ?: "default_targets.bed"
params.genomes_base = params.genomes_base ?: "/path/to/reference/genomes"
params.genome = params.genome ?: "mm39"
params.bwa_index_path = params.bwa_index_path ?: null

process count_targets {
    publishDir "${params.output_dir}/stats/ms_counts", mode: 'copy'
    //container "quay.io/biocontainers/bedtools:2.30.0--h94655ef_1"
    //conda "bioconda::"

    input:
    tuple val(sample_id), path(bam_reads), path(bam_index)

    output:
    path("${sample_id}_ms_counts.txt")
    
    script:
    println("${projectDir}")
    """
    python count_targets.py \
                --bam ${bam_reads} \
                --target_bed ${params.target_bed} \
                --min_depth 30 \
                --output ${sample_id}_ms_counts.txt
    """
}


workflow STATS {
    take:
    bam
    
    main:
    // Count microsatellites targets
    count_targets(bam)
    
    emit:
    ms_counts = count_targets.out
}