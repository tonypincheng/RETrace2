#!/usr/bin/env nextflow

nextflow.enable.dsl=2

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
    
}