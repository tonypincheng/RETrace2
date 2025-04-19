#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process count_targets {
    publishDir "${params.output_dir}/stats/ms_counts", mode: 'copy'
    container "quay.io/biocontainers/pysam:0.22.1--py39hdd5828d_3"
    conda "bioconda::pysam=0.22.1"

    input:
    tuple val(sample_id), path(bam_reads), path(bam_index)

    output:
    path("${sample_id}_ms_counts.txt")
    
    script:
    """
    python ${baseDir}/modules/stats/count_targets.py \
                --bam ${bam_reads} \
                --target_bed ${params.target_bed} \
                --min_depth ${params.min_reads} \
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