#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process count_targets {
    publishDir "${params.output_dir}/stats/ms_counts", mode: 'copy'


    input:
    tuple val(sample_id), path(bam_reads), path(bam_index)

    output:
    
    
    script:
    
    """
    python count_targets.py ${sample_id} ${bam_reads}

    python count_targets.py \
                --bam ${bam_reads} \
                --target_bed microsatellites.bed \
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