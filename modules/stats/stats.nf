#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process COUNT_TARGETS {
    publishDir "${params.output_dir}/stats/ms_counts", mode: 'copy'
    container "quay.io/biocontainers/pysam:0.22.1--py39hdd5828d_3"
    conda "bioconda::pysam=0.22.1"

    input:
    tuple val(sample_id), path(bam_reads), path(bam_index)

    output:
    path("${sample_id}_ms_counts.txt"), emit: ms_counts
    
    script:
    """
    python ${baseDir}/modules/stats/count_targets.py \
                --bam ${bam_reads} \
                --target_bed ${params.target_bed} \
                --output ${sample_id}_ms_counts.txt
    """
}

process COMBINE_STATS {
    publishDir "${params.output_dir}/stats", mode: 'copy'
    //container "quay.io/biocontainers/python:3.9"
    //conda "conda-forge::python=3.9 conda-forge::matplotlib=3.9.4 seaborn=0.13.2"
    
    input:
    path(ms_counts)
    path(allc_files)
    
    output:
    path("sample_stats.tsv"), emit: sample_stats
    path("summary_stats.txt"), emit: summary_stats
    path("figures"), emit: figures

    script:
    // If no methylation data is provided, pass an empty list to the script
    def allc_files_param = allc_files ? "--allc-files ${allc_files}" : ""
    """
    python ${baseDir}/modules/stats/generate_summary.py \
        --ms-counts ${ms_counts} \
        ${allc_files_param} \
        --min-reads-per-target ${params.min_reads} \
        --min-targets-per-sample ${params.min_targets} \
        --min-cpgs-per-sample ${params.min_cpgs} \
        --output-dir .
    """
}

workflow STATS {
    take:
    bam
    allc
    
    main:
    // Count microsatellites targets
    COUNT_TARGETS(bam)
    
    // Prepare allc channel for stats generation
    allc_files_ch = Channel.empty()
    if (params.run_methylation) {
        allc_files_ch = allc.map { it -> it[1] }.collect()
    }
    
    // Generate statistics incorporating both MS counts and methylation data if available
    // Use ifEmpty([]) to ensure the process runs even when allc_files_ch is empty
    COMBINE_STATS(COUNT_TARGETS.out.ms_counts.collect(), allc_files_ch.ifEmpty([]))

    emit:
    sample_stats = COMBINE_STATS.out.sample_stats
}