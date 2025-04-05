#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Process: Region-based Analysis
process region_analysis {
    publishDir "${params.output_dir}/region", mode: 'copy'
    
    input:
    path methylation_beds
    
    output:
    path "region_analysis_results.txt", emit: region_results
    path "region_analysis_summary.pdf", emit: region_summary
    
    script:
    """
    python $baseDir/scripts/region/region_analysis.py \
      --input "${methylation_beds}" \
      --output "region_analysis_results.txt" \
      --summary "region_analysis_summary.pdf"
    """
}

// Module workflow
workflow REGION {
    take:
    methylation_bed
    
    main:
    region_analysis(methylation_bed.collect())
    
    emit:
    region_results = region_analysis.out.region_results
    region_summary = region_analysis.out.region_summary
} 