#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Process: Visualization
process visualization {
    publishDir "${params.output_dir}/visualization", mode: 'copy'
    
    input:
    path methylation_beds
    
    output:
    path "methylation_plots.pdf", emit: plots
    path "methylation_summary.html", emit: summary
    
    script:
    """
    Rscript $baseDir/scripts/visualization/plot_methylation.R \
      --input "${methylation_beds}" \
      --plots "methylation_plots.pdf" \
      --summary "methylation_summary.html"
    """
}

// Module workflow
workflow VISUALIZATION {
    take:
    methylation_bed
    
    main:
    visualization(methylation_bed.collect())
    
    emit:
    plots = visualization.out.plots
    summary = visualization.out.summary
} 