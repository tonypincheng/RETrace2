#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Process: Differential Methylation Analysis
process differential_methylation {
    publishDir "${params.output_dir}/differential", mode: 'copy'
    
    input:
    path methylation_beds
    
    output:
    path "differential_methylation_results.txt", emit: diff_results
    path "differential_methylation_plots.pdf", emit: diff_plots
    
    script:
    """
    Rscript $baseDir/scripts/differential/diff_meth.R \
      --input "${methylation_beds}" \
      --output "differential_methylation_results.txt" \
      --plots "differential_methylation_plots.pdf"
    """
}

// Module workflow
workflow DIFFERENTIAL {
    take:
    methylation_bed
    
    main:
    differential_methylation(methylation_bed.collect())
    
    emit:
    diff_results = differential_methylation.out.diff_results
    diff_plots = differential_methylation.out.diff_plots
} 