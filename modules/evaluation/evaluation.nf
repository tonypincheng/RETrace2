#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process EVALUATE_TREE {
    publishDir "${params.output_dir}/evaluation", mode: 'copy'
    container "tonypincheng/retrace2-python:latest"
    
    input:
    path(tree)
    path(ground_truth)
    
    output:
    path("evaluation_results.txt"), emit: results
    path("evaluation_plot.pdf"), emit: plot
    
    script:
    """
    # Evaluate tree accuracy
    python $baseDir/scripts/evaluation/evaluate.py \
        --tree $tree \
        --ground_truth $ground_truth \
        --output evaluation_results.txt \
        --plot evaluation_plot.pdf
    """
}

workflow EVALUATION {
    take:
    tree
    ground_truth
    
    main:
    // Evaluate tree accuracy
    eval_results = EVALUATE_TREE(tree, ground_truth)
    
    emit:
    results = eval_results.results
    plot = eval_results.plot
} 