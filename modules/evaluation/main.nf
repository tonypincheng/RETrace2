#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Process: Evaluate tree accuracy
process evaluate_tree {
    publishDir "${params.output_dir}/evaluation", mode: 'copy'
    
    input:
    path tree
    path ground_truth
    
    output:
    path "tree_accuracy.txt", emit: accuracy
    path "tree_comparison.pdf", emit: comparison_plot
    
    script:
    """
    python $baseDir/scripts/evaluation/evaluate_tree_accuracy.py \
      --tree "${tree}" \
      --ground_truth "${ground_truth}" \
      --output "tree_accuracy.txt" \
      --plot "tree_comparison.pdf"
    """
}

// Module workflow
workflow EVALUATION {
    take:
    tree
    ground_truth
    
    main:
    evaluate_tree(tree, ground_truth)
    
    emit:
    accuracy = evaluate_tree.out.accuracy
    comparison_plot = evaluate_tree.out.comparison_plot
} 