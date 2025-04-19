#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process bootstrap {
    tag "bootstrap"
    publishDir "${params.output_dir}/bootstrap", mode: 'copy'
    
    input:
    path(tree)
    path(vcf_files)
    
    output:
    path("bootstrap_trees.nwk"), emit: trees
    path("bootstrap_support.txt"), emit: support
    path("bootstrap_plot.pdf"), emit: plot
    
    script:
    """
    # Run bootstrap analysis
    python $baseDir/scripts/bootstrap/bootstrap.py \
        --tree $tree \
        --vcf $vcf_files \
        --iterations ${params.bootstrap_iterations} \
        --output bootstrap_trees.nwk \
        --support bootstrap_support.txt \
        --plot bootstrap_plot.pdf
    """
}

workflow BOOTSTRAP {
    take:
    tree
    vcf_files
    
    main:
    // Run bootstrap analysis
    bootstrap_results = bootstrap(tree, vcf_files)
    
    emit:
    trees = bootstrap_results.trees
    support = bootstrap_results.support
    plot = bootstrap_results.plot
} 