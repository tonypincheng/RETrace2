#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Process: Bootstrap phylogenetic trees
process bootstrap_trees {
    publishDir "${params.output_dir}/bootstrap", mode: 'copy'
    
    input:
    path microsatellites
    path tree
    
    output:
    path "bootstrap_trees.nwk", emit: bootstrap_trees
    path "bootstrap_support.pdf", emit: support_plot
    
    script:
    """
    Rscript $baseDir/scripts/bootstrap/bootstrap_trees.R \
      --input "${microsatellites}" \
      --tree "${tree}" \
      --output "bootstrap_trees.nwk" \
      --plot "bootstrap_support.pdf" \
      --replicates ${params.bootstrap_replicates}
    """
}

// Module workflow
workflow BOOTSTRAP {
    take:
    microsatellites
    tree
    
    main:
    bootstrap_trees(microsatellites, tree)
    
    emit:
    bootstrap_trees = bootstrap_trees.out.bootstrap_trees
    support_plot = bootstrap_trees.out.support_plot
} 