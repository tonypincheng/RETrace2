#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process BUILD_TREE {
    tag "phylo"
    publishDir "${params.output_dir}/phylo", mode: 'copy'
    
    input:
    path(vcf_files)
    
    output:
    path("phylogenetic_tree.nwk"), emit: tree
    path("distance_matrix.txt"), emit: matrix
    path("tree_stats.txt"), emit: stats
    
    script:
    """
    # Convert VCF to distance matrix and build phylogenetic tree
    python $baseDir/scripts/phylo/build_tree.py \
        --vcf $vcf_files \
        --matrix distance_matrix.txt \
        --tree phylogenetic_tree.nwk \
        --stats tree_stats.txt
    """
}

workflow PHYLO {
    take:
    vcf_files
    
    main:
    // Build phylogenetic tree
    tree_results = BUILD_TREE(vcf_files)
    
    emit:
    tree = tree_results.tree
    matrix = tree_results.matrix
    stats = tree_results.stats
} 