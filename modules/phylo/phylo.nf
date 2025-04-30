#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process BUILD_TREE {
    publishDir "${params.output_dir}/phylo", mode: 'copy'
    
    input:
    path(alleleDict)
    path(sample_list)
    
    output:
    path("*.buildPhylo.newick-original.txt"), emit: newick_tree
    path("*.buildPhylo.stats.txt"), optional: true, emit: stats
    path("*.buildPhylo.distDict.pkl"), emit: dist_dict
    path("*.buildPhylo.newick-bootstrap.txt"), optional: true, emit: bootstrap_tree
    path("*_bootstrap_stats.txt"), optional: true, emit: bootstrap_stats
    
    script:
    """
    python ${baseDir}/modules/phylo/build_phylo.py \\
        --alleleDict ${alleleDict} \\
        --sample_list ${sample_list} \\
        --prefix ${params.output_prefix} \\
        --dist_metric ${params.dist_metric} \\
        --outgroup ${params.outgroup} \\
        ${params.run_bootstrap ? "--bootstrap --bootstrap_iterations ${params.bootstrap_iterations}" : ""}
    """
}

process VIEW_PHYLO {
    publishDir "${params.output_dir}/phylo", mode: 'copy'
    
    input:
    path(samplesheet)
    path(tree_file)
    
    output:
    path("*.viewPhylo.pdf"), emit: tree_pdf
    path("*.viewPhylo.png"), emit: tree_png
    
    script:
    """
    python ${baseDir}/modules/phylo/view_phylo.py \\
        --samplesheet ${samplesheet} \\
        --tree_file ${tree_file} \\
        --prefix ${params.output_prefix} \\
        ${params.run_bootstrap ? "--bootstrap" : ""}
    """
}

workflow PHYLO {
    take:
    alleleDict
    sample_list
    
    main:
    BUILD_TREE(alleleDict, sample_list)
    
    // Choose tree file based on bootstrap parameter
    tree_file = params.run_bootstrap ? BUILD_TREE.out.bootstrap_tree : BUILD_TREE.out.newick_tree
    
    // Create a channel from the samplesheet path
    samplesheet_ch = Channel.fromPath(params.samplesheet)
    
    VIEW_PHYLO(samplesheet_ch, tree_file)
    
    emit:
    newick_tree = BUILD_TREE.out.newick_tree
    stats = BUILD_TREE.out.stats
    dist_dict = BUILD_TREE.out.dist_dict
    bootstrap_tree = BUILD_TREE.out.bootstrap_tree
    bootstrap_stats = BUILD_TREE.out.bootstrap_stats
    tree_pdf = VIEW_PHYLO.out.tree_pdf
    tree_png = VIEW_PHYLO.out.tree_png
} 