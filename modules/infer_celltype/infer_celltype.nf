#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process CALCULATE_PD_MATRIX {
    publishDir "${params.output_dir}/infer_celltype/pd_matrix", mode: 'copy'

    input:
    path allc_files
    path celltype_ref_files

    output:
    path "pairwise_dissimilarity_matrix.csv", emit: pd_matrix
    path "shared_sites_matrix.csv", emit: sites_matrix

    script:
    """
    python ${baseDir}/modules/infer_celltype/calculate_pd_matrix_batched.py \
        --sc_files ${allc_files} \
        --ref_files ${celltype_ref_files} \
        --output_dir . \
        --min_reads ${params.min_reads_per_site} \
        --min_sites ${params.min_shared_sites} \
        --n_processes ${task.cpus} \
        ${params.all_cytosines ? '--all_cytosines' : ''}
    """
}

// process ASSIGN_CELLTYPE {
//     publishDir "${params.output_dir}/infer_celltype/assignments", mode: 'copy'

//     input:
//     path pd_matrix
//     path sites_matrix

//     output:
//     path "celltype_assignments.csv", emit: assignments
//     path "celltype_confidence_scores.csv", emit: confidence_scores

//     script:
//     """
//     python ${baseDir}/modules/infer_celltype/assign_celltype.py \
//         --pd_matrix ${pd_matrix} \
//         --sites_matrix ${sites_matrix} \
//         --output_dir . \
//         --min_confidence_score ${params.min_confidence_score ?: 0.7} \
//         --assignment_method ${params.celltype_assignment_method ?: 'nearest_neighbor'}
//     """
// }

// Create channel for reference files
if (params.celltype_ref_dir) {
    // Specify directory + pattern
    celltype_ref_path = "${params.celltype_ref_dir}/${params.celltype_ref_pattern}"
    celltype_ref_ch = Channel.fromPath(celltype_ref_path)
} else {
    // No reference files specified
    log.info "No celltype reference directory specified (celltype_ref_dir). Cell type inference will be skipped."
    celltype_ref_ch = Channel.empty()
}

workflow INFER_CELLTYPE {
    take:
    allc_files
    
    main:
    // Calculate pairwise dissimilarity matrix
    CALCULATE_PD_MATRIX(allc_files.collect(), celltype_ref_ch.collect())
    
    // Assign cell types based on the matrix
    //ASSIGN_CELLTYPE(CALCULATE_PD_MATRIX.out.pd_matrix)
    
    emit:
    pd_matrix = CALCULATE_PD_MATRIX.out.pd_matrix
    //assignments = ASSIGN_CELLTYPE.out.assignments

} 