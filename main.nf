#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { MAPPING } from './modules/mapping/main.nf'
include { HIPSTR } from './modules/hipstr/main.nf'
include { PHYLO } from './modules/phylo/main.nf'
include { BOOTSTRAP } from './modules/bootstrap/main.nf' optional params.run_bootstrap
include { EVALUATION } from './modules/evaluation/main.nf' optional params.run_evaluation
include { METHYLATION } from './modules/methylation/main.nf' optional params.run_methylation

// Pipeline parameters
params.input_dir = "data/MSH2"
params.output_dir = "results/"
params.fastq_pattern = "*.fastq.gz"
params.threads = 30
params.memory = '100.GB'
params.help = false

// Reference genome parameters
params.genomes_base = "/path/to/reference/genomes"
params.genome = "mm39"
params.download_reference = false
params.bwa_index_path = null

// Optional analysis parameters
params.run_bootstrap = false
params.run_evaluation = false
params.run_methylation = false
params.ground_truth = null
params.bootstrap_iterations = 100

// If methylation is enabled, set default methylation input parameters
if (params.run_methylation) {
    params.methylation_input_dir = "data/methylation/"
    params.methylation_fastq_pattern = "*.fastq.gz"
}

// HipSTR parameters
params.min_qual = 0.9
params.min_reads = 10
params.max_stutter = 1.0
params.target_bed = "${baseDir}/resources/${params.genome}_microsatellites.bed"

// Help message
def helpMessage() {
    def methylationHelp = params.run_methylation ? """
      --methylation_input_dir Directory containing methylation FASTQ files (default: ${params.methylation_input_dir})
      --methylation_fastq_pattern Pattern to match methylation FASTQ files (default: ${params.methylation_fastq_pattern})
    """ : ""
    
    log.info"""
    ===========================================
      RETrace2 Pipeline v1.0
    ===========================================
    
    Usage:
      nextflow run main.nf --input_dir /path/to/fastqs --output_dir /path/to/results
    
    Mandatory arguments:
      --input_dir      Directory containing input FASTQ files (default: ${params.input_dir})
    
    Optional arguments:
      --output_dir      Directory for output files (default: ${params.output_dir})
      --fastq_pattern   Pattern to match FASTQ files (default: ${params.fastq_pattern})
      --threads         Number of CPU threads to use (default: ${params.threads})
      --memory          Memory to allocate for processes (default: ${params.memory})
      
      --genome          Reference genome: 'mm39' or 'hg38' (default: ${params.genome})
      --genomes_base    Directory containing reference genomes (default: ${params.genomes_base})
      --bwa_index_path  Path to BWA index [optional]. If not specified, will use ${params.genomes_base}/${params.genome}/bwa-index/${params.genome}.fa
      --download_reference  Download reference genome if not available (default: ${params.download_reference})
      --target_bed      BED file with microsatellite regions (default: ${params.target_bed})
      
      --min_qual        Minimum quality score for HipSTR (default: ${params.min_qual})
      --min_reads       Minimum number of reads for HipSTR (default: ${params.min_reads})
      --max_stutter     Maximum stutter ratio for HipSTR (default: ${params.max_stutter})
      
      --run_bootstrap   Run bootstrap analysis (default: ${params.run_bootstrap})
      --bootstrap_iterations Number of bootstrap iterations (default: ${params.bootstrap_iterations})
      
      --run_evaluation  Run evaluation (requires ground truth) (default: ${params.run_evaluation})
      --ground_truth    Path to ground truth data (default: ${params.ground_truth})
      
      --run_methylation Run methylation analysis (default: ${params.run_methylation})
      --methylation_input_dir Directory containing methylation FASTQ files (default: ${params.methylation_input_dir})
      --methylation_fastq_pattern Pattern to match methylation FASTQ files (default: ${params.methylation_fastq_pattern})
      --help            Display this help message
    """
}

// Show help message if --help specified
if (params.help) {
    helpMessage()
    exit 0
}

// Check for required parameters
if (!params.input_dir) {
    log.error "Input directory not specified!"
    exit 1
}

// Log pipeline info
log.info"""
===========================================
 RETrace2 Pipeline v1.0
===========================================
Input directory    : ${params.input_dir}
Output directory   : ${params.output_dir}
Threads            : ${params.threads}
Memory             : ${params.memory}
Reference genome   : ${params.genome}
HipSTR parameters  :
  - Min quality    : ${params.min_qual}
  - Min reads      : ${params.min_reads}
  - Max stutter    : ${params.max_stutter}
Optional analyses  :
  - Bootstrap      : ${params.run_bootstrap}
  - Evaluation     : ${params.run_evaluation}
  - Methylation    : ${params.run_methylation}
===========================================
"""

// Input channel for FASTQ files
input_ch = Channel.fromPath("${params.input_dir}/${params.fastq_pattern}", checkIfExists: true)
                   .map { file -> tuple(file.simpleName, file) }

// Input channel for methylation FASTQ files (only if methylation is enabled)
if (params.run_methylation) {
    methylation_input_ch = Channel.fromPath("${params.methylation_input_dir}/${params.methylation_fastq_pattern}", checkIfExists: true)
                                 .map { file -> tuple(file.simpleName, file) }
}

// Main workflow
workflow {
    // Core pipeline
    MAPPING(input_ch)
    HIPSTR(MAPPING.out.bam)
    PHYLO(HIPSTR.out.vcf)
    
    // Capture main outputs
    tree_ch = PHYLO.out.tree
    matrix_ch = PHYLO.out.matrix
    stats_ch = PHYLO.out.stats
    
    // Optional analyses
    if (params.run_bootstrap) {
        BOOTSTRAP(tree_ch, HIPSTR.out.vcf)
        bootstrap_support_ch = BOOTSTRAP.out.support
        bootstrap_trees_ch = BOOTSTRAP.out.trees
    }
    
    if (params.run_evaluation) {
        if (params.ground_truth) {
            EVALUATION(tree_ch, file(params.ground_truth))
            evaluation_results_ch = EVALUATION.out.results
        } else {
            log.warn "Evaluation requested but no ground truth provided. Skipping evaluation."
        }
    }
    
    if (params.run_methylation) {
        METHYLATION(methylation_input_ch)
        methylation_results_ch = METHYLATION.out.results
        methylation_stats_ch = METHYLATION.out.detailed_stats
        methylation_summary_ch = METHYLATION.out.summary_stats
        methylation_plot_ch = METHYLATION.out.summary_plot
    }
    
    // Print workflow completion message
    tree_ch.view { "Phylogenetic tree completed: ${it}" }
}

// Handle workflow completion
workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'SUCCESS' : 'FAILED' }"
    log.info "Execution duration: $workflow.duration"
    
    if (workflow.success) {
        log.info """
        ===========================================
         RETrace2 Pipeline - COMPLETED SUCCESSFULLY
        ===========================================
        Results are available in: ${params.output_dir}
        
        Core results:
        - Phylogenetic tree: ${params.output_dir}/phylo/phylogenetic_tree.nwk
        - Distance matrix: ${params.output_dir}/phylo/distance_matrix.txt
        """
        
        if (params.run_bootstrap) {
            log.info """
        Bootstrap results:
        - Support values: ${params.output_dir}/bootstrap/bootstrap_support.txt
        - Bootstrap trees: ${params.output_dir}/bootstrap/bootstrap_trees.nwk
            """
        }
        
        if (params.run_evaluation && params.ground_truth) {
            log.info """
        Evaluation results:
        - Evaluation metrics: ${params.output_dir}/evaluation/evaluation_results.txt
            """
        }
        
        if (params.run_methylation) {
            log.info """
        Methylation results:
        - Methylation profiles: ${params.output_dir}/methylation/
        - Methylation stats: ${params.output_dir}/methylation/methylation_stats.txt
        - Methylation summary: ${params.output_dir}/methylation/methylation_summary.txt
        - Methylation plot: ${params.output_dir}/methylation/methylation_plot.png
            """
        }
        
        log.info "==========================================="
    }
}
