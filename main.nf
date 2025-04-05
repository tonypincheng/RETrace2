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
params.input_dir = "data/"
params.output_dir = "results/"
params.fastq_pattern = "*.fastq.gz"
params.threads = 30
params.memory = '100.GB'
params.help = false

// Reference genome parameters
params.genomes_base = "/path/to/reference/genomes"
params.genome = "mm39"
params.download_reference = false

// Optional analysis parameters
params.run_bootstrap = false
params.run_evaluation = false
params.run_methylation = false
params.ground_truth = null
params.bootstrap_iterations = 100

// HipSTR parameters
params.min_qual = 0.9
params.min_reads = 10
params.max_stutter = 1.0

// Help message
def helpMessage() {
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
      --download_reference  Download reference genome if not available (default: ${params.download_reference})
      
      --min_qual        Minimum quality score for HipSTR (default: ${params.min_qual})
      --min_reads       Minimum number of reads for HipSTR (default: ${params.min_reads})
      --max_stutter     Maximum stutter ratio for HipSTR (default: ${params.max_stutter})
      
      --run_bootstrap   Run bootstrap analysis (default: ${params.run_bootstrap})
      --bootstrap_iterations Number of bootstrap iterations (default: ${params.bootstrap_iterations})
      
      --run_evaluation  Run evaluation (requires ground truth) (default: ${params.run_evaluation})
      --ground_truth    Path to ground truth data (default: ${params.ground_truth})
      
      --run_methylation Run methylation analysis (default: ${params.run_methylation})
      
      --help            Display this help message
    """
}

// Show help message if --help specified
if (params.help) {
    helpMessage()
    exit 0
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

// Main workflow
workflow {
    // Run core pipeline
    MAPPING(input_ch)
    HIPSTR(MAPPING.out.bam)
    PHYLO(HIPSTR.out.vcf)
    
    // Optional analyses
    if (params.run_bootstrap) {
        BOOTSTRAP(PHYLO.out.tree, params.bootstrap_iterations)
    }
    
    if (params.run_evaluation) {
        EVALUATION(PHYLO.out.tree, params.ground_truth)
    }
    
    if (params.run_methylation) {
        METHYLATION(input_ch)
    }
    
    // Print workflow completion message
    PHYLO.out.tree.view { "Phylogenetic tree completed: ${it}" }
}

// Handle workflow completion
workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'SUCCESS' : 'FAILED' }"
    log.info "Execution duration: $workflow.duration"
}
