#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { METHYLATION } from './modules/methylation/main.nf'
include { DIFFERENTIAL } from './modules/differential/main.nf' optional params.run_differential
include { REGION } from './modules/region/main.nf' optional params.run_region_analysis
include { VISUALIZATION } from './modules/visualization/main.nf' optional params.run_visualization

// Pipeline parameters
params.input_dir = "data/"
params.output_dir = "results/"
params.fastq_pattern = "*.fastq"
params.threads = 30
params.memory = '100.GB'
params.help = false

// Reference genome parameters
params.genomes_base = "/path/to/reference/genomes"
params.genome = "mm39"
params.download_reference = false

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
      
      --run_differential    Run differential methylation analysis (default: ${params.run_differential})
      --run_region_analysis Run region-based analysis (default: ${params.run_region_analysis})
      --run_visualization   Run visualization (default: ${params.run_visualization})
      
      --help            Display this help message
    """
}

// Show help message if --help specified
if (params.help) {
    helpMessage()
    exit 0
}

// Check reference genome exists or download if needed
if (params.download_reference) {
    log.info "Checking for reference genome: ${params.genome}"
    
    // Check if reference files exist
    def refFasta = file("${params.genomes_base}/${params.genome}/genome.fa")
    
    if (!refFasta.exists()) {
        log.info "Reference genome not found. Downloading..."
        
        process downloadReference {
            output:
            path "${params.genome}", emit: ref_dir
            
            script:
            """
            bash $baseDir/bin/download_references.sh -g ${params.genome} -d ${params.genomes_base}
            """
        }
    } else {
        log.info "Reference genome found at: ${params.genomes_base}/${params.genome}"
    }
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
Optional analyses  :
  - Differential   : ${params.run_differential}
  - Region         : ${params.run_region_analysis}
  - Visualization  : ${params.run_visualization}
===========================================
"""

// Get reference paths from config
def getRefPath(ref_type) {
    if (params.genomes.containsKey(params.genome)) {
        if (params.genomes[params.genome].containsKey(ref_type)) {
            return params.genomes[params.genome][ref_type]
        }
    }
    return null
}

// Input channel for FASTQ files
input_ch = Channel.fromPath("${params.input_dir}/${params.fastq_pattern}", checkIfExists: true)
                   .map { file -> tuple(file.simpleName, file) }

// Main workflow
workflow {
    // Run core pipeline
    METHYLATION(input_ch)
    
    // Optional analyses
    if (params.run_differential) {
        DIFFERENTIAL(METHYLATION.out.methylation_bed)
    }
    
    if (params.run_region_analysis) {
        REGION(METHYLATION.out.methylation_bed)
    }
    
    if (params.run_visualization) {
        VISUALIZATION(METHYLATION.out.methylation_bed)
    }
    
    // Print workflow completion message
    METHYLATION.out.methylation_bed.view { "Methylation calling completed for sample: ${it}" }
}

// Handle workflow completion
workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'SUCCESS' : 'FAILED' }"
    log.info "Execution duration: $workflow.duration"
}
