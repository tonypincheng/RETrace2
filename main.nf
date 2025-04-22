#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Pipeline parameters
params.samplesheet = null
params.output_dir = "results"
params.target_bed = "${baseDir}/resources/targets/mm39/RETrace2.mm39.1nt10-30bp.92460targets169818probes.bed"
params.threads = 30
params.memory = '100.GB'

// Reference genome parameters
params.genome_base = "/path/to/reference/genome_base"
params.genome = "mm39"
params.ref_fasta = null
params.bwa_index_path = null

// HipSTR parameters
params.min_qual = 0.9
params.min_reads = 10
params.max_stutter = 1.0

// Optional analysis parameters
params.run_bootstrap = false
params.bootstrap_iterations = 100
params.run_evaluation = false
params.ground_truth = null
params.run_methylation = false

// Methylation parameters
params.methylpy_ref = null

// System parameters
params.help = false



// Help message
def helpMessage() {

    log.info"""
    ===========================================
      RETrace2 Pipeline v1.0
    ===========================================
    
    Usage:
      nextflow run main.nf --samplesheet samplesheet.csv --output_dir /path/to/results
    
    Mandatory arguments:
      --samplesheet     CSV file specifying samples and their details (see example in assets/samplesheet.csv)
      --genome_base     Directory containing reference genomes (default: ${params.genome_base})
      --genome          Reference genome: 'mm39' or 'hg38' (default: ${params.genome})
      --target_bed      BED file with targetmicrosatellite regions (default: ${params.target_bed})
    
    Optional arguments:
      --output_dir      Directory for output files (default: ${params.output_dir})
      --threads         Number of CPU threads to use (default: ${params.threads})
      --memory          Memory to allocate for processes (default: ${params.memory})

      --bwa_index_path  Path to BWA index [optional]. If not specified, will use ${params.genome_base}/${params.genome}/bwa-index/${params.genome}.fa
      --ref_fasta Path to reference FASTA [optional]. If not specified, will use ${params.genome_base}/${params.genome}/raw_fasta/${params.genome}.fa
      
      --min_qual        Minimum quality score for HipSTR (default: ${params.min_qual})
      --min_reads       Minimum number of reads for HipSTR (default: ${params.min_reads})
      --max_stutter     Maximum stutter ratio for HipSTR (default: ${params.max_stutter})
      
      --run_bootstrap   Run bootstrap analysis (default: ${params.run_bootstrap})
      --bootstrap_iterations Number of bootstrap iterations (default: ${params.bootstrap_iterations})
      
      --run_evaluation  Run evaluation (requires ground truth) (default: ${params.run_evaluation})
      --ground_truth    Path to ground truth data (default: ${params.ground_truth})
      
      --run_methylation Run methylation analysis (default: ${params.run_methylation})
      --methylpy_ref    Path prefix for methylpy reference files [optional]. If not specified, will use ${params.genome_base}/${params.genome}/methylpl-ref/${params.genome}
      
      --help            Display this help message
    """
}

// Show help message if --help specified
if (params.help) {
    helpMessage()
    exit 0
}

// Check for required parameters
if (!params.samplesheet) {
    log.error "No samplesheet CSV file specified! Please provide a valid samplesheet file using --samplesheet"
    exit 1
}

if (!file(params.samplesheet).exists()) {
    log.error "Samplesheet file '${params.samplesheet}' does not exist or is not accessible!"
    exit 1
}

if (!file(params.genome_base).exists()) {
    log.error "Reference genome base directory '${params.genome_base}' does not exist or is not accessible!"
    exit 1
}


// Log pipeline info
log.info"""
===========================================
 RETrace2 Pipeline v1.0
===========================================
Samplesheet        : ${params.samplesheet}
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


// Include modules
include { MAPPING } from './modules/mapping/mapping.nf'
include { STATS } from './modules/stats/stats.nf'
//include { HIPSTR } from './modules/hipstr/hipstr.nf'
//include { PHYLO } from './modules/phylo/phylo.nf'

// Conditionally include METHYLATION module
if (params.run_methylation) {
    include { METHYLATION } from './modules/methylation/methylation.nf'
}
//include { BOOTSTRAP } from './modules/bootstrap/bootstrap.nf' 
//include { EVALUATION } from './modules/evaluation/evaluation.nf' 


// Main workflow
workflow {
    // Create input channels directly from the samplesheet
    Channel.fromPath(params.samplesheet)
        .splitCsv(header:true)
        .filter { row -> !row.sample_id.startsWith('#') && row.ms_fastq_1 }
        .map { row -> 
            ms_fastq = file(row.ms_fastq_1)
            
            if (!ms_fastq.exists()) {
                log.error "ERROR: Microsatellite FASTQ file does not exist: ${row.ms_fastq_1}"
                exit 1
            }

            tuple(row.sample_id, ms_fastq)
        }
        .set { ms_input_ch }
    
    // Create channel for methylation inputs if run_methylation is enabled
    if (params.run_methylation) {
        Channel.fromPath(params.samplesheet)
            .splitCsv(header:true)
            .filter { row -> !row.sample_id.startsWith('#') && row.meth_fastq_1 }
            .map { row ->
                meth_fastq = file(row.meth_fastq_1)
                
                if (!meth_fastq.exists()) {
                    log.error "ERROR: Methylation FASTQ file does not exist: ${row.meth_fastq_1}"
                    exit 1
                }

                tuple(row.sample_id, meth_fastq)
            }
            .set { methylation_input_ch }
    }
    
    // Core pipeline
    MAPPING(ms_input_ch)
    STATS(MAPPING.out.bam)
    //HIPSTR(MAPPING.out.bam)
    //PHYLO(HIPSTR.out.vcf)    
    
    // // Optional analyses
    // if (params.run_bootstrap) {
    //     BOOTSTRAP(tree_ch, HIPSTR.out.vcf)
    //     bootstrap_support_ch = BOOTSTRAP.out.support
    //     bootstrap_trees_ch = BOOTSTRAP.out.trees
    // }
    
    // if (params.run_evaluation) {
    //     if (params.ground_truth) {
    //         EVALUATION(tree_ch, file(params.ground_truth))
    //         evaluation_results_ch = EVALUATION.out.results
    //     } else {
    //         log.warn "Evaluation requested but no ground truth provided. Skipping evaluation."
    //     }
    // }
    
    if (params.run_methylation) {
        METHYLATION(methylation_input_ch)
        // methylation_results_ch = METHYLATION.out.results
        // methylation_stats_ch = METHYLATION.out.detailed_stats
        // methylation_summary_ch = METHYLATION.out.summary_stats
        // methylation_plot_ch = METHYLATION.out.summary_plot
    }
    
    // Print workflow completion message
    // tree_ch.view { "Phylogenetic tree completed: ${it}" }
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
        - Microsatellite FASTQC: ${params.output_dir}/mapping/fastqc/
        - Microsatellite BAM files: ${params.output_dir}/mapping/bam/
        - Sample stats: ${params.output_dir}/stats/
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
        - Methylation FASTQC: ${params.output_dir}/methylation/fastqc/
        - Methylation allc files: ${params.output_dir}/methylation/allc/

            """
        }
        
        log.info "==========================================="
    }
}
