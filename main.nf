#!/usr/bin/env nextflow

nextflow.enable.dsl=2

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
      --output_prefix   Prefix for output files (default: ${params.output_prefix})
      --threads         Number of CPU threads per task (default: ${params.threads})
      --memory          Memory to allocate per task (default: ${params.memory})
      --paired_end      Specify if data is paired-end sequencing (default: ${params.paired_end}) [NOTE: Feature not fully implemented]

      --bwa_index_path  Path to BWA index [optional]. If not specified, will use ${params.genome_base}/${params.genome}/bwa-index/${params.genome}.fa
      --ref_fasta       Path to reference FASTA [optional]. If not specified, will use ${params.genome_base}/${params.genome}/raw_fasta/${params.genome}.fa

    Sample Quality parameters:
      --min_targets     Minimum number of microsatellite targets per sample (default: ${params.min_targets})
      --min_cpgs        Minimum number of CpGs per sample for methylation analysis (default: ${params.min_cpgs})

    HipSTR parameters:
      --hipstr_path     Path to HipSTR executable [optional]. If not specified, will use "HipSTR" from PATH
      --min_qual        Minimum quality score per target for HipSTR (default: ${params.min_qual})
      --min_reads       Minimum number of reads per target per sample for HipSTR (default: ${params.min_reads}). Note: this is also used in the STATS module to count targets with min coverage.
      --max_stutter     Maximum stutter ratio for per target for HipSTR (default: ${params.max_stutter}
      --by_chrom        Run HipSTR by chromosome in parallel (default: ${params.by_chrom})
      --snp_vcf         Optional SNP VCF file for HipSTR (default: ${params.snp_vcf ?: "Not used"})
      
    Phylogenetic parameters:
      --dist_metric     Distance metric for phylogenetic tree construction (default: ${params.dist_metric})
      --outgroup        Outgroup for tree rooting (default: ${params.outgroup})
      --color_background Apply color to background instead of node circles (default: ${params.color_background})
      --circular_tree   Render phylogenetic tree in circular layout (default: ${params.circular_tree})
      
    Optional analyses:
      --run_bootstrap   Run bootstrap analysis (default: ${params.run_bootstrap})
      --bootstrap_iterations Number of bootstrap iterations (default: ${params.bootstrap_iterations})
      
      --run_evaluation  Run evaluation (requires ground truth) (default: ${params.run_evaluation})
      --ground_truth    Path to ground truth data (default: ${params.ground_truth})
      
      --run_methylation Run methylation analysis (default: ${params.run_methylation})
      --methylpy_ref    Path prefix for methylpy reference files [optional]. If not specified, will use ${params.genome_base}/${params.genome}/methylpl-ref/${params.genome}      
      --min_reads_per_site Minimum number of reads required per cytosine site (default: 1)
      --min_shared_sites Minimum number of shared sites required for comparison (default: 100)
      --all_cytosines   Use all methylation contexts, not just CpG sites (default: false)
      --celltype_ref_dir Directory containing reference cell type files (required for cell type inference)
      --celltype_ref_pattern Pattern to match files in celltype_ref_dir (default: *.tsv.gz)
      --zscore_threshold Z-score threshold for cell type inference (default: ${params.zscore_threshold})

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

if (!file(params.target_bed).exists()) {
    log.error "Target BED file '${params.target_bed}' does not exist or is not accessible!"
    exit 1
}

// Resolve target_bed to absolute path to ensure modules can access it
params.target_bed_resolved = file(params.target_bed).toString()


// Log pipeline info
log.info"""
===========================================
 RETrace2 Pipeline v1.0
===========================================
Samplesheet        : ${params.samplesheet}
Output directory   : ${params.output_dir}
Reference genome   : ${params.genome}
Profile            : ${workflow.profile ?: 'default'}
Sequencing mode    : ${params.paired_end ? 'Paired-end' : 'Single-end'}
Threads            : ${params.threads}
Memory             : ${params.memory}
Per Sample parameters:
  - Min targets    : ${params.min_targets}
  - Min CpGs       : ${params.min_cpgs}
HipSTR parameters  :
  - Path           : ${params.hipstr_path ?: "HipSTR (from PATH)"}
  - Min quality    : ${params.min_qual}
  - Min reads      : ${params.min_reads}
  - Max stutter    : ${params.max_stutter}
Phylogenetic parameters:
  - Distance metric: ${params.dist_metric}
  - Outgroup       : ${params.outgroup}
  - Color background: ${params.color_background}
  - Circular layout: ${params.circular_tree}
Optional analyses  :
  - Bootstrap      : ${params.run_bootstrap}
  - Evaluation     : ${params.run_evaluation}
  - Methylation    : ${params.run_methylation}
===========================================
"""


// Include modules
include { MAPPING } from './modules/mapping/mapping.nf'
include { STATS } from './modules/stats/stats.nf'
include { HIPSTR } from './modules/hipstr/hipstr.nf'
include { PHYLO } from './modules/phylo/phylo.nf'

// Conditionally include METHYLATION module
if (params.run_methylation) {
    include { METHYLATION } from './modules/methylation/methylation.nf'
    include { INFER_CELLTYPE } from './modules/infer_celltype/infer_celltype.nf'
}

//include { EVALUATION } from './modules/evaluation/evaluation.nf' 


// Main workflow
workflow {
    // Create input channels directly from the samplesheet
    Channel.fromPath(params.samplesheet)
        .splitCsv(header:true)
        .filter { row -> !row.sample_id.startsWith('#') && row.ms_fastq_1 }
        .map { row -> 
            ms_fastq = file(row.ms_fastq_1, checkIfExists: true)
            tuple(row.sample_id, ms_fastq)
        }
        .set { ms_input_ch }
    
    // Create channel for methylation inputs if run_methylation is enabled
    if (params.run_methylation) {
        Channel.fromPath(params.samplesheet)
            .splitCsv(header:true)
            .filter { row -> !row.sample_id.startsWith('#') && row.meth_fastq_1 }
            .map { row ->
                meth_fastq = file(row.meth_fastq_1, checkIfExists: true)
                tuple(row.sample_id, meth_fastq)
            }
            .set { methylation_input_ch }
    }
    
    // Core pipeline
    MAPPING(ms_input_ch)
    
    // Create an empty channel for optional methylation files
    methylation_allc_ch = Channel.empty()
    
    // Run methylation analysis if enabled
    if (params.run_methylation) {
        METHYLATION(methylation_input_ch)
        methylation_allc_ch = METHYLATION.out.allc
    }
    
    
    // Pass both microsatellite BAM and methylation files to STATS
    STATS(MAPPING.out.bam, methylation_allc_ch)
    
    // Run HipSTR for microsatellite genotyping
    HIPSTR(MAPPING.out.bam, STATS.out.sample_stats)
    
    // Run PHYLO for phylogenetic tree construction
    PHYLO(HIPSTR.out.alleleDict, HIPSTR.out.sample_list)
    
    // Run cell type inference if enabled
    if (params.run_methylation && params.celltype_ref_dir) {
        INFER_CELLTYPE(methylation_allc_ch.map { tuple -> tuple[1] })
    }

    // if (params.run_evaluation) {
    //     if (params.ground_truth) {
    //         EVALUATION(tree_ch, file(params.ground_truth))
    //         evaluation_results_ch = EVALUATION.out.results
    //     } else {
    //         log.warn "Evaluation requested but no ground truth provided. Skipping evaluation."
    //     }
    // }
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
        - Microsatellite statistics: ${params.output_dir}/stats/ms_counts/
        - Summary statistics: ${params.output_dir}/stats/
        - HipSTR genotyping: ${params.output_dir}/hipstr/
        - Phylogenetic tree: ${params.output_dir}/phylo/newick-orginal.txt
        """
        
        if (params.run_bootstrap) {
            log.info """
        Bootstrap results:
        - Support values: ${params.output_dir}/phylo/bootstrap_stats.txt
        - Bootstrap trees: ${params.output_dir}/phylo/newwick-bootstrap.txt
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
        - Methylation allc files: ${params.output_dir}/methylation/methylpy/
        - Combined statistics: ${params.output_dir}/stats/ (includes both MS and CpG data)
            """
            
            if (params.celltype_ref_dir) {
                log.info """
        Cell type inference results:
        - Pairwise dissimilarity matrix: ${params.output_dir}/infer_celltype/pd_matrix/
        - Cell type assignments: ${params.output_dir}/infer_celltype/assignments/
        - Z-score matrix: ${params.output_dir}/infer_celltype/assignments/
        - Cell type plots: ${params.output_dir}/infer_celltype/assignments/
                """
            }
        }
        
        log.info "==========================================="
    }
}
