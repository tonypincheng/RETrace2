#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Process 1: Quality Control
process METH_FASTQC {
    publishDir "${params.output_dir}/methylation/fastqc", mode: 'copy'
    container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
    conda "bioconda::fastqc=0.12.1"
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path("${reads.simpleName}_fastqc.{zip,html}"), emit: fastqc_results
    
    script:
    """
    fastqc -t ${task.cpus} -o . ${reads}
    """
}

// Process 2: MultiQC Report
process METH_MULTIQC {
    publishDir "${params.output_dir}/methylation/fastqc", mode: 'copy'
    container "quay.io/biocontainers/multiqc:1.28--pyhdfd78af_0"
    conda "bioconda::multiqc=1.28"
   
    input:
    path('fastqc/*')
    
    output:
    path("multiqc_report.html")
    path("multiqc_data")
    
    script:
    """
    multiqc .
    """
}

// Process 3: Read Preprocessing
process METH_TRIM_GALORE {
    publishDir "${params.output_dir}/methylation/trimmed", mode: 'copy'
    container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"
    conda "bioconda::trim-galore=0.6.10"

    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fq.gz"), emit: trimmed_reads
    
    script:
    """
    trim_galore \
        --quality 20 \
        --phred33 \
        --stringency 3 \
        --length 36 \
        --rrbs \
        --cores ${task.cpus} \
        --output_dir . \
        --basename ${sample_id} \
        ${reads}
    """
}

// Process 4: Run Methylpy
process METHYLPY {
    publishDir "${params.output_dir}/methylation/methylpy", mode: 'copy'
    container "community.wave.seqera.io/library/pip_methylpy:ae44180dc4227f32"
    conda "bioconda::methylpy=1.4.7"
    
    input:
    tuple val(sample_id), path(trimmed_reads)
    
    output:
    path("log/${sample_id}_methylpy.log"), emit: log_file
    tuple val(sample_id), path("allc/*bam"), path("allc/*bam.bai"), emit: methylpy_bam
    tuple val(sample_id), path("allc/allc*tsv.gz"), emit: allc
    
    script:
    def ref_prefix = params.methylpy_ref ?: "${params.genome_base}/${params.genome}/methylpl-ref/${params.genome}"
    def ref_fasta = params.ref_fasta ?: "${params.genome_base}/${params.genome}/raw_fasta/${params.genome}.fa"
    """
    mkdir -p log

    methylpy single-end-pipeline \
        --read-files ${trimmed_reads} \
        --sample ${sample_id} \
        --forward-ref ${ref_prefix}_f \
        --reverse-ref ${ref_prefix}_r \
        --ref-fasta ${ref_fasta} \
        --num-procs ${task.cpus} \
        --remove-clonal False \
        --min-qual-score 30 \
        --trim-reads False \
        --path-to-picard="picard" \
        --path-to-output allc/ \
        > log/${sample_id}_methylpy.log 2>&1
    """
}

// Process 5: Analyze methylpy output and generate summary statistics
// process CALCUALE_PD_MATRIX {
//     publishDir "${params.output_dir}/methylation/infer_celltype", mode: 'copy'

    
//     input:
//     path(tsv_files)
    
//     output:

    
//     script:
//     """

//     """
// }

// Module workflow
workflow METHYLATION {
    take:
    reads
    
    main:
    // Run QC on methylation reads
    METH_FASTQC(reads)
    
    // Generate QC report
    METH_MULTIQC(METH_FASTQC.out.fastqc_results.collect().ifEmpty([]))
    
    // Run preprocessing on methylation reads
    METH_TRIM_GALORE(reads)
    
    // Run methylpy
    METHYLPY(METH_TRIM_GALORE.out.trimmed_reads)
    
    // // Analyze methylpy output and stats
    // stats = ANALYZE_METHYLPY_STATS(methylpy_results.log_file.collect(), methylpy_results.results)
    
    emit:
    allc = METHYLPY.out.allc
}