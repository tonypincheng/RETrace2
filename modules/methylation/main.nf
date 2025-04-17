#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Process 1: Quality Control
process fastqc {
    publishDir "${params.output_dir}/fastqc", mode: 'copy'
    container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
    conda "bioconda::fastqc=0.12.1"
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path("${sample_id}_fastqc.{zip,html}"), emit: fastqc_results
    
    script:
    """
    fastqc -t ${task.cpus} -o . ${reads}
    """
}

// Process 2: MultiQC Report
process multiqc {
    publishDir "${params.output_dir}/fastqc", mode: 'copy'
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
process preprocess {
    publishDir "${params.output_dir}/trimmed", mode: 'copy'
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
        ${reads}
    """
}

// Process 4: Run Methylpy
process methylpy {
    publishDir "${params.output_dir}/methylpy", mode: 'copy'
    container "community.wave.seqera.io/library/pip_methylpy:ae44180dc4227f32"
    conda "bioconda::methylpy=1.4.7"
    
    input:
    tuple val(sample_id), path(trimmed_reads)
    
    output:
    path("log/${sample_id}_methylpy.log"), emit: log_file
    path("allc/*"), emit: allc
    
    script:
    """
    mkdir -p log

    methylpy single-end-pipeline \
        --read-files ${trimmed_reads} \
        --sample ${sample_id} \
        --forward-ref ${params.genomes_base}/${params.genome}/methylpl-ref/${params.genome}_f \
        --reverse-ref ${params.genomes_base}/${params.genome}/methylpl-ref/${params.genome}_r \
        --ref-fasta ${params.genomes_base}/${params.genome}/raw_fasta/${params.genome}.fa \
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
process analyze_methylpy_stats {
    publishDir "${params.output_dir}/stats", mode: 'copy'
    container "quay.io/biocontainers/python:3.9"
    conda "conda-forge::python=3.9 pandas matplotlib seaborn"
    
    input:
    path(log_files)
    path(tsv_files)
    
    output:
    path("detailed_stats.csv"), emit: detailed_stats
    path("summary_stats.txt"), emit: summary_stats
    path("methylpy_summary_plot.pdf"), emit: summary_plot
    
    script:
    """
    python modules/methylation/analyze_methylpy_stats.py \
        --log ${log_files} \
        --tsv-dir . \
        --output-dir .
    """
}

// Module workflow
workflow METHYLATION {
    take:
    reads
    
    main:
    // Run QC on methylation reads
    fastqc(reads)
    
    // Generate QC report
    multiqc(fastqc.out.fastqc_results.collect().ifEmpty([]))
    
    // Run preprocessing on methylation reads
    preprocess(reads)
    
    // Run methylpy
    methylpy(preprocess.out.trimmed_reads)
    
    // // Analyze methylpy output and stats
    // stats = analyze_methylpy_stats(methylpy_results.log_file.collect(), methylpy_results.results)
    
    // emit:
    // log_file = methylpy_results.log_file
    // results = methylpy_results.results
    // detailed_stats = stats.detailed_stats
    // summary_stats = stats.summary_stats
    // summary_plot = stats.summary_plot
}

// Add this at the end of the file
workflow.onComplete {
    println "Methylation module completed"
} 