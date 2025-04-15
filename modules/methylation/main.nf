#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Process 1: Quality Control
process fastqc {
    publishDir "${params.output_dir}/fastqc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "${sample_id}_fastqc.{zip,html}", emit: fastqc_results
    
    script:
    """
    fastqc -t ${task.cpus} -o . ${reads}
    """
}

// Process 2: MultiQC Report
process multiqc {
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    path('fastqc/*')
    
    output:
    path "multiqc_report.html", emit: multiqc_report
    
    script:
    """
    multiqc .
    """
}

// Process 3: Read Preprocessing
process preprocess {
    publishDir "${params.output_dir}/trimmed", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_R1_trimmed.fq.gz"), emit: trimmed_reads
    
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
    publishDir "${params.output_dir}/methylation", mode: 'copy'
    
    input:
    tuple val(sample_id), path(trimmed_reads)
    
    output:
    path("${sample_id}.log"), emit: log
    path("${sample_name}/*"), emit: results
    
    script:
    """
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
        --path-to-output . \
        > ${sample_id}.log 2>&1
    """
}

// Process 5: Analyze methylpy output and generate summary statistics
process analyze_methylpy_stats {
    publishDir "${params.output_dir}/stats", mode: 'copy'
    
    input:
    path(log_files)
    path(tsv_files)
    
    output:
    path("detailed_stats.csv"), emit: detailed_stats
    path("summary_stats.txt"), emit: summary_stats
    path("methylpy_summary_plot.pdf"), emit: summary_plot
    
    script:
    """
    python $baseDir/scripts/methylation/analyze_methylpy_stats.py \
        --log ${log_files} \
        --tsv-dir . \
        --output-dir .
    """
}


// Module workflow
workflow METHYLATION {
    take:
    methylation_reads
    
    main:
    // Run QC on methylation reads
    fastqc(methylation_reads)
    
    // Generate QC report
    multiqc(
        fastqc.out.fastqc_results.collect().ifEmpty([])
    )
    
    // Only run methylation analysis if specified
    if (params.run_methylation) {
        // Run preprocessing on methylation reads
        preprocess(methylation_reads)
        
        // Run methylpy
        methylpy_results = methylpy(preprocess.out.trimmed_reads)
        
        // Analyze methylpy output and stats
        stats = analyze_methylpy_stats(methylpy_results.log.collect(), methylpy_results.results)
        
        emit:
        log = methylpy_results.log
        results = methylpy_results.results
        detailed_stats = stats.detailed_stats
        summary_stats = stats.summary_stats
        summary_plot = stats.summary_plot
    } else {
        emit:
        log = []
        results = []
        detailed_stats = []
        summary_stats = []
        summary_plot = []
    }
} 