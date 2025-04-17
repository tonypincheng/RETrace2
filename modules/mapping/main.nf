#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process FASTQC {
    tag "${sample_id}"
    publishDir "${params.output_dir}/fastqc", mode: 'copy'
    container "quay.io/biocontainers/fastqc:${params.versions.fastqc}"
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path("${sample_id}_fastqc.html"), emit: fastqc
    
    script:
    """
    fastqc $reads -o .
    """
}

process TRIM_GALORE {
    tag "${sample_id}"
    publishDir "${params.output_dir}/trimmed", mode: 'copy'
    container "quay.io/biocontainers/trim-galore:${params.versions.trim_galore}"
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fq.gz"), emit: trimmed
    
    script:
    """
    trim_galore \
        --quality 30 \
        --phred33 \
        --stringency 3 \
        --length 36 \
        --output_dir . \
        $reads
    """
}

process BWA_MEM {
    tag "${sample_id}"
    publishDir "${params.output_dir}/bam", mode: 'copy'
    container "quay.io/biocontainers/bwa:${params.versions.bwa}"
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), emit: bam
    path("${sample_id}.stats"), emit: stats
    
    script:
    // Use direct reference path if provided, otherwise build from genomes_base structure
    def reference = params.bwa_index_path ?: "${params.genomes_base}/${params.genome}/bwa-index/${params.genome}.fa"
    """
    bwa mem -t ${params.threads} ${reference} $reads | \
    samtools sort -@${params.threads} -o ${sample_id}.sorted.bam
    
    samtools flagstat ${sample_id}.sorted.bam > ${sample_id}.stats
    samtools index ${sample_id}.sorted.bam
    """
}

process MULTIQC {
    publishDir "${params.output_dir}/multiqc", mode: 'copy'
    container "quay.io/biocontainers/multiqc:${params.versions.multiqc}"
    
    input:
    path(fastqc_files)
    path(stats_files)
    
    output:
    path("multiqc_report.html"), emit: report
    
    script:
    """
    multiqc --force --filename multiqc_report.html .
    """
}

workflow MAPPING {
    take:
    reads
    
    main:
    // Run FastQC
    fastqc_results = FASTQC(reads)
    
    // Generate QC report
    multiqc_report = MULTIQC(
        fastqc_results.fastqc.collect()
    )
    
    // Trim reads
    trimmed_reads = TRIM_GALORE(reads)
    
    // Align reads
    bam_files = BWA_MEM(trimmed_reads)
    
    emit:
    bam = bam_files.bam
    stats = bam_files.stats
    fastqc = fastqc_results.fastqc
    multiqc = multiqc_report.report
} 