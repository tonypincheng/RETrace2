#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Process 1: Quality Control
process FASTQC {
    publishDir "${params.output_dir}/mapping/fastqc", mode: 'copy'
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
process MULTIQC {
    publishDir "${params.output_dir}/mapping/fastqc", mode: 'copy'
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

process TRIM_GALORE {
    publishDir "${params.output_dir}/mapping/trimmed", mode: 'copy'
    container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"
    conda "bioconda::trim-galore=0.6.10"
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fq.gz"), emit: trimmed_reads
    
    script:
    """
    trim_galore \
        --quality 30 \
        --phred33 \
        --stringency 3 \
        --length 36 \
        --cores ${task.cpus} \
        --output_dir . \
        ${reads}
    """
}

process BWA_MEM {
    publishDir "${params.output_dir}/mapping/bam", mode: 'copy'
    container "quay.io/biocontainers/bwa:0.7.19--h577a1d6_0"
    conda "bioconda::bwa=0.7.19"

    input:
    tuple val(sample_id), path(trimmed_reads)
    
    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: bam
    path("${sample_id}.stats"), emit: stats
    
    script:
    // Use direct reference path if provided, otherwise build from genomes_base structure
    def reference = params.bwa_index_path ?: "${params.genomes_base}/${params.genome}/bwa-index/${params.genome}.fa"
    """
    bwa mem -t ${task.cpus} ${reference} ${trimmed_reads} | \
    samtools sort -@${task.cpus} -o ${sample_id}.sorted.bam
    
    samtools index ${sample_id}.sorted.bam
    samtools flagstat ${sample_id}.sorted.bam > ${sample_id}.stats
    """
}


workflow MAPPING {
    take:
    reads
    
    main:
    // Run FastQC
    FASTQC(reads)
    
    // Generate QC report
    MULTIQC(FASTQC.out.fastqc_results.collect())
    
    // Trim reads
    TRIM_GALORE(reads)
    
    // Align reads
    BWA_MEM(TRIM_GALORE.out.trimmed_reads)
} 