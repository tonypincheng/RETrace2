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

// Process 2: Read Preprocessing
process preprocess {
    publishDir "${params.output_dir}/trimmed", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq"), emit: trimmed_reads
    
    script:
    """
    trimmomatic SE -threads ${task.cpus} \
      ${reads} ${sample_id}_trimmed.fastq \
      ILLUMINACLIP:$baseDir/assets/adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

// Process 3: Alignment
process alignment {
    publishDir "${params.output_dir}/aligned", mode: 'copy'
    
    input:
    tuple val(sample_id), path(trimmed_reads)
    
    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: aligned_reads
    
    script:
    def fasta = getRefPath("fasta")
    def bwa_index = getRefPath("bwa_index")
    
    """
    bwa mem -t ${task.cpus} ${bwa_index}/${params.genome} ${trimmed_reads} | \
    samtools sort -@ ${task.cpus} -o ${sample_id}.bam -
    samtools index ${sample_id}.bam
    """
}

// Process 4: Methylation calling with methylpl
process methylation {
    publishDir "${params.output_dir}/methylation", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    
    output:
    path "${sample_id}_methylation.bed", emit: methylation_bed
    
    script:
    def fasta = getRefPath("fasta")
    def methylpl_index = getRefPath("methylpl_index")
    
    """
    methylpl call \
      --bam ${bam} \
      --reference ${fasta} \
      --index-dir ${methylpl_index} \
      --threads ${task.cpus} \
      --output ${sample_id}_methylation.bed
    """
}

// Process 5: MultiQC Report
process multiqc {
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    path('fastqc/*')
    path('methylation/*')
    
    output:
    path "multiqc_report.html", emit: multiqc_report
    
    script:
    """
    multiqc .
    """
}

process METHYLATION_ANALYSIS {
    tag "${sample_id}"
    publishDir "${params.output_dir}/methylation", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path("${sample_id}.bed"), emit: bed
    path("${sample_id}.stats"), emit: stats
    
    script:
    """
    # Run methylation analysis
    python $baseDir/scripts/methylation/analyze_methylation.py \
        --input $reads \
        --genome ${params.genomes_base}/${params.genome}/raw_fasta/${params.genome}.fa \
        --output ${sample_id}.bed \
        --stats ${sample_id}.stats
    """
}

process CELL_TYPE_INFERENCE {
    tag "cell_type"
    publishDir "${params.output_dir}/cell_type", mode: 'copy'
    
    input:
    path(bed_files)
    
    output:
    path("cell_type_predictions.txt"), emit: predictions
    path("cell_type_plot.pdf"), emit: plot
    
    script:
    """
    # Infer cell types from methylation patterns
    python $baseDir/scripts/methylation/infer_cell_types.py \
        --input $bed_files \
        --output cell_type_predictions.txt \
        --plot cell_type_plot.pdf
    """
}

// Module workflow
workflow METHYLATION {
    take:
    reads
    
    main:
    // Run QC
    fastqc(reads)
    
    // Run preprocessing
    preprocess(reads)
    
    // Run alignment and methylation calling
    alignment(preprocess.out.trimmed_reads)
    methylation(alignment.out.aligned_reads)
    
    // Run methylation analysis
    meth_results = METHYLATION_ANALYSIS(reads)
    
    // Infer cell types
    cell_type_results = CELL_TYPE_INFERENCE(meth_results.bed)
    
    // Generate MultiQC report
    multiqc(
        fastqc.out.fastqc_results.collect().ifEmpty([]),
        methylation.out.methylation_bed.collect().ifEmpty([])
    )
    
    emit:
    bed = meth_results.bed
    stats = meth_results.stats
    predictions = cell_type_results.predictions
    plot = cell_type_results.plot
} 