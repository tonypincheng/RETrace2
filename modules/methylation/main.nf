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
    tuple val(sample_id), path("${sample_id}_trimmed.fq.gz"), emit: trimmed_reads
    
    script:
    """
    trim_galore \
        --quality 30 \
        --phred33 \
        --stringency 3 \
        --length 36 \
        --output_dir . \
        -j ${task.cpus} \
        ${reads}
    """
}

// Process 4: Alignment
process alignment {
    publishDir "${params.output_dir}/aligned", mode: 'copy'
    
    input:
    tuple val(sample_id), path(trimmed_reads)
    
    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), path("${sample_id}.stats"), emit: aligned_reads
    
    script:
    """
    bwa mem -t ${task.cpus} /home/tcheng/Projects/GenomeDB/mm39/bwa-index/mm39.fa ${trimmed_reads} | \
    samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam -
    samtools flagstat ${sample_id}.sorted.bam > ${sample_id}.stats
    samtools index ${sample_id}.sorted.bam
    """
}

// Process 5: Generate Summary Stats
process generate_stats {
    publishDir "${params.output_dir}/stats", mode: 'copy'
    
    input:
    path(stats_files)
    
    output:
    path("all_stats.txt"), emit: summary_stats
    
    script:
    """
    if [ -f all_stats.txt ]; then rm all_stats.txt; fi
    
    for f in ${stats_files}; do
        fbname=\$(basename "\$f" .stats)
        echo "\$fbname" >> all_stats.txt
        
        # Extract the total number of primary reads
        n_read=\$(grep " in total " "\$f" | awk '{print \$1}')
        echo "Total Reads: \$n_read" >> all_stats.txt
        
        # Extract the primary mapped percentage
        align_rate=\$(grep " mapped (" "\$f" | awk -F '[()]' '{print \$2}')
        echo "Alignment Rate: \$align_rate" >> all_stats.txt
    done
    """
}

// Process 5: Methylation calling with methylpl
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
    
    // Generate QC report
    multiqc(
        fastqc.out.fastqc_results.collect().ifEmpty([])
    )
    
    // Run preprocessing
    preprocess(reads)
    
    // Run alignment
    alignment(preprocess.out.trimmed_reads)
    
    // Generate summary stats
    stats = generate_stats(alignment.out.aligned_reads.collect{it[2]})
    
    // Run methylation calling
    methylation(alignment.out.aligned_reads.collect{it[0,1]})
    
    // Run methylation analysis
    meth_results = METHYLATION_ANALYSIS(reads)
    
    // Infer cell types
    cell_type_results = CELL_TYPE_INFERENCE(meth_results.bed)
    
    emit:
    bed = meth_results.bed
    stats = meth_results.stats
    predictions = cell_type_results.predictions
    plot = cell_type_results.plot
    summary_stats = stats.summary_stats
} 