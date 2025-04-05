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

// Process 2: Alignment
process alignment {
    publishDir "${params.output_dir}/aligned", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: aligned_reads
    
    script:
    def fasta = getRefPath("fasta")
    def bwa_index = getRefPath("bwa_index")
    
    """
    bwa mem -t ${task.cpus} ${bwa_index}/${params.genome} ${reads} | \
    samtools sort -@ ${task.cpus} -o ${sample_id}.bam -
    samtools index ${sample_id}.bam
    """
}

// Process 3: HipSTR calling
process call_hipstr {
    publishDir "${params.output_dir}/hipstr", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    
    output:
    path "${sample_id}_microsatellites.vcf", emit: microsatellites
    
    script:
    def fasta = getRefPath("fasta")
    
    """
    python $baseDir/scripts/core/call_hipstr.py \
      --bam ${bam} \
      --reference ${fasta} \
      --output ${sample_id}_microsatellites.vcf \
      --threads ${task.cpus}
    """
}

// Process 4: Build phylogenetic tree
process build_tree {
    publishDir "${params.output_dir}/phylogeny", mode: 'copy'
    
    input:
    path microsatellites
    
    output:
    path "phylogenetic_tree.nwk", emit: tree
    path "tree_metadata.tsv", emit: metadata
    
    script:
    """
    Rscript $baseDir/scripts/core/build_tree.R \
      --input "${microsatellites}" \
      --tree "phylogenetic_tree.nwk" \
      --metadata "tree_metadata.tsv"
    """
}

// Module workflow
workflow CORE {
    take:
    input_ch
    
    main:
    // Run QC
    fastqc(input_ch)
    
    // Run alignment
    alignment(input_ch)
    
    // Call microsatellites with HipSTR
    call_hipstr(alignment.out.aligned_reads)
    
    // Build phylogenetic tree
    build_tree(call_hipstr.out.microsatellites.collect())
    
    emit:
    tree = build_tree.out.tree
    metadata = build_tree.out.metadata
} 