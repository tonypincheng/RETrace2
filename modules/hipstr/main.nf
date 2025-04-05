#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process HIPSTR_CALLING {
    tag "${sample_id}"
    publishDir "${params.output_dir}/hipstr", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    path("${sample_id}.vcf"), emit: vcf
    path("${sample_id}.stats"), emit: stats
    
    script:
    """
    HipSTR \
        --bams $bam \
        --fasta ${params.genomes_base}/${params.genome}/raw_fasta/${params.genome}.fa \
        --regions ${params.target_bed} \
        --min-reads ${params.min_reads} \
        --min-qual ${params.min_qual} \
        --max-stutter ${params.max_stutter} \
        --output ${sample_id}.vcf \
        --log ${sample_id}.log
    
    # Generate statistics
    python $baseDir/scripts/hipstr/stats.py ${sample_id}.vcf > ${sample_id}.stats
    """
}

workflow HIPSTR {
    take:
    bam_files
    
    main:
    // Call microsatellites with HipSTR
    hipstr_results = HIPSTR_CALLING(bam_files)
    
    emit:
    vcf = hipstr_results.vcf
    stats = hipstr_results.stats
} 