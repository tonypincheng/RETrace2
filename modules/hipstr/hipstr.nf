#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process HIPSTR_CALLING {
    publishDir "${params.output_dir}/hipstr", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    
    output:
    // path("${sample_id}.vcf"), emit: vcf
    // path("${sample_id}.stats"), emit: stats
    
    script:
    """
    # Call microsatellites using Python wrapper
    python $baseDir/scripts/hipstr/call_hipstr.py \
        --bam ${bam} \
        --reference ${params.genomes_base}/${params.genome}/raw_fasta/${params.genome}.fa \
        --regions ${params.target_bed} \
        --output ${sample_id}.vcf \
        --min-reads ${params.min_reads} \
        --min-allele-freq ${params.min_qual}
    
    # Generate statistics
    python $baseDir/scripts/hipstr/stats.py ${sample_id}.vcf > ${sample_id}.stats
    """
}

workflow HIPSTR {
    take:
    bam
    
    main:
    // Call microsatellite genotypes with HipSTR
    HIPSTR_CALLING(bam)
    
    // emit:
    // vcf = hipstr_results.vcf
    // stats = hipstr_results.stats
} 