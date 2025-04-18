#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process stats {
    input:
    val vcf = file('results/hipstr/output.vcf')

    output:
    
}


workflow STATS{
    take
    
    
    
    val vcf = MAPPING.vcf
    stats(vcf)
}