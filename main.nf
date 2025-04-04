#!/usr/bin/env nextflow

params.input_file = "data/input.fastq"
params.output_dir = "results/"

// Create input channel
input_ch = Channel.fromPath(params.input_file, checkIfExists: true)

process processData {
    publishDir params.output_dir, mode: 'copy'

    input:
    path input_file
    
    output:
    path "output.txt", emit: output_ch

    script:
    """
    echo "Processing ${input_file}" > output.txt
    """
}

workflow {
    processData(input_ch)
    processData.out.output_ch.view { "Generated output file: $it" }
}
