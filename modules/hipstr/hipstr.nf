#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process filter_bams_by_stats {
    input:
    tuple val(sample_id), path(bam_file), path(bam_index)
    path(sample_stats)
    
    output:
    tuple val(sample_id), path("*.pass.bam"), path("*.pass.bam.bai"), emit: passing_bams, optional: true
    
    script:
    """
    # Parse sample stats to check if this sample passes
    pass_status=\$(awk -F'\\t' -v sample="${sample_id}" '\$1 == sample {print \$4}' ${sample_stats})
    
    if [ "\$pass_status" = "True" ]; then
        # Create symlinks with .pass suffix for passing sample
        ln -s ${bam_file} \$(basename ${bam_file} .bam).pass.bam
        ln -s ${bam_index} \$(basename ${bam_index} .bam.bai).pass.bam.bai
    else
        # Exit with code 0 but produce no output for failing samples
        echo "Sample ${sample_id} did not pass QC (status: \$pass_status). Skipping."
        exit 0
    fi
    """
}

process hipstr_per_chrom {
    tag "$chrom"
    
    input:
    val(chrom)
    path(bam_files)
    path(bam_indices)
    
    output:
    path("${params.output_prefix}.${chrom}.vcf.gz"), emit: vcf
    path("${params.output_prefix}.${chrom}.vcf.gz.tbi"), emit: vcf_index
    path("${params.output_prefix}.${chrom}.log"), emit: log
    
    script:
    // Create comma-separated list for HipSTR
    def bam_list = bam_files.join(',')
    
    // Add SNP VCF option if provided
    def snp_option = params.snp_vcf ? "--snp-vcf ${params.snp_vcf}" : ""
    
    // Handle paired-end vs single-end
    def use_unpaired = params.paired_end ? "" : "--use-unpaired"
    
    // Determine fasta path
    def fasta_path = params.ref_fasta ?: "${params.genome_base}/${params.genome}/raw_fasta/${params.genome}.fa"
    
    // Get HipSTR path from params or use command directly if not specified
    def hipstr_path = params.hipstr_path ?: "HipSTR"
    
    """
    # Run HipSTR for this chromosome
    ${hipstr_path} \\
        --bams ${bam_list} \\
        --fasta ${fasta_path} \\
        --regions ${params.target_bed} \\
        --str-vcf ${params.output_prefix}.${chrom}.vcf.gz \\
        --log ${params.output_prefix}.${chrom}.log \\
        ${use_unpaired} \\
        --no-rmdup \\
        --chrom ${chrom} \\
        ${snp_option}
        
    # Index for merging
    tabix -p vcf ${params.output_prefix}.${chrom}.vcf.gz
    """
}

process merge_vcfs {
    publishDir "${params.output_dir}/hipstr", mode: 'copy'
    
    input:
    path(vcfs)
    path(vcf_indices)
    path(logs)
    
    output:
    path("${params.output_prefix}.vcf"), emit: vcf
    path("${params.output_prefix}.log"), emit: hipstr_log
    
    script:
    """
    # Merge VCFs
    bcftools concat -a -o ${params.output_prefix}.unsorted.vcf ${vcfs}
    
    # Sort the merged VCF
    bcftools sort -o ${params.output_prefix}.vcf ${params.output_prefix}.unsorted.vcf
    
    # Remove the unsorted intermediate file
    rm ${params.output_prefix}.unsorted.vcf
    
    # Combine logs
    cat ${logs} > ${params.output_prefix}.log
    """
}

process hipstr_calling {
    publishDir "${params.output_dir}/hipstr", mode: 'copy'
    
    input:
    path(bam_files)
    path(bam_indices)
    path(sample_stats)
    
    output:
    path("${params.output_prefix}.vcf"), emit: vcf
    path("${params.output_prefix}.log"), emit: hipstr_log
    
    script:
    // Create comma-separated list for HipSTR
    def bam_list = bam_files.join(',')
    
    // Add SNP VCF option if provided
    def snp_option = params.snp_vcf ? "--snp-vcf ${params.snp_vcf}" : ""
    
    // Handle paired-end vs single-end
    def use_unpaired = params.paired_end ? "" : "--use-unpaired"
    
    // Determine fasta path
    def fasta_path = params.ref_fasta ?: "${params.genome_base}/${params.genome}/raw_fasta/${params.genome}.fa"
    
    // Get HipSTR path from params or use command directly if not specified
    def hipstr_path = params.hipstr_path ?: "HipSTR"
    
    """
    # Run HipSTR in standard mode
    ${hipstr_path} \\
        --bams ${bam_list} \\
        --fasta ${fasta_path} \\
        --regions ${params.target_bed} \\
        --str-vcf ${params.output_prefix}.vcf.gz \\
        --log ${params.output_prefix}.log \\
        ${use_unpaired} \\
        --no-rmdup \\
        ${snp_option}
    
    # Unzip the VCF
    gunzip ${params.output_prefix}.vcf.gz
    """
}

workflow HIPSTR {
    take:
    bam
    sample_stats
    
    main:
    // Filter BAM files based on sample stats
    filter_bams_by_stats(bam, sample_stats)
    
    // Extract BAM files and BAI files from passing samples
    filtered_bam_channel = filter_bams_by_stats.out.passing_bams
    
    // Get all filtered BAM and BAI files
    all_bams = filtered_bam_channel.map { sample_id, bam_files, bam_indices -> bam_files }.flatten().collect()
    all_bais = filtered_bam_channel.map { sample_id, bam_files, bam_indices -> bam_indices }.flatten().collect()
    
    if (params.by_chrom) {
        // Extract unique chromosomes from BED file
        chroms_ch = Channel.fromPath(params.target_bed)
            .splitCsv(sep: '\t')
            .map { row -> row[0] }
            .unique()
        
        // Run HipSTR per chromosome in parallel
        hipstr_per_chrom(chroms_ch, all_bams, all_bais)
        
        // Merge results
        merge_vcfs(
            hipstr_per_chrom.out.vcf.collect(),
            hipstr_per_chrom.out.vcf_index.collect(),
            hipstr_per_chrom.out.log.collect()
        )
        
        vcf_output = merge_vcfs.out.vcf
        log_output = merge_vcfs.out.hipstr_log
    } else {
        // Run standard HipSTR
        hipstr_calling(all_bams, all_bais, sample_stats)
        
        vcf_output = hipstr_calling.out.vcf
        log_output = hipstr_calling.out.hipstr_log
    }
    
    emit:
    vcf = vcf_output
    hipstr_log = log_output
} 