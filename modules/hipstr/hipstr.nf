#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process hipstr_calling {
    publishDir "${params.output_dir}/hipstr", mode: 'copy'
    //container "community.wave.seqera.io/library/hipstr:latest"
    //conda "bioconda::hipstr bioconda::bcftools=1.21 bioconda::samtools=1.21 bioconda::tabix=1.21 conda-forge::pandas=2.2.1"
    

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
    // By default, we'll use --use-unpaired unless paired_end=true is specified
    def use_unpaired = params.paired_end ? "" : "--use-unpaired"
    
    // Determine fasta path
    def fasta_path = params.ref_fasta ?: "${params.genome_base}/${params.genome}/raw_fasta/${params.genome}.fa"
    
    // Get HipSTR path from params or use command directly if not specified
    def hipstr_path = params.hipstr_path ?: "HipSTR"
    
    // Choose between running by chromosome or standard mode
    if (params.by_chrom) {
        """
        #!/bin/bash
        # Extract unique chromosome names from the BED file
        cut -f1 ${params.target_bed} | sort | uniq > chroms.txt
        
        # Create directory for chromosome VCFs
        mkdir -p chrom_vcfs
        
        # Run HipSTR on each chromosome
        cat chroms.txt | while read chrom; do
          echo "Processing chromosome: \$chrom"
          
          # Run HipSTR for this chromosome
          ${hipstr_path} \\
            --bams ${bam_list} \\
            --fasta ${fasta_path} \\
            --regions ${params.target_bed} \\
            --str-vcf chrom_vcfs/${params.output_prefix}.\$chrom.vcf.gz \\
            --log chrom_vcfs/${params.output_prefix}.\$chrom.log \\
            ${use_unpaired} \\
            --no-rmdup \\
            --chrom \$chrom \\
            ${snp_option}
            
          # Index for merging
          tabix -p vcf chrom_vcfs/${params.output_prefix}.\$chrom.vcf.gz
        done
        
        # Create file list for bcftools
        find chrom_vcfs -name "*.vcf.gz" | sort > vcf_list.txt
        
        # Combine VCFs
        bcftools concat -a --file-list vcf_list.txt -o ${params.output_prefix}.vcf
        
        # Combine logs
        cat chrom_vcfs/*.log > ${params.output_prefix}.log
        
        # Clean up temporary files
        rm -rf chrom_vcfs
        rm chroms.txt vcf_list.txt
        """
    } else {
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
}

workflow HIPSTR {
    take:
    bam
    sample_stats
    
    main:
    // Extract BAM files and BAI files
    all_bams = bam.map { sample_id, bam_files, bam_indices -> bam_files }.flatten().collect()
    all_bais = bam.map { sample_id, bam_files, bam_indices -> bam_indices }.flatten().collect()
    
    // Run HipSTR calling process with all BAMs and their indices
    hipstr_calling(all_bams, all_bais, sample_stats)
    
    emit:
    vcf = hipstr_calling.out.vcf
    hipstr_log = hipstr_calling.out.log
} 