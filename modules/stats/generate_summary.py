#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import numpy as np
import gzip
import matplotlib.pyplot as plt
import seaborn as sns

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Generate summary statistics from MS and CpG counts')
    parser.add_argument('--ms-counts', nargs='+', help='Microsatellite count files')
    parser.add_argument('--allc-files', nargs='*', help='ALLC methylation files (optional)')
    parser.add_argument('--min-reads-per-target', type=int, default=10, help='Minimum reads per microsatellite target')
    parser.add_argument('--min-targets-per-sample', type=int, default=10, help='Minimum targets per sample')
    parser.add_argument('--min-cpgs-per-sample', type=int, default=1000, help='Minimum CpGs per sample')
    parser.add_argument('--output-dir', required=True, help='Output directory for summary files')
    parser.add_argument('--no-plots', action='store_true', help='Disable plot generation')
    return parser.parse_args()

def load_ms_counts(ms_files, min_reads_per_target=10):
    """Load microsatellite count data from files"""
    data = {}
    
    for file_path in ms_files:
        sample_id = os.path.basename(file_path).replace('_ms_counts.txt', '')
        data[sample_id] = {'ms_data': {}}
        
        target_count = 0
        target_with_sufficient_reads = 0
        
        with open(file_path, 'r') as f:
            for line in f:
                # Skip empty lines or header
                if not line.strip() or line.startswith('target_id'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    target_id = parts[0]
                    depth = int(parts[1])
                    
                    data[sample_id]['ms_data'][target_id] = depth
                    target_count += 1
                    
                    if depth >= min_reads_per_target:
                        target_with_sufficient_reads += 1
        
        # Add summary statistics
        data[sample_id]['ms_data']['total_targets'] = target_count
        data[sample_id]['ms_data']['targets_with_min_reads'] = target_with_sufficient_reads
    
    return data

def load_cpg_counts(allc_files):
    """Load CpG methylation data from ALLC files"""
    data = {}
    
    for file_path in allc_files:
        # Extract sample ID from file name
        basename = os.path.basename(file_path)
        sample_id = basename.replace('allc_', '').replace('.tsv.gz', '')
        
        data[sample_id] = {'cpg_data': {}}
        
        cpg_count = 0
        methylated_count = 0
        
        # Process ALLC file to count CpGs
        # Check if file is gzipped
        is_gzipped = file_path.endswith('.gz')
        
        if is_gzipped:
            f = gzip.open(file_path, 'rt')  # Open in text mode
        else:
            f = open(file_path, 'r')
            
        try:
            for line in f:
                if line.startswith('#'):  # Skip header lines
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) >= 7:  # ALLC format typically has at least 7 columns
                    context = fields[3]
                    # Count only CpG context (typically marked as 'CG')
                    if context.startswith('CG'):
                        cpg_count += 1
                        # Typically columns are: chr, pos, strand, context, mc, cov, methylated
                        # Where mc = methylated read count, cov = coverage
                        try:
                            mc = int(fields[4])
                            cov = int(fields[5])
                            if cov > 0 and mc/cov > 0.5:  # Consider methylated if >50% reads show methylation
                                methylated_count += 1
                        except (ValueError, IndexError):
                            pass  # Skip malformed lines
        finally:
            f.close()
        
        # Store the counts
        data[sample_id]['cpg_data']['total_cpgs'] = cpg_count
        data[sample_id]['cpg_data']['methylated_cpgs'] = methylated_count
        if cpg_count > 0:
            data[sample_id]['cpg_data']['methylation_rate'] = methylated_count / cpg_count
        else:
            data[sample_id]['cpg_data']['methylation_rate'] = 0
    
    return data

def merge_data(ms_data, cpg_data=None):
    """
    Merge microsatellite and CpG data into a single dataset.
    
    Returns a combined dictionary with all samples and their respective data.
    """
    # Start with a copy of the microsatellite data
    all_data = ms_data.copy()
    
    # If CpG data is provided, merge it with the microsatellite data
    if cpg_data:
        for sample_id, sample_data in cpg_data.items():
            if sample_id in all_data:
                # Sample exists in both datasets - add CpG data
                all_data[sample_id].update(sample_data)
            else:
                # Sample only exists in CpG data - add it with empty ms_data
                all_data[sample_id] = {'ms_data': {}}
                all_data[sample_id].update(sample_data)
    
    return all_data

def create_plots(df, output_dir):
    """Create visualization plots for summary statistics"""
    # Set the style
    sns.set_theme(style="whitegrid")
    
    # Create figures directory if it doesn't exist
    figures_dir = os.path.join(output_dir, 'figures')
    os.makedirs(figures_dir, exist_ok=True)
    
    # Sort dataframe by sample_id
    df_sorted = df.sort_values('sample_id')
    
    # Plot 1: Bar plot of MS targets per sample
    plt.figure(figsize=(max(10, len(df)*0.3), 8))
    ax = sns.barplot(data=df_sorted, x='sample_id', y='ms_targets_with_min_reads', hue='pass', hue_order=[True, False])
    plt.axhline(y=df_sorted.iloc[0]['min_targets_per_sample'], color='red', linestyle='--',
               label=f'Min threshold: {df_sorted.iloc[0]["min_targets_per_sample"]}')
    plt.title('Microsatellite Targets per Sample', fontsize=18)
    plt.xlabel('Sample ID', fontsize=16)
    plt.ylabel('Number of Microsatellite Targets\n', fontsize=16)
    plt.xticks(rotation=45, ha='right')
    plt.legend()
    plt.tight_layout()
    # Save in both formats
    plt.savefig(os.path.join(figures_dir, 'ms_targets_barplot.png'), dpi=300)
    plt.savefig(os.path.join(figures_dir, 'ms_targets_barplot.pdf'))
    plt.close()
    
    # Check if CpG data is available
    if 'cpg_count' in df.columns and df['cpg_count'].sum() > 0:
        # Plot 2: Bar plot of CpG counts per sample (using the same sample_id sorted order)
        plt.figure(figsize=(max(10, len(df)*0.3), 8))
        ax = sns.barplot(data=df_sorted, x='sample_id', y='cpg_count', hue='pass', hue_order=[True, False])
        plt.axhline(y=df_sorted.iloc[0]['min_cpgs_per_sample'], color='red', linestyle='--',
                   label=f'Min threshold: {df_sorted.iloc[0]["min_cpgs_per_sample"]}')
        plt.title('CpG Counts per Sample', fontsize=18)
        plt.xlabel('Sample ID', fontsize=16)
        plt.ylabel('Number of CpGs\n', fontsize=16)
        plt.xticks(rotation=45, ha='right')
        plt.legend()
        plt.tight_layout()
        # Save in both formats
        plt.savefig(os.path.join(figures_dir, 'cpg_count_barplot.png'), dpi=300)
        plt.savefig(os.path.join(figures_dir, 'cpg_count_barplot.pdf'))
        plt.close()
    
    return figures_dir

def main():
    args = parse_args()
    
    # Load MS count data
    ms_data = load_ms_counts(args.ms_counts, args.min_reads_per_target)
    
    # Load CpG count data if available
    cpg_data = None
    if args.allc_files and len(args.allc_files) > 0:
        cpg_data = load_cpg_counts(args.allc_files)
    
    # Merge data
    merged_data = merge_data(ms_data, cpg_data)
    
    # Process data for output - creating a simple DataFrame
    rows = []
    for sample_id, values in merged_data.items():
        row = {'sample_id': sample_id}
        
        # Get microsatellite targets count with min reads
        ms_targets_with_min_reads = 0
        if 'ms_data' in values:
            ms_targets_with_min_reads = values['ms_data'].get('targets_with_min_reads', 0)
        
        # Get CPG count
        cpg_count = 0
        methylation_rate = 0
        if 'cpg_data' in values and 'total_cpgs' in values['cpg_data']:
            cpg_count = values['cpg_data']['total_cpgs']
            methylation_rate = values['cpg_data'].get('methylation_rate', 0)
        
        # Determine if sample passes criteria
        passes_ms = ms_targets_with_min_reads >= args.min_targets_per_sample
        passes_cpg = cpg_count >= args.min_cpgs_per_sample if args.allc_files else True
        passes = passes_ms and passes_cpg
        
        # Store data
        row['ms_targets_with_min_reads'] = ms_targets_with_min_reads
        row['cpg_count'] = cpg_count
        row['methylation_rate'] = methylation_rate
        row['pass'] = passes
        row['min_targets_per_sample'] = args.min_targets_per_sample
        row['min_cpgs_per_sample'] = args.min_cpgs_per_sample
        
        rows.append(row)
    
    # Create DataFrame
    df = pd.DataFrame(rows)
    
    # Calculate statistics
    summary = {}
    
    # All samples stats
    summary['total_samples'] = len(df)
    summary['mean_ms_targets_with_min_reads'] = df['ms_targets_with_min_reads'].mean()
    summary['median_ms_targets_with_min_reads'] = df['ms_targets_with_min_reads'].median()
    
    if args.allc_files:
        summary['mean_cpg_count'] = df['cpg_count'].mean()
        summary['median_cpg_count'] = df['cpg_count'].median()
        summary['mean_methylation_rate'] = df['methylation_rate'].mean()
        summary['median_methylation_rate'] = df['methylation_rate'].median()
    
    # Filtered samples stats
    filtered_df = df[df['pass']]
    summary['passing_samples'] = len(filtered_df)
    
    if len(filtered_df) > 0:
        summary['filtered_mean_ms_targets_with_min_reads'] = filtered_df['ms_targets_with_min_reads'].mean()
        summary['filtered_median_ms_targets_with_min_reads'] = filtered_df['ms_targets_with_min_reads'].median()
        
        if args.allc_files:
            summary['filtered_mean_cpg_count'] = filtered_df['cpg_count'].mean()
            summary['filtered_median_cpg_count'] = filtered_df['cpg_count'].median()
            summary['filtered_mean_methylation_rate'] = filtered_df['methylation_rate'].mean()
            summary['filtered_median_methylation_rate'] = filtered_df['methylation_rate'].median()
    
    # Add filter criteria to summary
    summary['min_reads_per_target'] = args.min_reads_per_target
    summary['min_targets_per_sample'] = args.min_targets_per_sample
    summary['min_cpgs_per_sample'] = args.min_cpgs_per_sample
    
    # Write main output TSV
    output_file = os.path.join(args.output_dir, 'sample_stats.tsv')
    output_df = df[['sample_id', 'ms_targets_with_min_reads', 'cpg_count', 'pass']]
    output_df = output_df.sort_values('sample_id')
    output_df.to_csv(output_file, sep='\t', index=False)
    
    # Write summary file
    summary_file = os.path.join(args.output_dir, 'summary_stats.txt')
    with open(summary_file, 'w') as f:
        # Filter parameters
        f.write("Filter Parameters:\n")
        f.write(f"Minimum reads per target: {args.min_reads_per_target}\n")
        f.write(f"Minimum targets per sample: {args.min_targets_per_sample}\n")
        f.write(f"Minimum CpGs per sample: {args.min_cpgs_per_sample}\n\n")
        
        # Sample counts
        f.write("Sample Counts:\n")
        f.write(f"Total samples: {summary['total_samples']}\n")
        f.write(f"Passing samples: {summary['passing_samples']}\n")
        f.write(f"Failed samples: {summary['total_samples'] - summary['passing_samples']}\n\n")
        
        # All samples stats
        f.write("All Samples Statistics:\n")
        f.write(f"Mean microsatellite targets (with min reads): {summary['mean_ms_targets_with_min_reads']:.2f}\n")
        f.write(f"Median microsatellite targets (with min reads): {summary['median_ms_targets_with_min_reads']:.2f}\n")
        
        if args.allc_files:
            f.write(f"Mean CpG count: {summary['mean_cpg_count']:.2f}\n")
            f.write(f"Median CpG count: {summary['median_cpg_count']:.2f}\n")
            f.write(f"Mean methylation rate: {summary['mean_methylation_rate']:.4f}\n")
            f.write(f"Median methylation rate: {summary['median_methylation_rate']:.4f}\n\n")
        
        # Filtered samples stats
        if len(filtered_df) > 0:
            f.write("Filtered Samples Statistics:\n")
            f.write(f"Mean microsatellite targets (with min reads): {summary['filtered_mean_ms_targets_with_min_reads']:.2f}\n")
            f.write(f"Median microsatellite targets (with min reads): {summary['filtered_median_ms_targets_with_min_reads']:.2f}\n")
            
            if args.allc_files:
                f.write(f"Mean CpG count: {summary['filtered_mean_cpg_count']:.2f}\n")
                f.write(f"Median CpG count: {summary['filtered_median_cpg_count']:.2f}\n")
                f.write(f"Mean methylation rate: {summary['filtered_mean_methylation_rate']:.4f}\n")
                f.write(f"Median methylation rate: {summary['filtered_median_methylation_rate']:.4f}\n")
    
    # Generate plots if not disabled
    if not args.no_plots:
        figures_dir = create_plots(df, args.output_dir)
        print(f"Plots generated in {figures_dir}")
    
    print(f"Summary generation completed successfully! {summary['passing_samples']}/{summary['total_samples']} samples passed filtering criteria.")
    print(f"Results written to {output_file} and {summary_file}")

if __name__ == "__main__":
    main() 