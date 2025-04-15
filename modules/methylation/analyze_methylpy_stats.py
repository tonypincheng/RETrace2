#!/usr/bin/env python3

import re
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import argparse
import gzip

def extract_stats_from_logs(log_file_path):
    """Extract statistics from methylpy log file."""
    stats = []
    
    with open(log_file_path, 'r') as file:
        log_text = file.read()
    
    # Split log into sessions for each sample
    file_sessions = re.split(r'\r?\n(?=Processing file:)', log_text)
    
    # Regular expressions for extracting stats
    file_name_regex = r"Processing file:\s*(.+)"
    paired_reads_regex = r"There are (\d+) total input read pairs"
    alignment_rate_regex = r"There are \d+ uniquely mapping read pairs, ([\d\.]+) percent remaining"
    
    for session in file_sessions:
        file_name_match = re.search(file_name_regex, session)
        paired_reads_match = re.search(paired_reads_regex, session)
        alignment_rate_match = re.search(alignment_rate_regex, session)
        
        if file_name_match and paired_reads_match and alignment_rate_match:
            stats.append({
                "sample_name": file_name_match.group(1).strip(),
                "Number of Paired Reads": int(paired_reads_match.group(1)),
                "Alignment Rate": float(alignment_rate_match.group(1))
            })
    
    return pd.DataFrame(stats)

def count_cpgs_from_tsv(tsv_file):
    """Count number of CpGs from methylpy output TSV file."""
    try:
        with gzip.open(tsv_file, 'rt') as f:
            # Skip header
            next(f)
            # Count only lines where context type is CG
            cpg_count = 0
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 4:  # Ensure we have enough fields
                    ctype = fields[3]  # Context type is in the 4th column
                    if ctype.startswith('CG'):
                        cpg_count += 1
            return cpg_count
    except Exception as e:
        print(f"Error reading {tsv_file}: {e}")
        return 0

def generate_summary_stats(df):
    """Generate summary statistics from the dataframe."""
    summary = {
        "Total Reads (millions)": sum(df["Number of Paired Reads"]) / 1e6,
        "Median Reads": df["Number of Paired Reads"].median(),
        "Mean Reads": df["Number of Paired Reads"].mean(),
        "Median Alignment Rate": df["Alignment Rate"].median(),
        "Mean Alignment Rate": df["Alignment Rate"].mean(),
        "Total CpGs": sum(df["Number of CpGs"]),
        "Median CpGs": df["Number of CpGs"].median(),
        "Mean CpGs": df["Number of CpGs"].mean()
    }
    return summary

def create_plots(df, output_dir):
    """Create plots of read counts, alignment rates, and CpG counts."""
    # Create figure
    fig, (ax1, ax2, ax3) = plt.subplots(
        figsize=(26, 15),
        dpi=150,
        nrows=3,
        gridspec_kw={'height_ratios': [2, 1, 1]}
    )
    
    # Plot 1: Number of Reads
    ax1.bar(df.sample_name, df["Number of Paired Reads"]/1e6, color='tab:red')
    ax1.set_xticks(range(len(df.sample_name)))
    ax1.set_xticklabels(df.sample_name, rotation=55, ha='right', fontsize=12)
    ax1.set_ylabel("Number of Reads (Millions)", fontsize=20)
    ax1.tick_params(axis='x', labelbottom=False)
    
    # Plot 2: Alignment Rate
    ax2.bar(df.sample_name, df['Alignment Rate'], color='tab:purple')
    ax2.set_xticks(range(len(df.sample_name)))
    ax2.set_xticklabels(df.sample_name, rotation=55, ha='right', fontsize=12)
    ax2.set_ylabel("Alignment Rate (%)", fontsize=20)
    ax2.set_ylim(0, 100)
    ax2.tick_params(axis='x', labelbottom=False)
    
    # Plot 3: Number of CpGs
    ax3.bar(df.sample_name, df['Number of CpGs']/1e6, color='tab:blue')
    ax3.set_xticks(range(len(df.sample_name)))
    ax3.set_xticklabels(df.sample_name, rotation=55, ha='right', fontsize=12)
    ax3.set_xlabel("Samples", fontsize=20)
    ax3.set_ylabel("Number of CpGs (Millions)", fontsize=20)
    
    # Style adjustments
    for ax in (ax1, ax2, ax3):
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.tick_params(axis='y', which='major', labelsize=15)
    
    plt.tight_layout()
    plt.savefig(Path(output_dir) / "methylpy_summary_plot.pdf")
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Analyze methylpy log files and CpG TSV files to generate summary statistics')
    parser.add_argument('--log', required=True, help='Path to methylpy log file')
    parser.add_argument('--tsv-dir', required=True, help='Directory containing CpG TSV files')
    parser.add_argument('--output-dir', required=True, help='Directory to save output files')
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Extract stats from log
    df = extract_stats_from_logs(args.log)
    
    # Add CpG counts from TSV files
    df['Number of CpGs'] = df['sample_name'].apply(
        lambda x: count_cpgs_from_tsv(Path(args.tsv_dir) / f"{x}_CpG.tsv.gz")
    )
    
    # Generate summary statistics
    summary = generate_summary_stats(df)
    
    # Save detailed stats to CSV
    df.to_csv(output_dir / "detailed_stats.csv", index=False)
    
    # Save summary stats to text file
    with open(output_dir / "summary_stats.txt", 'w') as f:
        for key, value in summary.items():
            f.write(f"{key}: {value}\n")
    
    # Create and save plots
    create_plots(df, output_dir)

if __name__ == "__main__":
    main() 