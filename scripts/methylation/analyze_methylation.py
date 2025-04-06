#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description='Analyze methylation patterns and infer cell types')
    parser.add_argument('--input', required=True, help='Input BAM file')
    parser.add_argument('--genome', required=True, help='Reference genome FASTA')
    parser.add_argument('--output', required=True, help='Output BED file')
    parser.add_argument('--stats', required=True, help='Output statistics file')
    return parser.parse_args()

def call_methylation(bam_file, genome_file, output_bed):
    """Call methylation using methylpl"""
    cmd = [
        'methylpl',
        'call',
        '--bam', bam_file,
        '--genome', genome_file,
        '--output', output_bed
    ]
    
    try:
        subprocess.run(cmd, check=True)
        print(f"Successfully generated {output_bed}")
    except subprocess.CalledProcessError as e:
        print(f"Error running methylpl: {e}", file=sys.stderr)
        sys.exit(1)

def analyze_methylation(bed_file):
    """Analyze methylation patterns"""
    # Read BED file
    df = pd.read_csv(bed_file, sep='\t', header=None,
                    names=['chrom', 'start', 'end', 'name', 'score', 'strand',
                           'meth_level', 'coverage', 'p_value'])
    
    # Calculate statistics
    stats = {
        'mean_methylation': df['meth_level'].mean(),
        'median_methylation': df['meth_level'].median(),
        'coverage_mean': df['coverage'].mean(),
        'coverage_median': df['coverage'].median(),
        'total_sites': len(df)
    }
    
    return df, stats

def infer_cell_types(meth_data, n_clusters=3):
    """Infer cell types using clustering"""
    # Prepare data for clustering
    X = meth_data[['meth_level', 'coverage']].values
    
    # Perform PCA
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X)
    
    # Cluster using K-means
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    clusters = kmeans.fit_predict(X_pca)
    
    return X_pca, clusters

def plot_results(X_pca, clusters, output_plot):
    """Plot clustering results"""
    plt.figure(figsize=(10, 8))
    sns.scatterplot(x=X_pca[:, 0], y=X_pca[:, 1], hue=clusters, palette='viridis')
    plt.title('Methylation Pattern Clustering')
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.savefig(output_plot)
    plt.close()

def write_stats(stats, stats_file):
    """Write statistics to file"""
    with open(stats_file, 'w') as f:
        for key, value in stats.items():
            f.write(f"{key}\t{value}\n")

def main():
    args = parse_args()
    
    # Call methylation
    print("Calling methylation...")
    call_methylation(args.input, args.genome, args.output)
    
    # Analyze methylation patterns
    print("Analyzing methylation patterns...")
    meth_data, stats = analyze_methylation(args.output)
    
    # Infer cell types
    print("Inferring cell types...")
    X_pca, clusters = infer_cell_types(meth_data)
    
    # Plot results
    print("Generating plots...")
    plot_results(X_pca, clusters, f"{args.output}.plot.pdf")
    
    # Write statistics
    print("Writing statistics...")
    write_stats(stats, args.stats)
    
    print("Done!")

if __name__ == '__main__':
    main() 