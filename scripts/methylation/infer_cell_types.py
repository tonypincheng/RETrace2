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
    parser = argparse.ArgumentParser(description='Infer cell types from methylation patterns')
    parser.add_argument('--input', required=True, help='Input BED files')
    parser.add_argument('--output', required=True, help='Output cell type predictions file')
    parser.add_argument('--plot', required=True, help='Output plot file')
    parser.add_argument('--num_types', type=int, default=3, help='Number of cell types to infer')
    return parser.parse_args()

def read_methylation(bed_files):
    """Read methylation data from BED files"""
    # Combine all methylation data
    all_data = []
    sample_names = []
    
    for bed_file in bed_files.split():
        sample_name = Path(bed_file).stem
        sample_names.append(sample_name)
        
        df = pd.read_csv(bed_file, sep='\t', header=None,
                        names=['chrom', 'start', 'end', 'name', 'score', 'strand',
                               'meth_level', 'coverage', 'p_value'])
        
        # Calculate average methylation for this sample
        avg_meth = df['meth_level'].mean()
        coverage = df['coverage'].mean()
        total_sites = len(df)
        density = df.groupby('chrom').size().mean()
        
        all_data.append({
            'sample': sample_name,
            'avg_methylation': avg_meth,
            'coverage': coverage,
            'total_sites': total_sites,
            'density': density
        })
    
    return pd.DataFrame(all_data)

def infer_cell_types(meth_data, n_clusters):
    """Infer cell types using clustering"""
    # Feature matrix for clustering
    X = meth_data[['avg_methylation', 'coverage', 'density']].values
    
    # Normalize features
    X = (X - X.mean(axis=0)) / X.std(axis=0)
    
    # Perform PCA
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X)
    
    # Cluster using K-means
    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    clusters = kmeans.fit_predict(X_pca)
    
    # Add cluster assignments to data
    meth_data['cell_type'] = clusters
    
    return meth_data, X_pca, clusters

def plot_results(X_pca, clusters, samples, output_plot):
    """Plot clustering results"""
    plt.figure(figsize=(10, 8))
    scatter = plt.scatter(X_pca[:, 0], X_pca[:, 1], c=clusters, cmap='viridis', s=100)
    
    # Add sample labels
    for i, sample in enumerate(samples):
        plt.annotate(sample, (X_pca[i, 0], X_pca[i, 1]), fontsize=8)
    
    plt.colorbar(scatter, label='Cell Type')
    plt.title('Cell Type Inference from Methylation Patterns')
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.tight_layout()
    plt.savefig(output_plot, dpi=300)
    plt.close()

def main():
    args = parse_args()
    
    # Read methylation data
    print("Reading methylation data...")
    meth_data = read_methylation(args.input)
    
    # Infer cell types
    print(f"Inferring {args.num_types} cell types...")
    results, X_pca, clusters = infer_cell_types(meth_data, args.num_types)
    
    # Write results
    print("Writing results...")
    results.to_csv(args.output, sep='\t', index=False)
    
    # Plot results
    print("Generating plot...")
    plot_results(X_pca, clusters, results['sample'], args.plot)
    
    print("Done!")

if __name__ == '__main__':
    main() 