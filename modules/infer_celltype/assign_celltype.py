#!/usr/bin/env python3
"""
Cell Type Assignment Script for Methylation Data

This script assigns cell types to samples based on pairwise dissimilarity matrices:
1. Reads pairwise dissimilarity matrix and optional sites matrix
2. Calculates z-score transformed matrix for each sample
3. Assigns cell types based on minimum z-score (best match) and threshold
4. Generates visualization plots and outputs results

The sites matrix (optional) contains the number of shared CpG sites between each 
sample and reference cell type, providing context for the reliability of 
dissimilarity calculations.

Usage:
    python assign_celltype.py --pd_matrix matrix.csv --output_dir results/
    python assign_celltype.py --pd_matrix matrix.csv --sites_matrix sites.csv --output_dir results/
"""

import argparse
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from pathlib import Path

def setup_matplotlib():
    """Set up matplotlib parameters for high-quality plots"""
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams['font.size'] = 10
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['axes.titlesize'] = 14
    plt.rcParams['legend.fontsize'] = 10
    plt.rcParams['xtick.labelsize'] = 8
    plt.rcParams['ytick.labelsize'] = 8

def read_matrix(filepath):
    """
    Read matrix from file (supports CSV, TSV, Excel).
    
    Args:
        filepath (str): Path to the matrix file
        
    Returns:
        pd.DataFrame: Matrix with proper index/columns
    """
    try:
        # Try to read as CSV first
        if filepath.endswith('.csv'):
            df = pd.read_csv(filepath, index_col=0)
        elif filepath.endswith('.tsv') or filepath.endswith('.txt'):
            df = pd.read_csv(filepath, sep='\t', index_col=0)
        elif filepath.endswith('.xlsx'):
            df = pd.read_excel(filepath, index_col=0)
        else:
            # Try to auto-detect format
            df = pd.read_csv(filepath, index_col=0)
        
        print(f"Successfully loaded matrix with shape: {df.shape}")
        print(f"Sample names: {list(df.index[:5])}...")
        print(f"Reference types: {list(df.columns[:5])}...")
        
        return df
    
    except Exception as e:
        print(f"Error reading matrix file: {e}")
        sys.exit(1)

def calculate_adaptive_figsize(df, base_width=6, base_height=10, min_width=4, min_height=6, max_width=10, max_height=15):
    """
    Calculate adaptive figure size based on dataframe dimensions.
    
    Args:
        df (pd.DataFrame): Data matrix
        base_width (int): Base width for scaling
        base_height (int): Base height for scaling
        min_width (int): Minimum width
        min_height (int): Minimum height
        max_width (int): Maximum width
        max_height (int): Maximum height
        
    Returns:
        tuple: (width, height) for figure size
    """
    n_samples = len(df.index)
    n_references = len(df.columns)
    
    # Scale based on number of samples and references
    # Use log scaling to prevent extremely large plots
    width_scale = min(n_references / 10, 2.0)  # Cap at 2x scaling
    height_scale = min(n_samples / 10, 1.5)   # Cap at 1.5x scaling
    
    # Calculate adaptive size
    width = base_width * width_scale
    height = base_height * height_scale
    
    # Apply bounds
    width = max(min_width, min(width, max_width))
    height = max(min_height, min(height, max_height))
    
    return (width, height)

def plot_sites_matrix(df_sites, output_dir):
    """
    Generate heatmap for sites matrix showing shared CpG sites between samples and references.
    
    The sites matrix contains the number of shared CpG sites used for calculating 
    pairwise dissimilarity between each sample and reference cell type. Higher values 
    indicate more shared sites available for comparison, which typically leads to 
    more reliable dissimilarity calculations.
    """
    figsize = calculate_adaptive_figsize(df_sites)
    plt.figure(figsize=figsize)
    
    # Create heatmap
    sns.heatmap(df_sites, 
                cmap='RdBu_r',
                cbar_kws={'label': 'Number of Shared CpG Sites'},
                xticklabels=True,
                yticklabels=True,
                linewidths=0.5)
    
    plt.title('Sites Matrix - Shared CpG Sites Count', fontsize=16, pad=20)
    plt.xlabel('Reference Cell Types', fontsize=12)
    plt.ylabel('Sample IDs', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    
    # Save plot as PDF
    output_path = os.path.join(output_dir, 'sites_matrix.pdf')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Sites matrix plot saved to: {output_path}")
    plt.close()

def plot_raw_matrix(df, output_dir):
    """
    Generate heatmap for raw pairwise dissimilarity matrix.
    
    Args:
        df (pd.DataFrame): Raw PD matrix
        output_dir (str): Output directory for saving plots
    """
    figsize = calculate_adaptive_figsize(df)
    plt.figure(figsize=figsize)
    
    # Create heatmap
    sns.heatmap(df, 
                cmap='RdBu_r',
                cbar_kws={'label': 'Pairwise Dissimilarity'},
                xticklabels=True,
                yticklabels=True,
                linewidths=0.5)
    
    plt.title('Raw Pairwise Dissimilarity Matrix', fontsize=16, pad=20)
    plt.xlabel('Reference Cell Types', fontsize=12)
    plt.ylabel('Sample IDs', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    
    # Save plot as PDF
    output_path = os.path.join(output_dir, 'raw_pd_matrix.pdf')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Raw PD matrix plot saved to: {output_path}")
    plt.close()

def plot_zscore_matrix(df_zscore, output_dir):
    """
    Generate heatmap for z-score transformed pairwise dissimilarity matrix.
    
    Args:
        df_zscore (pd.DataFrame): Z-score transformed PD matrix
        output_dir (str): Output directory for saving plots
        threshold (float): Z-score threshold for assignment
    """
    figsize = calculate_adaptive_figsize(df_zscore)
    plt.figure(figsize=figsize)
    
    # Create heatmap with diverging colormap centered at 0
    sns.heatmap(df_zscore,
                cmap='RdBu_r',
                center=0,
                cbar_kws={'label': 'Z-score'},
                xticklabels=True,
                yticklabels=True,
                linewidths=0.5)
    
    plt.title('Pairwise Dissimilarity Matrix (Z-score of each row)', 
              fontsize=16, pad=20)
    plt.xlabel('Reference Cell Types', fontsize=12)
    plt.ylabel('Sample IDs', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    
    # Save plot as PDF
    output_path = os.path.join(output_dir, 'zscore_pd_matrix.pdf')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Z-score PD matrix plot saved to: {output_path}")
    plt.close()

def calculate_zscore_matrix(df):
    """
    Calculate z-score transformation for each sample (row-wise).
    
    Args:
        df (pd.DataFrame): Raw PD matrix
        
    Returns:
        pd.DataFrame: Z-score transformed matrix
    """
    # Calculate z-score for each row (sample)
    df_zscore = df.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
    
    print(f"Z-score matrix statistics:")
    print(f"Mean: {df_zscore.values.mean():.4f}")
    print(f"Std: {df_zscore.values.std():.4f}")
    print(f"Min: {df_zscore.values.min():.4f}")
    print(f"Max: {df_zscore.values.max():.4f}")
    
    return df_zscore

def create_assignment_table(df_zscore, threshold=-1.2):
    """
    Create assignment table based on z-score threshold.
    
    Args:
        df_zscore (pd.DataFrame): Z-score transformed matrix
        threshold (float): Z-score threshold for assignment
        
    Returns:
        pd.DataFrame: Assignment table
    """
    # Find the minimum z-score (best match) for each sample
    best_matches = df_zscore.idxmin(axis=1)
    min_zscores = df_zscore.min(axis=1)
    
    # Calculate second best matches and z-scores
    second_best_matches = df_zscore.apply(
        lambda row: row.nsmallest(2).index[1], axis=1
    )
    second_best_zscores = df_zscore.apply(
        lambda row: row.nsmallest(2).iloc[1], axis=1
    )
    
    # Check if passes threshold
    passes_threshold = min_zscores < threshold
    
    # Create assignment table
    assignment_table = pd.DataFrame({
        'Sample_ID': df_zscore.index,
        'Best_Match': best_matches,
        'Min_Zscore': min_zscores,
        'Second_Best_Match': second_best_matches,
        'Second_Best_Zscore': second_best_zscores,
        'Passes_Threshold': passes_threshold,
        'Assignment': np.where(passes_threshold, best_matches, 'Unassigned')
    })
    
    return assignment_table

def print_assignment_summary(assignment_table, threshold):
    """
    Print summary statistics for assignments.
    
    Args:
        assignment_table (pd.DataFrame): Assignment table
        threshold (float): Z-score threshold used
    """
    total_samples = len(assignment_table)
    assigned_samples = assignment_table['Passes_Threshold'].sum()
    unassigned_samples = total_samples - assigned_samples
    
    print(f"\nAssignment Summary (z-score threshold: {threshold}):")
    print(f"Total samples: {total_samples}")
    print(f"Assigned samples: {assigned_samples} ({assigned_samples/total_samples*100:.1f}%)")
    print(f"Unassigned samples: {unassigned_samples} ({unassigned_samples/total_samples*100:.1f}%)")
    
    if assigned_samples > 0:
        print(f"\nAssignment distribution:")
        assignment_counts = assignment_table[assignment_table['Passes_Threshold']]['Best_Match'].value_counts()
        for cell_type, count in assignment_counts.items():
            print(f"  {cell_type}: {count} samples")

def main():
    parser = argparse.ArgumentParser(
        description='Assign cell types to samples based on pairwise dissimilarity matrices'
    )
    parser.add_argument('--pd_matrix', required=True,
                        help='Path to pairwise dissimilarity matrix file')
    parser.add_argument('--sites_matrix', required=False,
                        help='Path to shared CpG sites matrix file (optional) - contains counts of shared sites between samples and references')
    parser.add_argument('--output_dir', required=True,
                        help='Output directory for results')
    parser.add_argument('--threshold', '-t', type=float, default=-1.2,
                        help='Z-score threshold for assignment (default: -1.2)')
    parser.add_argument('--figsize', nargs=2, type=int, default=[12, 10],
                        help='Figure size for plots (width height, default: 12 10)')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.pd_matrix):
        print(f"Error: Pairwise dissimilarity matrix file {args.pd_matrix} does not exist")
        sys.exit(1)
    
    if args.sites_matrix and not os.path.exists(args.sites_matrix):
        print(f"Error: Sites matrix file {args.sites_matrix} does not exist")
        sys.exit(1)
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Setup matplotlib
    setup_matplotlib()
    
    print(f"Reading pairwise dissimilarity matrix from: {args.pd_matrix}")
    if args.sites_matrix:
        print(f"Reading shared CpG sites matrix from: {args.sites_matrix}")
    else:
        print("No sites matrix provided - proceeding with PD matrix only")
    print(f"Output directory: {args.output_dir}")
    print(f"Z-score threshold: {args.threshold}")
    
    # Read the matrices
    pd_matrix = read_matrix(args.pd_matrix)
    sites_matrix = None
    if args.sites_matrix:
        sites_matrix = read_matrix(args.sites_matrix)
        print(f"Shared CpG sites matrix loaded (shape: {sites_matrix.shape})")
        print("Proceeding with cell type assignment using PD matrix and sites matrix for visualization")
    
    # Generate sites matrix plot if provided
    if sites_matrix is not None:
        print("\nGenerating shared CpG sites matrix plot...")
        plot_sites_matrix(sites_matrix, args.output_dir)
    
    # Generate raw matrix plot
    print("\nGenerating raw PD matrix plot...")
    plot_raw_matrix(pd_matrix, args.output_dir)
    
    # Calculate z-score transformation
    print("\nCalculating z-score transformation...")
    df_zscore = calculate_zscore_matrix(pd_matrix)
    
    # Generate z-score matrix plot
    print("\nGenerating z-score PD matrix plot...")
    plot_zscore_matrix(df_zscore, args.output_dir)
    
    # Create assignment table
    print(f"\nCreating cell type assignments...")
    assignment_table = create_assignment_table(df_zscore, args.threshold)
    
    # Save assignment table
    assignments_path = os.path.join(args.output_dir, 'celltype_assignments.tsv')
    assignment_table.to_csv(assignments_path, sep='\t', index=False)
    print(f"Cell type assignments saved to: {assignments_path}")
    
    # Save z-score matrix for reference
    zscore_path = os.path.join(args.output_dir, 'zscore_matrix.csv')
    df_zscore.to_csv(zscore_path)
    print(f"Z-score matrix saved to: {zscore_path}")
    
    # Print summary
    print_assignment_summary(assignment_table, args.threshold)
    
    print(f"\nAnalysis complete! All outputs saved to: {args.output_dir}")

if __name__ == "__main__":
    main() 