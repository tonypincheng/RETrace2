#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
import sys
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description='Build phylogenetic tree from HipSTR VCF files')
    parser.add_argument('--vcf', required=True, help='Input VCF file')
    parser.add_argument('--matrix', required=True, help='Output distance matrix file')
    parser.add_argument('--tree', required=True, help='Output tree file (Newick format)')
    parser.add_argument('--stats', required=True, help='Output statistics file')
    return parser.parse_args()

def read_vcf(vcf_file):
    """Read VCF file and extract genotype information"""
    genotypes = {}
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    samples = line.strip().split('\t')[9:]
                    for sample in samples:
                        genotypes[sample] = []
                continue
            
            fields = line.strip().split('\t')
            for i, sample in enumerate(samples):
                gt = fields[9+i].split(':')[0]
                genotypes[sample].append(gt)
    
    return pd.DataFrame(genotypes)

def calculate_distance_matrix(genotypes):
    """Calculate distance matrix between samples"""
    n_samples = len(genotypes.columns)
    distances = np.zeros((n_samples, n_samples))
    
    for i in range(n_samples):
        for j in range(i+1, n_samples):
            # Calculate Hamming distance between samples
            dist = sum(genotypes.iloc[:,i] != genotypes.iloc[:,j])
            distances[i,j] = distances[j,i] = dist
    
    return distances

def build_tree(distances, samples):
    """Build phylogenetic tree using neighbor-joining"""
    # Create distance matrix object
    matrix = DistanceMatrix(names=samples, matrix=distances)
    
    # Build tree
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(matrix)
    
    return tree

def write_stats(tree, stats_file):
    """Write tree statistics to file"""
    with open(stats_file, 'w') as f:
        f.write(f"Number of leaves: {len(tree.get_terminals())}\n")
        f.write(f"Tree height: {tree.total_branch_length()}\n")

def main():
    args = parse_args()
    
    # Read VCF and extract genotypes
    print("Reading VCF file...")
    genotypes = read_vcf(args.vcf)
    
    # Calculate distance matrix
    print("Calculating distance matrix...")
    distances = calculate_distance_matrix(genotypes)
    
    # Save distance matrix
    print("Saving distance matrix...")
    np.savetxt(args.matrix, distances, delimiter='\t')
    
    # Build tree
    print("Building phylogenetic tree...")
    tree = build_tree(distances, genotypes.columns)
    
    # Save tree
    print("Saving tree...")
    Phylo.write(tree, args.tree, 'newick')
    
    # Write statistics
    print("Writing statistics...")
    write_stats(tree, args.stats)
    
    print("Done!")

if __name__ == '__main__':
    main() 