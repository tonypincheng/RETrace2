#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
import sys
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description='Perform bootstrap analysis on phylogenetic trees')
    parser.add_argument('--tree', required=True, help='Input tree file (Newick format)')
    parser.add_argument('--vcf', required=True, help='Input VCF file')
    parser.add_argument('--iterations', type=int, default=100, help='Number of bootstrap iterations')
    parser.add_argument('--output', required=True, help='Output bootstrap trees file')
    parser.add_argument('--support', required=True, help='Output support values file')
    parser.add_argument('--plot', required=True, help='Output plot file')
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

def bootstrap_sample(genotypes):
    """Generate bootstrap sample of genotypes"""
    n_sites = len(genotypes)
    indices = np.random.choice(n_sites, size=n_sites, replace=True)
    return genotypes.iloc[indices]

def calculate_support(original_tree, bootstrap_trees):
    """Calculate support values for each branch"""
    support = {}
    for clade in original_tree.get_nonterminals():
        clade_str = str(sorted([leaf.name for leaf in clade.get_terminals()]))
        support[clade_str] = 0
    
    for tree in bootstrap_trees:
        for clade in tree.get_nonterminals():
            clade_str = str(sorted([leaf.name for leaf in clade.get_terminals()]))
            if clade_str in support:
                support[clade_str] += 1
    
    # Convert to percentages
    for clade in support:
        support[clade] = (support[clade] / len(bootstrap_trees)) * 100
    
    return support

def write_support(support, support_file):
    """Write support values to file"""
    with open(support_file, 'w') as f:
        for clade, value in support.items():
            f.write(f"{clade}\t{value:.2f}\n")

def main():
    args = parse_args()
    
    # Read original tree
    print("Reading original tree...")
    original_tree = Phylo.read(args.tree, 'newick')
    
    # Read VCF
    print("Reading VCF file...")
    genotypes = read_vcf(args.vcf)
    
    # Perform bootstrap
    print(f"Performing {args.iterations} bootstrap iterations...")
    bootstrap_trees = []
    
    for i in range(args.iterations):
        print(f"Iteration {i+1}/{args.iterations}")
        # Generate bootstrap sample
        boot_genotypes = bootstrap_sample(genotypes)
        
        # Calculate distance matrix
        distances = np.zeros((len(genotypes.columns), len(genotypes.columns)))
        for i in range(len(genotypes.columns)):
            for j in range(i+1, len(genotypes.columns)):
                dist = sum(boot_genotypes.iloc[:,i] != boot_genotypes.iloc[:,j])
                distances[i,j] = distances[j,i] = dist
        
        # Build tree
        matrix = DistanceMatrix(names=genotypes.columns, matrix=distances)
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(matrix)
        bootstrap_trees.append(tree)
    
    # Calculate support values
    print("Calculating support values...")
    support = calculate_support(original_tree, bootstrap_trees)
    
    # Write results
    print("Writing results...")
    with open(args.output, 'w') as f:
        for tree in bootstrap_trees:
            Phylo.write(tree, f, 'newick')
    
    write_support(support, args.support)
    
    print("Done!")

if __name__ == '__main__':
    main() 