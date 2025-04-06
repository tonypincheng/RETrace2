#!/usr/bin/env python3

import argparse
import numpy as np
from Bio import Phylo
import sys
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description='Evaluate phylogenetic tree accuracy')
    parser.add_argument('--tree', required=True, help='Input tree file (Newick format)')
    parser.add_argument('--ground_truth', required=True, help='Ground truth tree file (Newick format)')
    parser.add_argument('--output', required=True, help='Output results file')
    parser.add_argument('--plot', required=True, help='Output plot file')
    return parser.parse_args()

def get_clades(tree):
    """Extract all clades from a tree"""
    clades = set()
    for clade in tree.get_nonterminals():
        clade_str = str(sorted([leaf.name for leaf in clade.get_terminals()]))
        clades.add(clade_str)
    return clades

def calculate_metrics(test_tree, true_tree):
    """Calculate accuracy metrics"""
    test_clades = get_clades(test_tree)
    true_clades = get_clades(true_tree)
    
    # Calculate true positives, false positives, and false negatives
    tp = len(test_clades.intersection(true_clades))
    fp = len(test_clades - true_clades)
    fn = len(true_clades - test_clades)
    
    # Calculate metrics
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    
    return {
        'precision': precision,
        'recall': recall,
        'f1': f1,
        'true_positives': tp,
        'false_positives': fp,
        'false_negatives': fn
    }

def write_results(metrics, output_file):
    """Write evaluation results to file"""
    with open(output_file, 'w') as f:
        f.write("Metric\tValue\n")
        for metric, value in metrics.items():
            f.write(f"{metric}\t{value}\n")

def main():
    args = parse_args()
    
    # Read trees
    print("Reading trees...")
    test_tree = Phylo.read(args.tree, 'newick')
    true_tree = Phylo.read(args.ground_truth, 'newick')
    
    # Calculate metrics
    print("Calculating metrics...")
    metrics = calculate_metrics(test_tree, true_tree)
    
    # Write results
    print("Writing results...")
    write_results(metrics, args.output)
    
    print("Done!")

if __name__ == '__main__':
    main() 