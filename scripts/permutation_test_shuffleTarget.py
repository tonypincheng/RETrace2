#!/usr/bin/env python3

import os
import sys
import subprocess
import pickle
import glob
import json
import random
import time
import pandas as pd
import numpy as np
from Bio import Phylo
from multiprocessing import Pool
import matplotlib.pyplot as plt
import argparse

# Add the parent directory to the Python path to find the modules package
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import parsing functions from modules/hipstr/parse_vcf.py
from modules.hipstr.parse_vcf import parseVCF, load_target_bed

# ----------------------------
# Set up argument parser
# ----------------------------
def parse_arguments():
    parser = argparse.ArgumentParser(description="Perform permutation test by shuffling microsatellite (target) genotypes from a VCF file")
    
    # Required arguments
    parser.add_argument("--samplesheet", required=True, help="Path to samplesheet file")
    parser.add_argument("--target_bed", required=True, help="Path to target BED file")
    parser.add_argument("--input_vcf", required=True, help="Path to input VCF file")
    parser.add_argument("--sample_list", required=True, help="Path to sample list file. This should be generated with the vcf file in the hipstr_output directory")
    parser.add_argument("--observed_pkl", required=True, help="Path to observed pickle file. This should be generated with the vcf file in the hipstr_output directory")
    
    # Optional arguments with defaults (Use the same setting as how you processed your observed data)
    parser.add_argument("--output_prefix", default="permutation_test", help="Prefix for output files (default: permutation_test)")
    parser.add_argument("--build_phylo_path", default="modules/phylo/build_phylo.py", 
                       help="Path to build_phylo.py script (default: modules/phylo/build_phylo.py)")
    parser.add_argument("--permuted_vcfs_dir", default="permuted_vcfs", help="Directory for permuted VCFs (default: permuted_vcfs)")
    parser.add_argument("--permuted_pkls_dir", default="permuted_pkls", help="Directory for permuted pickle files (default: permuted_pkls)")
    parser.add_argument("--phylo_output_dir", default="permuted_phylo", help="Directory for phylogeny output (default: permuted_phylo)")
    parser.add_argument("--n_permutations", type=int, default=100, help="Number of permutations (default: 100)")
    parser.add_argument("--min_qual", type=float, default=0.9, help="Minimum quality threshold (default: 0.9)")
    parser.add_argument("--min_reads", type=int, default=10, help="Minimum number of reads (default: 10)")
    parser.add_argument("--max_stutter", type=float, default=1.0, help="Maximum stutter fraction (default: 1.0)")
    parser.add_argument("--dist_metric", default="EqorNot_minComp", 
                       help="Distance metric to use for phylogeny construction (default: EqorNot_minComp)")
    parser.add_argument("--outgroup", default="None", 
                       help="Outgroup to use for rooting the phylogenetic tree (default: None)")
    parser.add_argument("--n_processes", type=int, default=4, help="Number of processes for multiprocessing (default: 4)")
    
    return parser.parse_args()

# ----------------------------
# STEP 1: Shuffle VCF functions
# ----------------------------

def extract_gt(row, sample):
    """Extract GT field from sample data (for checking missingness)"""
    format_fields = row['FORMAT'].split(':')
    sample_data = row[sample].split(':')
    gt_idx = format_fields.index('GT')
    return sample_data[gt_idx]

def process_permutation(perm, df, sample_names, header_lines, output_dir):
    """Process a single permutation of the VCF file"""
    permuted_df = df.copy()
    
    for index, row in permuted_df.iterrows():
        # Get full sample info for all samples in this row
        sample_infos = []
        valid_samples = []
        for sample in sample_names:
            gt = extract_gt(row, sample)
            if gt != '.' and gt != './.':  # Only include non-missing genotypes
                sample_infos.append(row[sample])
                valid_samples.append(sample)
        
        # Shuffle full sample info among valid samples
        if len(sample_infos) > 1:  # Only shuffle if there are multiple valid entries
            shuffled_infos = sample_infos.copy()
            random.shuffle(shuffled_infos)
            
            # Assign shuffled info back to valid samples
            for sample, new_info in zip(valid_samples, shuffled_infos):
                permuted_df.at[index, sample] = new_info
    
    # Write permuted VCF
    output_file = os.path.join(output_dir, f"permuted_{perm+1}.vcf")
    with open(output_file, 'w') as f:
        f.writelines(header_lines)
        for _, row in permuted_df.iterrows():
            f.write('\t'.join(row) + '\n')
    
    return f"Generated {output_file}"

# ----------------------------
# STEP 2: Parse VCF functions
# ----------------------------

def process_single_vcf(args_tuple):
    """Process a single VCF file (for multiprocessing)"""
    target_bed, vcf_file, output_dir, min_qual, min_reads, max_stutter = args_tuple
    alleleDict, _ = parseVCF(vcf_file, target_bed, min_qual, min_reads, max_stutter)
    base_name = os.path.basename(vcf_file).replace('.vcf', '.pkl')
    output_file = os.path.join(output_dir, base_name)
    with open(output_file, 'wb') as f:
        pickle.dump(alleleDict, f)
    return f"Saved parsed data to {output_file}"

def read_group_mapping(file_path):
    """
    Read the samplesheet file and extract sample IDs and group information.
    Expects columns named 'sample_id' and 'group'.
    """
    group_dc = {}
    try:
        with open(file_path, 'r') as file:
            header = None
            sample_idx = None
            group_idx = None
            
            for line in file:
                # Skip comments
                if line.startswith('#'):
                    continue
                    
                columns = line.strip().split(',')  # Assuming CSV format
                
                # Process header
                if header is None:
                    header = columns
                    try:
                        sample_idx = header.index('sample_id')
                        group_idx = header.index('group')
                    except ValueError:
                        print("Error: Could not find 'sample_id' and 'group' columns in samplesheet.")
                        print(f"Header found: {header}")
                        return {}
                    continue
                
                # Process data rows
                if len(columns) > max(sample_idx, group_idx):
                    sample_id = columns[sample_idx].strip()
                    group = columns[group_idx].strip()
                    if sample_id and group:
                        group_dc[sample_id] = group
        
        if not group_dc:
            print("Warning: No group mappings found in samplesheet.")
        return group_dc
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        return {}

# ----------------------------
# STEP 3: Process Trees functions
# ----------------------------

def calculate_cophenetic_distance(tree, tip1_name, tip2_name):
    """Calculate cophenetic distance between two tips in the tree"""
    tip1 = [t for t in tree.get_terminals() if t.name == tip1_name]
    tip2 = [t for t in tree.get_terminals() if t.name == tip2_name]
    if not tip1 or not tip2:
        return -1
    tip1, tip2 = tip1[0], tip2[0]
    
    path1 = []
    current = tip1
    while current:
        path1.append(current)
        if current == tree.root:
            break
        for clade in tree.find_clades():
            if current in clade.clades:
                current = clade
                break
    
    path2 = []
    current = tip2
    while current:
        path2.append(current)
        if current == tree.root:
            break
        for clade in tree.find_clades():
            if current in clade.clades:
                current = clade
                break
    
    mrca = None
    for node1 in path1:
        for node2 in path2:
            if node1 == node2:
                mrca = node1
                break
        if mrca:
            break
    
    if not mrca:
        return -1
    
    dist1 = 0
    current = tip1
    while current != mrca:
        dist1 += 1
        for clade in tree.find_clades():
            if current in clade.clades:
                current = clade
                break
    
    dist2 = 0
    current = tip2
    while current != mrca:
        dist2 += 1
        for clade in tree.find_clades():
            if current in clade.clades:
                current = clade
                break
    
    return dist1 + dist2

def find_siblings(tree, tip):
    """Find sibling tips of a given tip"""
    if tip == tree.root:
        return []
    parent = None
    for clade in tree.find_clades():
        if tip in clade.clades:
            parent = clade
            break
    if not parent:
        return []
    return [sibling for sibling in parent.clades if sibling.is_terminal() and sibling != tip]

def count_same_tissue_pairs(tree, tissues, distance_threshold=2):
    """Count pairs of tips from the same tissue that are close in the tree"""
    tip_names = [tip.name for tip in tree.get_terminals()]
    unique_pairs = set()
    unique_comparisons = set()
    
    for tip1_name in tip_names:
        tissue1 = tissues.get(tip1_name)
        if tissue1 is None:
            continue
        tip1 = [t for t in tree.get_terminals() if t.name == tip1_name][0]
        siblings = find_siblings(tree, tip1)
        
        for sibling in siblings:
            tip2_name = sibling.name
            comparison = tuple(sorted([tip1_name, tip2_name]))
            unique_comparisons.add(comparison)
            
            tissue2 = tissues.get(tip2_name)
            if tissue2 is None:
                continue
            
            dist = calculate_cophenetic_distance(tree, tip1_name, tip2_name)
            if dist == -1 or dist > distance_threshold:
                continue
            
            if tissue1 == tissue2:
                pair = tuple(sorted([tip1_name, tip2_name]))
                unique_pairs.add(pair)
    
    pairs = len(unique_pairs)
    comparisons = len(unique_comparisons)
    percentage = (pairs / comparisons * 100) if comparisons > 0 else 0.0
    return percentage, comparisons, pairs

def process_pkl(pkl_file):
    """Process a single pickle file to build a phylogenetic tree and analyze tissue clustering"""
    prefix = os.path.join(args.phylo_output_dir, f"{args.output_prefix}.perm_{os.path.basename(pkl_file).replace('.pkl', '')}")
    
    # Build phylogeny using modules/phylo/build_phylo.py
    cmd_build = [
        "python", args.build_phylo_path,
        "--alleleDict", pkl_file,
        "--sample_list", args.sample_list,
        "--prefix", prefix,
        "--dist_metric", args.dist_metric,
        "--outgroup", args.outgroup
    ]
    subprocess.run(cmd_build, check=True)
    
    # Load tree (original, not bootstrap)
    tree_file = f"{prefix}.newick"
    tree = Phylo.read(tree_file, "newick")
    tissues = read_group_mapping(args.samplesheet)
    percentage, _, _ = count_same_tissue_pairs(tree, tissues)
    return percentage

def process_observed():
    """Process the observed data to build a phylogenetic tree and analyze tissue clustering"""
    prefix = os.path.join(args.phylo_output_dir, f"{args.output_prefix}")
    cmd_build = [
        "python", args.build_phylo_path,
        "--alleleDict", args.observed_pkl,
        "--sample_list", args.sample_list,
        "--prefix", prefix,
        "--dist_metric", args.dist_metric,
        "--outgroup", args.outgroup
    ]
    subprocess.run(cmd_build, check=True)
    
    tree_file = f"{prefix}.newick"
    tree = Phylo.read(tree_file, "newick")
    tissues = read_group_mapping(args.samplesheet)
    percentage, comparisons, pairs = count_same_tissue_pairs(tree, tissues)
    return percentage, comparisons, pairs

# ----------------------------
# Main function to run all steps
# ----------------------------

def run_permutation_test():
    start_time = time.time()
    
    # Ensure all required directories exist
    for directory in [args.permuted_vcfs_dir, args.permuted_pkls_dir, args.phylo_output_dir]:
        os.makedirs(directory, exist_ok=True)
    
    # Step 1: Shuffle VCF
    print("\n" + "="*80)
    print("STEP 1: SHUFFLING VCF FILES")
    print("="*80)
    
    # Read the VCF file
    print(f"Reading input VCF file: {args.input_vcf}")
    with open(args.input_vcf, 'r') as f:
        lines = f.readlines()
    
    # Separate headers and data
    header_lines = [line for line in lines if line.startswith('#')]
    data_lines = [line.strip().split('\t') for line in lines if not line.startswith('#')]
    
    # Extract sample names from the last header line
    sample_names = header_lines[-1].strip().split('\t')[9:]
    print(f"Found {len(sample_names)} samples in VCF")
    
    # Convert data to a DataFrame
    df = pd.DataFrame(data_lines, columns=header_lines[-1].strip().split('\t'))
    print(f"Read {len(df)} variants from VCF")
    
    # Generate permuted VCFs
    print(f"Generating {args.n_permutations} permuted VCFs using {args.n_processes} processes...")
    with Pool(processes=args.n_processes) as pool:
        results = pool.starmap(process_permutation, [(perm, df, sample_names, header_lines, args.permuted_vcfs_dir) for perm in range(args.n_permutations)])
    
    print(f"Completed generating {args.n_permutations} permuted VCFs")
    step1_time = time.time()
    print(f"Step 1 completed in {(step1_time - start_time)/60:.2f} minutes")
    
    # Step 2: Parse VCFs
    print("\n" + "="*80)
    print("STEP 2: PARSING VCF FILES TO EXTRACT ALLELOTYPES")
    print("="*80)
    
    # Import target information
    print(f"Loading target information from {args.target_bed}")
    # Load the target information using the imported function
    target_bed = args.target_bed
    print(f"Target BED file: {target_bed}")
    
    # Find all permuted VCF files
    permuted_vcf_files = glob.glob(os.path.join(args.permuted_vcfs_dir, "permuted_*.vcf"))
    print(f"Found {len(permuted_vcf_files)} permuted VCF files to process.")
    
    # Prepare arguments for each process
    args_list = [(target_bed, vcf_file, args.permuted_pkls_dir, args.min_qual, args.min_reads, args.max_stutter) 
                 for vcf_file in permuted_vcf_files]
    
    # Use multiprocessing Pool
    print(f"Parsing VCFs using {args.n_processes} processes...")
    with Pool(processes=args.n_processes) as pool:
        results = pool.map(process_single_vcf, args_list)
    
    print(f"Completed parsing {len(permuted_vcf_files)} VCFs to pickle files")
    step2_time = time.time()
    print(f"Step 2 completed in {(step2_time - step1_time)/60:.2f} minutes")
    
    # Step 3: Process Trees
    print("\n" + "="*80)
    print("STEP 3: BUILDING PHYLOGENETIC TREES AND ANALYZING TISSUE CLUSTERING")
    print("="*80)
    
    # Process observed data
    print("Processing observed data...")
    observed_percentage, observed_comparisons, observed_pairs = process_observed()
    print(f"Observed percentage: {observed_percentage:.2f}%, Comparisons: {observed_comparisons}, Pairs: {observed_pairs}")
    
    # Process permuted data
    permuted_pkl_files = glob.glob(os.path.join(args.permuted_pkls_dir, "permuted_*.pkl"))
    print(f"Found {len(permuted_pkl_files)} permuted pkl files to process.")
    
    print(f"Building trees and calculating tissue clustering using {args.n_processes} processes...")
    with Pool(processes=args.n_processes) as pool:
        permuted_percentages = pool.map(process_pkl, permuted_pkl_files)
    
    # Calculate p-value
    count_extreme = sum(1 for p in permuted_percentages if p >= observed_percentage)
    p_value = count_extreme / len(permuted_percentages)
    
    # Get group/tissue counts from samplesheet
    tissues = read_group_mapping(args.samplesheet)
    tissue_counts = {}
    for tissue in tissues.values():
        if tissue in tissue_counts:
            tissue_counts[tissue] += 1
        else:
            tissue_counts[tissue] = 1
    
    # Print results in a nicer format
    print("\n" + "="*80)
    print("PERMUTATION TEST RESULTS")
    print("="*80)
    print(f"Observed percentage of same-tissue pairs: {observed_percentage:.2f}%")
    print(f"P-value: {p_value:.4f} (based on {len(permuted_percentages)} permutations)")
    print(f"Mean permuted percentage: {np.mean(permuted_percentages):.2f}%")
    print(f"Median permuted percentage: {np.median(permuted_percentages):.2f}%")
    print(f"Max permuted percentage: {np.max(permuted_percentages):.2f}%")
    
    # Save statistics to a text file
    stats_file = f"{args.output_prefix}_permutation_stats.txt"
    with open(stats_file, "w") as f:
        f.write("Permutation Test (Shuffle Target MS) Statistics:\n")
        f.write("\nTissue distribution in the samplesheet:\n")
        for tissue, count in tissue_counts.items():
            f.write(f"  {tissue}: {count} samples\n")
        f.write(f"\nNumber of permutations: {args.n_permutations}\n")
        f.write(f"Observed Percentage of Same-Tissue Pairs: {observed_percentage:.2f}%\n")
        f.write(f"Unique Comparisons Made: {observed_comparisons}\n")
        f.write(f"Unique Same-Tissue Pairs: {observed_pairs}\n")
        f.write(f"P-value: {p_value:.4f}\n")
        f.write(f"Random Percentages - Mean: {np.mean(permuted_percentages):.2f}%, Median: {np.median(permuted_percentages):.2f}%, Max: {np.max(permuted_percentages):.2f}%\n")
    
    print(f"\nStatistics saved to '{stats_file}'")
    
    # Save permutation test data to JSON
    json_file = f"{args.output_prefix}_permutation_data.json"
    permutation_data = {
        "observed_percentage": observed_percentage,
        "random_percentages": permuted_percentages,
        "comparisons": observed_comparisons,
        "pairs": observed_pairs,
        "p_value": p_value
    }
    with open(json_file, "w") as f:
        json.dump(permutation_data, f, indent=4)
    print(f"Permutation test data saved to '{json_file}'")
    
    # Generate and save histogram
    plot_file = f"{args.output_prefix}_histogram.png"
    plt.figure(figsize=(10, 5), dpi=300)
    
    bins = np.arange(0, 100 + 5, 5)
    plt.hist(permuted_percentages, bins=bins, edgecolor='black')
    plt.axvline(x=observed_percentage, color='red', linestyle='dashed', linewidth=2)
    plt.title("Permutation Test: Shuffle Microsatellite Genotypes\n", fontsize=14)
    plt.xlabel("\nPercentage of Same-Tissue Pairs (%)")
    plt.ylabel("Frequency\n")
    plt.xlim(0, 100)
    plt.xticks(np.arange(0, 101, 10))
    
    # Remove right and top spines
    ax = plt.gca()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"Histogram saved to '{plot_file}'")
    plt.close()
    
    end_time = time.time()
    print(f"Step 3 completed in {(end_time - step2_time)/60:.2f} minutes")
    print(f"Total run time: {(end_time - start_time)/60:.2f} minutes")
    print("\n========== Permutation Test Complete ==========\n")

if __name__ == "__main__":
    # Parse command-line arguments
    args = parse_arguments()
    
    # Run the permutation test
    run_permutation_test() 