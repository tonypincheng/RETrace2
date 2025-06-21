#!/usr/bin/env python3
import pandas as pd
import numpy as np
import os
import glob
import argparse
from tqdm import tqdm
import gzip
import multiprocessing
from functools import partial
import json

def read_allc(file_path, cpg_only=True, has_header=False):
    """
    Reads an ALLC file (gzipped TSV, no header by default).
    Returns a dictionary with (chr, pos) as keys and (mc_count, total_count, mc_frac) as values.
    
    Parameters:
    -----------
    file_path : str
        Path to the ALLC file
    cpg_only : bool
        If True, only consider CpG sites (methylated cytosines in CG context)
        If False, consider all methylation contexts
    has_header : bool
        If True, skips the first line (neither single cell nor reference files currently have headers)
    """
    result = {}
    context_counts = {'CG': 0, 'CH': 0, 'Other': 0}
    
    try:
        with gzip.open(file_path, 'rt') as f:
            # Skip header if specified
            if has_header:
                header = f.readline()
                
            for line in f:
                if line.startswith('#'):
                    continue
                    
                parts = line.strip().split('\t')
                if len(parts) < 6:  # Need at least 6 columns
                    continue
                
                chrom = parts[0].replace("chr", "")
                pos = int(parts[1])
                strand = parts[2]  # Reading strand but not using in key
                mc_class = parts[3]
                mc_count = int(parts[4])
                total_count = int(parts[5])
                
                # Track context counts
                if mc_class.startswith('CG'):
                    context_counts['CG'] += 1
                elif mc_class.startswith('C'):
                    context_counts['CH'] += 1
                else:
                    context_counts['Other'] += 1
                
                # Filter based on cpg_only parameter
                if cpg_only and not mc_class.startswith('CG'):
                    continue
                
                # Calculate methylation fraction
                mc_frac = mc_count / total_count if total_count > 0 else 0
                result[(chrom, pos)] = (mc_count, total_count, mc_frac)
        
        print(f"{os.path.basename(file_path)}")
        print(f"  Contexts found: CG={context_counts['CG']}, CH={context_counts['CH']}, Other={context_counts['Other']}")
        print(f"  Sites used: {len(result)}")
        
        return result
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return {}

def calculate_pairwise_dissimilarity(sc_data, ref_data, min_reads=1, min_sites=300):
    """
    Calculate pairwise dissimilarity between a single cell and a reference cell type.
    Uses the absolute difference between methylation fractions.
    Sites are identified by chromosome and position only (strand-agnostic).
    Returns PD score, number of shared sites, and list of dissimilarity values.
    
    Parameters:
    -----------
    sc_data : dict
        Dictionary with (chr, pos) as keys and (mc_count, total_count, mc_frac) as values
    ref_data : dict
        Dictionary with (chr, pos) as keys and (mc_count, total_count, mc_frac) as values
    min_reads : int
        Minimum read coverage for a methylation site
    min_sites : int
        Minimum number of shared sites required
        
    Returns:
    --------
    avg_dissimilarity : float
        Average pairwise dissimilarity score
    len(filtered_sites) : int
        Number of shared sites used for calculation
    dissimilarities : list
        List of dissimilarity values for each shared site
    """
    # Find overlapping sites
    overlapping_sites = set(sc_data.keys()).intersection(set(ref_data.keys()))
    
    # Filter sites by minimum read coverage
    filtered_sites = [site for site in overlapping_sites 
                     if sc_data[site][1] >= min_reads]
    
    # If too few overlapping sites, return NaN
    if len(filtered_sites) < min_sites:
        return np.nan, len(filtered_sites), []
    
    # Calculate dissimilarity for each site
    dissimilarities = []
    for site in filtered_sites:
        sc_meth = sc_data[site][2]  # methylation fraction
        ref_meth = ref_data[site][2]  # methylation fraction
        dissimilarity = abs(sc_meth - ref_meth) * 100  # Dissimilarity as percentage
        dissimilarities.append(dissimilarity)
    
    # Calculate average dissimilarity
    avg_dissimilarity = np.mean(dissimilarities)
    
    return avg_dissimilarity, len(filtered_sites), dissimilarities

def process_single_comparison(sc_file, sc_name, ref_file, ref_name, min_reads=1, min_sites=10, cpg_only=True):
    """
    Process a single cell against a single reference file.
    This creates many small tasks that can run in parallel without memory issues.
    """
    # Load single cell data
    sc_data = read_allc(sc_file, cpg_only)
    if not sc_data:
        return {
            'sc_name': sc_name,
            'ref_name': ref_name,
            'pd_score': np.nan,
            'shared_sites': 0
        }
    
    # Load single reference
    ref_data = read_allc(ref_file, cpg_only)
    if not ref_data:
        return {
            'sc_name': sc_name,
            'ref_name': ref_name,
            'pd_score': np.nan,
            'shared_sites': 0
        }
    
    # Calculate pairwise dissimilarity
    pd_score, shared_sites, dissimilarities = calculate_pairwise_dissimilarity(
        sc_data, ref_data, min_reads, min_sites)
    
    # Print progress
    print(f"  {sc_name} vs {ref_name}: PD={pd_score:.2f}, Shared sites={shared_sites}")
    
    return {
        'sc_name': sc_name,
        'ref_name': ref_name,
        'pd_score': pd_score,
        'shared_sites': shared_sites
    }

def process_files(sc_files, ref_files, min_reads=1, min_sites=300, n_processes=1, cpg_only=True):
    """
    Process all single-cell files against all reference files using task-per-comparison.
    Returns a DataFrame with dissimilarity scores and a dictionary with detailed results.
    """
    # Use the requested number of processes
    n_processes = min(32, n_processes)
    
    # Extract names
    sc_names = [os.path.basename(f).replace('allc_', '').replace('.tsv.gz', '') for f in sc_files]
    ref_names = [os.path.basename(f).split('.')[0] for f in ref_files]
    
    # Initialize results
    pd_matrix = pd.DataFrame(index=sc_names, columns=ref_names, dtype=float)
    sites_matrix = pd.DataFrame(index=sc_names, columns=ref_names, dtype=int)
    
    # Create tasks: one task per (single_cell, reference) comparison
    tasks = []
    for sc_file, sc_name in zip(sc_files, sc_names):
        for ref_file, ref_name in zip(ref_files, ref_names):
            tasks.append((sc_file, sc_name, ref_file, ref_name))
    
    print(f"Processing {len(sc_files)} single cells against {len(ref_files)} references...")
    print(f"Total comparisons: {len(tasks)} (using {n_processes} processes)")
    print(f"Using {'CpG sites only' if cpg_only else 'all methylation contexts'}")
    
    # Create partial function
    process_func = partial(process_single_comparison, 
                          min_reads=min_reads, 
                          min_sites=min_sites,
                          cpg_only=cpg_only)
    
    # Process all comparisons in parallel
    with multiprocessing.Pool(processes=n_processes, maxtasksperchild=10) as pool:
        results = list(tqdm(
            pool.starmap(process_func, tasks),
            total=len(tasks),
            desc="Processing comparisons"
        ))
    
    # Consolidate results
    for result in results:
        sc_name = result['sc_name']
        ref_name = result['ref_name']
        
        # Store results
        pd_matrix.loc[sc_name, ref_name] = result['pd_score']
        sites_matrix.loc[sc_name, ref_name] = result['shared_sites']
    
    return pd_matrix, sites_matrix

def expand_file_patterns(file_patterns):
    """
    Expand file patterns into a list of matching files.
    Handles both direct file paths and wildcard patterns.
    """
    expanded_files = []
    for pattern in file_patterns:
        if '*' in pattern or '?' in pattern:
            # Handle wildcard pattern
            expanded_files.extend(sorted(glob.glob(pattern)))
        else:
            # Handle direct file path
            expanded_files.append(pattern)
    return sorted(expanded_files)

def calculate_pairwise_dissimilarity_matrix(sc_files, ref_files, output_dir='.', min_reads=1, min_sites=300, n_processes=None, cpg_only=True):
    """
    Calculate the pairwise dissimilarity matrix between single cells and reference cell types.
    
    Parameters:
    -----------
    sc_files : list
        List of paths to single-cell ALLC files or file patterns (e.g., "*.allc.tsv.gz")
    ref_files : list
        List of paths to reference ALLC files or file patterns (e.g., "*.tsv.gz")
    output_dir : str
        Directory to save results
    min_reads : int
        Minimum read coverage for a methylation site
    min_sites : int
        Minimum number of shared sites required
    n_processes : int, optional
        Number of processes to use for parallel processing. 
        If None, uses half of available CPU cores.
    cpg_only : bool
        If True, only consider CpG sites (methylated cytosines in CG context)
        If False, consider all methylation contexts
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Expand file patterns
    sc_files = expand_file_patterns(sc_files)
    ref_files = expand_file_patterns(ref_files)
    
    # Define output files
    pd_matrix_file = os.path.join(output_dir, 'pairwise_dissimilarity_matrix.csv')
    sites_matrix_file = os.path.join(output_dir, 'shared_sites_matrix.csv')
    
    print(f"Found {len(sc_files)} single-cell files")
    print(f"Found {len(ref_files)} reference files")
    
    if not sc_files or not ref_files:
        print("Error: No input files found")
        return None, None
    
    # Process files to calculate pairwise dissimilarity
    print("Calculating pairwise dissimilarity...")
    pd_matrix, sites_matrix = process_files(
        sc_files, ref_files, min_reads, min_sites, n_processes, cpg_only)
    
    # Save the raw PD matrix
    pd_matrix.to_csv(pd_matrix_file)
    sites_matrix.to_csv(sites_matrix_file)
    
    print(f"Pairwise dissimilarity matrix saved to {pd_matrix_file}")
    print(f"Shared sites matrix saved to {sites_matrix_file}")
    
    return pd_matrix, sites_matrix

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate pairwise dissimilarity matrix based on methylation data')
    parser.add_argument('--sc_files', type=str, nargs='+', required=True, 
                       help='List of single-cell ALLC files or patterns (e.g., "*.allc.tsv.gz")')
    parser.add_argument('--ref_files', type=str, nargs='+', required=True, 
                       help='List of reference ALLC files or patterns (e.g., "*.tsv.gz")')
    parser.add_argument('--output_dir', type=str, default='.', 
                       help='Output directory for results')
    parser.add_argument('--min_reads', type=int, default=1, 
                       help='Minimum number of reads for a methylation site')
    parser.add_argument('--min_sites', type=int, default=300, 
                       help='Minimum number of shared sites required')
    parser.add_argument('--n_processes', type=int, default=1, 
                       help='Number of processes to use for parallel processing. Default: 1, Maximum: 32')
    parser.add_argument('--all_cytosines', action='store_true', 
                       help='Use all methylation contexts (not just CpG sites)')
    
    args = parser.parse_args()
    
    calculate_pairwise_dissimilarity_matrix(
        args.sc_files, 
        args.ref_files, 
        args.output_dir,
        args.min_reads,
        args.min_sites,
        args.n_processes,
        not args.all_cytosines  # cpg_only is True when all_cytosines is False
    ) 