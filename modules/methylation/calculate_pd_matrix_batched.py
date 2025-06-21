#!/usr/bin/env python3
"""
Batched processing version that preloads reference data for maximum speed.
This approach minimizes I/O and prevents hangs with high process counts.
"""

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
import psutil
import resource
import gc

def get_optimal_process_count(n_processes=None):
    """Intelligently determine optimal process count."""
    if n_processes is not None and n_processes <= 8:
        return n_processes
    
    cpu_count = multiprocessing.cpu_count()
    memory_gb = psutil.virtual_memory().total / (1024**3)
    
    # More conservative for batched approach
    memory_limited = max(1, int(memory_gb * 0.6))  # Use 60% of RAM
    cpu_limited = max(1, min(cpu_count - 1, 12))   # Max 12 processes
    
    optimal = min(memory_limited, cpu_limited, 12)
    
    if n_processes is not None:
        optimal = min(optimal, n_processes)
    
    print(f"System: {cpu_count} CPUs, {memory_gb:.1f}GB RAM")
    print(f"Optimal processes: {optimal}")
    
    return optimal

def read_allc_fast(file_path, cpg_only=True, has_header=False):
    """Fast ALLC reader with memory optimization."""
    try:
        # Always read only the first 6 columns (the 7th column 'methylated' is never used)
        col_names = ['chromosome', 'position', 'strand', 'context', 'mc_count', 'total_count']
        col_types = {'chromosome': 'category', 'position': 'int32', 'strand': 'category', 
                    'context': 'category', 'mc_count': 'int16', 'total_count': 'int16'}
        
        # Read file efficiently - only first 6 columns
        df = pd.read_csv(file_path, sep='\t', header=None, names=col_names, 
                        comment='#', skiprows=1 if has_header else 0,
                        dtype=col_types, compression='gzip', engine='c',
                        usecols=range(6))  # Only read first 6 columns
        
        # Clean and filter
        df['chromosome'] = df['chromosome'].str.replace("chr", "", regex=False)
        
        if cpg_only:
            df = df[df['context'].str.startswith('CG', na=False)]
        
        # Calculate methylation fraction
        df['mc_frac'] = df['mc_count'] / np.maximum(df['total_count'], 1)
        
        # Create site key
        df['site_key'] = df['chromosome'].astype(str) + '_' + df['position'].astype(str)
        
        result = df[['site_key', 'mc_count', 'total_count', 'mc_frac']].copy()
        del df
        gc.collect()
        
        print(f"{os.path.basename(file_path)}: {len(result)} sites")
        return result
        
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return pd.DataFrame()

def calculate_dissimilarity_fast(sc_df, ref_df, min_reads=1, min_sites=300):
    """Fast vectorized dissimilarity calculation."""
    if sc_df.empty or ref_df.empty:
        return np.nan, 0
    
    # Fast merge
    merged = pd.merge(sc_df, ref_df, on='site_key', suffixes=('_sc', '_ref'), sort=False)
    
    if min_reads > 1:
        merged = merged[merged['total_count_sc'] >= min_reads]
    
    if len(merged) < min_sites:
        return np.nan, len(merged)
    
    # Vectorized calculation
    dissimilarities = np.abs(merged['mc_frac_sc'].values - merged['mc_frac_ref'].values) * 100
    return np.mean(dissimilarities), len(merged)

def preload_references(ref_files, cpg_only=True):
    """Preload all reference data to minimize I/O."""
    print("Preloading reference data...")
    ref_data = {}
    
    for ref_file in tqdm(ref_files, desc="Loading references"):
        ref_name = os.path.basename(ref_file).split('.')[0]
        ref_df = read_allc_fast(ref_file, cpg_only)
        if not ref_df.empty:
            ref_data[ref_name] = ref_df
    
    print(f"Preloaded {len(ref_data)} references")
    return ref_data

def process_single_cell_batch(sc_file, sc_name, ref_data_dict, min_reads=1, min_sites=300, cpg_only=True):
    """Process one single cell against all preloaded references."""
    try:
        # Load single cell once
        sc_data = read_allc_fast(sc_file, cpg_only)
        if sc_data.empty:
            return [(sc_name, ref_name, np.nan, 0) for ref_name in ref_data_dict.keys()]
        
        results = []
        for ref_name, ref_data in ref_data_dict.items():
            pd_score, shared_sites = calculate_dissimilarity_fast(
                sc_data, ref_data, min_reads, min_sites)
            results.append((sc_name, ref_name, pd_score, shared_sites))
            print(f"  {sc_name} vs {ref_name}: PD={pd_score:.2f}, sites={shared_sites}")
        
        del sc_data
        gc.collect()
        return results
        
    except Exception as e:
        print(f"Error processing {sc_name}: {e}")
        return [(sc_name, ref_name, np.nan, 0) for ref_name in ref_data_dict.keys()]

def process_files_batched(sc_files, ref_files, min_reads=1, min_sites=300, n_processes=None, cpg_only=True):
    """Batched processing with preloaded references."""
    n_processes = get_optimal_process_count(n_processes)
    
    # Preload all references in main process
    ref_data_dict = preload_references(ref_files, cpg_only)
    if not ref_data_dict:
        print("Error: No reference data loaded")
        return None, None
    
    # Extract names
    sc_names = [os.path.basename(f).replace('allc_', '').replace('.tsv.gz', '') for f in sc_files]
    ref_names = list(ref_data_dict.keys())
    
    # Initialize results
    pd_matrix = pd.DataFrame(index=sc_names, columns=ref_names, dtype=float)
    sites_matrix = pd.DataFrame(index=sc_names, columns=ref_names, dtype=int)
    
    print(f"Processing {len(sc_files)} single cells against {len(ref_names)} references...")
    print(f"Using {n_processes} processes")
    
    # Create tasks: one per single cell
    tasks = [(sc_file, sc_name, ref_data_dict, min_reads, min_sites, cpg_only) 
             for sc_file, sc_name in zip(sc_files, sc_names)]
    
    # Process with conservative settings
    with multiprocessing.Pool(processes=n_processes, maxtasksperchild=1) as pool:
        all_results = list(tqdm(
            pool.starmap(process_single_cell_batch, tasks),
            total=len(tasks),
            desc="Processing single cells"
        ))
    
    # Consolidate results
    for cell_results in all_results:
        for sc_name, ref_name, pd_score, shared_sites in cell_results:
            pd_matrix.loc[sc_name, ref_name] = pd_score
            sites_matrix.loc[sc_name, ref_name] = shared_sites
    
    return pd_matrix, sites_matrix

def expand_file_patterns(file_patterns):
    """Expand file patterns."""
    expanded_files = []
    for pattern in file_patterns:
        if '*' in pattern or '?' in pattern:
            expanded_files.extend(sorted(glob.glob(pattern)))
        else:
            expanded_files.append(pattern)
    return sorted(expanded_files)

def calculate_pairwise_dissimilarity_matrix_batched(sc_files, ref_files, output_dir='.', 
                                                   min_reads=1, min_sites=300, 
                                                   n_processes=None, cpg_only=True):
    """Main function for batched processing."""
    os.makedirs(output_dir, exist_ok=True)
    
    sc_files = expand_file_patterns(sc_files)
    ref_files = expand_file_patterns(ref_files)
    
    print(f"Found {len(sc_files)} single-cell files")
    print(f"Found {len(ref_files)} reference files")
    
    if not sc_files or not ref_files:
        print("Error: No input files found")
        return None, None
    
    # Process files
    pd_matrix, sites_matrix = process_files_batched(
        sc_files, ref_files, min_reads, min_sites, n_processes, cpg_only)
    
    if pd_matrix is None:
        return None, None
    
    # Save results
    pd_matrix_file = os.path.join(output_dir, 'pairwise_dissimilarity_matrix.csv')
    sites_matrix_file = os.path.join(output_dir, 'shared_sites_matrix.csv')
    
    pd_matrix.to_csv(pd_matrix_file)
    sites_matrix.to_csv(sites_matrix_file)
    
    print(f"Results saved to {pd_matrix_file} and {sites_matrix_file}")
    return pd_matrix, sites_matrix

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Batched pairwise dissimilarity matrix calculation with memory optimization')
    parser.add_argument('--sc_files', type=str, nargs='+', required=True,
                       help='List of single-cell ALLC files or patterns (e.g., "*.allc.tsv.gz")')
    parser.add_argument('--ref_files', type=str, nargs='+', required=True,
                       help='List of reference ALLC files or patterns (e.g., "*.tsv.gz")')
    parser.add_argument('--output_dir', type=str, default='.',
                       help='Output directory for results')
    parser.add_argument('--min_reads', type=int, default=1,
                       help='Minimum number of reads for a methylation site (applies to both single cell and reference files)')
    parser.add_argument('--min_sites', type=int, default=300,
                       help='Minimum number of shared sites required')
    parser.add_argument('--n_processes', type=int, default=None,
                       help='Number of processes to use for parallel processing. If None, automatically determines optimal count based on system resources (capped at 12 processes)')
    parser.add_argument('--all_cytosines', action='store_true',
                       help='Use all methylation contexts (not just CpG sites). Default: CpG sites only')
    
    args = parser.parse_args()
    
    calculate_pairwise_dissimilarity_matrix_batched(
        args.sc_files, args.ref_files, args.output_dir,
        args.min_reads, args.min_sites, args.n_processes,
        not args.all_cytosines  # cpg_only is True when all_cytosines is False (default: CpG only)
    ) 