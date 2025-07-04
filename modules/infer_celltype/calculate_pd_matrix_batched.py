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
import multiprocessing
import psutil
import gc
import time
from datetime import datetime

# Global variable to track report file path
REPORT_FILE = None

def write_to_report(message):
    """Write a message to the report file with timestamp."""
    global REPORT_FILE
    if REPORT_FILE:
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        with open(REPORT_FILE, 'a') as f:
            f.write(f"[{timestamp}] {message}\n")

def initialize_report_file(output_dir):
    """Initialize the report file for data validation issues."""
    global REPORT_FILE
    REPORT_FILE = os.path.join(output_dir, 'data_validation_report.txt')
    # Clear existing report file
    with open(REPORT_FILE, 'w') as f:
        f.write("RETrace2 Data Validation Report\n")
        f.write("=" * 50 + "\n\n")
    print(f"Data validation report will be written to: {REPORT_FILE}")

def get_memory_usage_mb():
    """Get current memory usage in MB."""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / 1024 / 1024

def estimate_reference_memory(ref_files, sample_size=3):
    """Estimate memory usage of reference data by sampling."""
    if not ref_files:
        return 0
    
    print("Estimating reference data memory usage...")
    sample_files = ref_files[:min(sample_size, len(ref_files))]
    
    total_estimated_mb = 0
    for ref_file in sample_files:
        initial_mem = get_memory_usage_mb()
        try:
            ref_df = read_allc_fast(ref_file, cpg_only=True)
            if not ref_df.empty:
                current_mem = get_memory_usage_mb()
                file_mem = current_mem - initial_mem
                total_estimated_mb += file_mem
                del ref_df
                gc.collect()
        except Exception as e:
            print(f"Warning: Could not estimate memory for {ref_file}: {e}")
            # Use a conservative estimate
            total_estimated_mb += 50  # 50MB per file as fallback
    
    # Extrapolate to all files
    avg_mem_per_file = total_estimated_mb / len(sample_files)
    total_estimated_mb = avg_mem_per_file * len(ref_files)
    
    print(f"Estimated reference data memory: {total_estimated_mb:.1f} MB")
    return total_estimated_mb

def get_optimal_process_count(n_processes=None, ref_memory_mb=0):
    """Intelligently determine optimal process count with memory considerations."""
    if n_processes is not None and n_processes <= 4:
        return n_processes
    
    cpu_count = multiprocessing.cpu_count()
    memory_gb = psutil.virtual_memory().total / (1024**3)
    available_memory_gb = psutil.virtual_memory().available / (1024**3)
    
    print(f"System: {cpu_count} CPUs, {memory_gb:.1f}GB RAM, {available_memory_gb:.1f}GB available")
    
    # Conservative memory calculation
    # Each process will hold a copy of reference data plus single cell data
    ref_memory_gb = ref_memory_mb / 1024
    memory_per_process_gb = ref_memory_gb + 0.5  # 0.5GB buffer per process
    
    # Use only 50% of available memory to be very conservative
    usable_memory_gb = available_memory_gb * 0.5
    memory_limited = max(1, int(usable_memory_gb / memory_per_process_gb))
    
    # Conservative CPU limit - especially important for large datasets
    if ref_memory_gb > 10:  # If references are large
        cpu_limited = max(1, min(cpu_count // 2, 12))  # Conservative
    else:
        cpu_limited = max(1, min(cpu_count - 1, 20))   # Set max cpu to 20 to avoid memory issues
    
    optimal = min(memory_limited, cpu_limited)
    
    if n_processes is not None:
        optimal = min(optimal, n_processes)
    
    print(f"Memory per process: {memory_per_process_gb:.1f}GB")
    print(f"Memory-limited processes: {memory_limited}")
    print(f"CPU-limited processes: {cpu_limited}")
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
                        dtype=col_types, compression='infer', engine='c',
                        usecols=range(6))  # Only read first 6 columns
        
        # Clean and filter
        df['chromosome'] = df['chromosome'].str.replace("chr", "", regex=False)
        
        if cpg_only:
            df = df[df['context'].str.startswith('CG', na=False)]
        
        # Validate data integrity before calculating methylation fraction
        initial_sites = len(df)
        
        # Check for impossible cases: mc_count > total_count
        invalid_counts = df['mc_count'] > df['total_count']
        if invalid_counts.any():
            n_invalid = invalid_counts.sum()
            write_to_report(f"File {file_path}: Found {n_invalid} sites where mc_count > total_count. These sites will be excluded.")
            df = df[~invalid_counts]
        
        # Remove sites with zero total_count (no reads)
        zero_coverage = df['total_count'] == 0
        if zero_coverage.any():
            n_zero = zero_coverage.sum()
            write_to_report(f"File {file_path}: Found {n_zero} sites with zero total_count. These sites will be excluded.")
            df = df[~zero_coverage]
        
        if df.empty:
            write_to_report(f"File {file_path}: No valid sites remaining after data validation.")
            return pd.DataFrame()
        
        # Calculate methylation fraction (now safe - no division by zero or invalid ratios)
        df['mc_frac'] = df['mc_count'] / df['total_count']
        
        # Final validation: ensure mc_frac is between 0 and 1
        invalid_fractions = (df['mc_frac'] < 0) | (df['mc_frac'] > 1)
        if invalid_fractions.any():
            n_invalid_frac = invalid_fractions.sum()
            write_to_report(f"File {file_path}: Found {n_invalid_frac} sites with invalid mc_frac (not between 0-1). These sites will be excluded.")
            df = df[~invalid_fractions]
        
        if df.empty:
            write_to_report(f"File {file_path}: No valid sites remaining after methylation fraction validation.")
            return pd.DataFrame()
        
        # Create site key and consider adding strand and GC context in future
        df['site_key'] = df['chromosome'].astype(str) + '_' + df['position'].astype(str)
        
        result = df[['site_key', 'mc_count', 'total_count', 'mc_frac']].copy()
        del df
        gc.collect()
        
        filtered_sites = len(result)
        if initial_sites != filtered_sites:
            write_to_report(f"File {file_path}: Filtered from {initial_sites} to {filtered_sites} sites ({initial_sites - filtered_sites} sites excluded).")
        
        print(f"{os.path.basename(file_path)}: {filtered_sites} valid sites (filtered from {initial_sites})")
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
    
    # Validate dissimilarity values should be between 0 and 100
    invalid_dissim = (dissimilarities < 0) | (dissimilarities > 100)
    if invalid_dissim.any():
        n_invalid = invalid_dissim.sum()
        max_invalid = np.max(dissimilarities[invalid_dissim]) if n_invalid > 0 else 0
        write_to_report(f"Found {n_invalid} invalid dissimilarity values (not between 0-100). Max invalid value: {max_invalid:.2f}. These will be excluded from calculation.")
        # Keep only valid dissimilarities
        dissimilarities = dissimilarities[~invalid_dissim]
        
        # If no valid dissimilarities remain, return NaN
        if len(dissimilarities) == 0:
            write_to_report("No valid dissimilarity values remaining after filtering.")
            return np.nan, len(merged)
    
    return np.mean(dissimilarities), len(merged)

def preload_references(ref_files, cpg_only=True):
    """Preload all reference data to minimize I/O with memory monitoring."""
    print("Preloading reference data...")
    ref_data = {}
    initial_memory = get_memory_usage_mb()
    
    for i, ref_file in enumerate(tqdm(ref_files, desc="Loading references")):
        ref_name = os.path.basename(ref_file).split('.')[0]
        ref_df = read_allc_fast(ref_file, cpg_only)
        if not ref_df.empty:
            ref_data[ref_name] = ref_df
        
        # Monitor memory every 10 files
        if (i + 1) % 10 == 0:
            current_memory = get_memory_usage_mb()
            memory_used = current_memory - initial_memory
            print(f"  Loaded {i+1}/{len(ref_files)} refs, using {memory_used:.1f} MB")
            
            # Check if we're approaching memory limits
            available_memory_gb = psutil.virtual_memory().available / (1024**3)
            if available_memory_gb < 2:  # Less than 2GB available
                print(f"Warning: Low memory detected ({available_memory_gb:.1f}GB available)")
    
    final_memory = get_memory_usage_mb()
    total_memory_used = final_memory - initial_memory
    print(f"Preloaded {len(ref_data)} references using {total_memory_used:.1f} MB")
    
    return ref_data

def process_single_cell_batch(sc_file, sc_name, ref_data_dict, min_reads=1, min_sites=300, cpg_only=True):
    """Process one single cell against all preloaded references with memory monitoring."""
    try:
        # Monitor memory at start
        initial_mem = get_memory_usage_mb()
        
        # Load single cell once
        sc_data = read_allc_fast(sc_file, cpg_only)
        if sc_data.empty:
            return [(sc_name, ref_name, np.nan, 0) for ref_name in ref_data_dict.keys()]
        
        results = []
        for ref_name, ref_data in ref_data_dict.items():
            pd_score, shared_sites = calculate_dissimilarity_fast(
                sc_data, ref_data, min_reads, min_sites)
            results.append((sc_name, ref_name, pd_score, shared_sites))
        
        # Clean up and report memory usage
        del sc_data
        gc.collect()
        
        final_mem = get_memory_usage_mb()
        if final_mem - initial_mem > 100:  # If process used more than 100MB
            print(f"  {sc_name}: memory used {final_mem - initial_mem:.1f} MB")
        
        return results
        
    except Exception as e:
        print(f"Error processing {sc_name}: {e}")
        return [(sc_name, ref_name, np.nan, 0) for ref_name in ref_data_dict.keys()]

def process_files_batched(sc_files, ref_files, min_reads=1, min_sites=300, n_processes=None, cpg_only=True):
    """Batched processing with preloaded references and memory optimization."""
    # Estimate memory requirements first
    ref_memory_mb = estimate_reference_memory(ref_files)
    
    # Get optimal process count based on memory constraints
    n_processes = get_optimal_process_count(n_processes, ref_memory_mb)
    
    # Check if we have enough memory for this approach
    available_memory_gb = psutil.virtual_memory().available / (1024**3)
    required_memory_gb = (ref_memory_mb / 1024) * n_processes * 1.5  # 1.5x safety factor
    
    if required_memory_gb > available_memory_gb:
        print(f"WARNING: Estimated memory requirement ({required_memory_gb:.1f}GB) exceeds available memory ({available_memory_gb:.1f}GB)")
        print("Consider reducing the number of processes or using a machine with more memory")
        
        # Force very conservative settings
        n_processes = max(1, min(n_processes, int(available_memory_gb / (ref_memory_mb / 1024))))
        print(f"Reducing processes to {n_processes} to fit in available memory")
    
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
    
    # Process with conservative settings and memory monitoring
    try:
        with multiprocessing.Pool(processes=n_processes, maxtasksperchild=1) as pool:
            all_results = list(tqdm(
                pool.starmap(process_single_cell_batch, tasks),
                total=len(tasks),
                desc="Processing single cells"
            ))
    except Exception as e:
        print(f"Error during multiprocessing: {e}")
        print("Falling back to sequential processing...")
        
        # Fallback to sequential processing
        all_results = []
        for task in tqdm(tasks, desc="Processing single cells (sequential)"):
            result = process_single_cell_batch(*task)
            all_results.append(result)
    
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

def process_single_cell_streaming(sc_file, sc_name, ref_files, min_reads=1, min_sites=300, cpg_only=True):
    """Process one single cell against references without preloading all references."""
    try:
        # Load single cell once
        sc_data = read_allc_fast(sc_file, cpg_only)
        if sc_data.empty:
            return [(sc_name, os.path.basename(ref_file).split('.')[0], np.nan, 0) for ref_file in ref_files]
        
        results = []
        for ref_file in ref_files:
            ref_name = os.path.basename(ref_file).split('.')[0]
            ref_data = read_allc_fast(ref_file, cpg_only)
            
            if ref_data.empty:
                results.append((sc_name, ref_name, np.nan, 0))
            else:
                pd_score, shared_sites = calculate_dissimilarity_fast(
                    sc_data, ref_data, min_reads, min_sites)
                results.append((sc_name, ref_name, pd_score, shared_sites))
            
            # Clean up reference data immediately
            del ref_data
            gc.collect()
        
        del sc_data
        gc.collect()
        return results
        
    except Exception as e:
        print(f"Error processing {sc_name}: {e}")
        return [(sc_name, os.path.basename(ref_file).split('.')[0], np.nan, 0) for ref_file in ref_files]

def process_files_streaming(sc_files, ref_files, min_reads=1, min_sites=300, n_processes=None, cpg_only=True):
    """Memory-efficient streaming approach that doesn't preload references."""
    n_processes = get_optimal_process_count(n_processes, 0)  # No reference memory preloaded
    
    # Extract names
    sc_names = [os.path.basename(f).replace('allc_', '').replace('.tsv.gz', '') for f in sc_files]
    ref_names = [os.path.basename(f).split('.')[0] for f in ref_files]
    
    # Initialize results
    pd_matrix = pd.DataFrame(index=sc_names, columns=ref_names, dtype=float)
    sites_matrix = pd.DataFrame(index=sc_names, columns=ref_names, dtype=int)
    
    print(f"Processing {len(sc_files)} single cells against {len(ref_names)} references (streaming mode)...")
    print(f"Using {n_processes} processes")
    
    # Create tasks: one per single cell
    tasks = [(sc_file, sc_name, ref_files, min_reads, min_sites, cpg_only) 
             for sc_file, sc_name in zip(sc_files, sc_names)]
    
    # Process with multiprocessing
    with multiprocessing.Pool(processes=n_processes, maxtasksperchild=1) as pool:
        all_results = list(tqdm(
            pool.starmap(process_single_cell_streaming, tasks),
            total=len(tasks),
            desc="Processing single cells (streaming)"
        ))
    
    # Consolidate results
    for cell_results in all_results:
        for sc_name, ref_name, pd_score, shared_sites in cell_results:
            pd_matrix.loc[sc_name, ref_name] = pd_score
            sites_matrix.loc[sc_name, ref_name] = shared_sites
    
    return pd_matrix, sites_matrix

def calculate_pairwise_dissimilarity_matrix_batched(sc_files, ref_files, output_dir='.', 
                                                   min_reads=1, min_sites=300, 
                                                   n_processes=None, cpg_only=True,
                                                   force_streaming=False, max_ref_memory_gb=8):
    """Main function for batched processing with automatic memory management."""
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize report file for data validation issues
    initialize_report_file(output_dir)
    
    sc_files = expand_file_patterns(sc_files)
    ref_files = expand_file_patterns(ref_files)
    
    print("="*60)
    print("RETrace2 Pairwise Dissimilarity Matrix Calculator")
    print("="*60)
    print(f"Found {len(sc_files)} single-cell files")
    print(f"Found {len(ref_files)} reference files")
    
    if not sc_files or not ref_files:
        print("Error: No input files found")
        return None, None
    
    print("\nMEMORY OPTIMIZATION:")
    print("- BATCHED MODE: Preloads all references for maximum speed")
    print("- STREAMING MODE: Loads references on-demand to save memory")
    
    # Decide processing approach based on memory constraints
    if force_streaming:
        print(f"\n→ Using STREAMING MODE (forced via --force_streaming)")
        pd_matrix, sites_matrix = process_files_streaming(
            sc_files, ref_files, min_reads, min_sites, n_processes, cpg_only)
    else:
        # Estimate memory and decide approach
        ref_memory_mb = estimate_reference_memory(ref_files)
        ref_memory_gb = ref_memory_mb / 1024
        
        if ref_memory_gb > max_ref_memory_gb:
            print(f"\n→ Reference data ({ref_memory_gb:.1f}GB) exceeds limit ({max_ref_memory_gb}GB)")
            print("→ Automatically switching to STREAMING MODE to avoid memory issues")
            pd_matrix, sites_matrix = process_files_streaming(
                sc_files, ref_files, min_reads, min_sites, n_processes, cpg_only)
        else:
            print(f"\n→ Using BATCHED MODE with preloaded references ({ref_memory_gb:.1f}GB)")
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
    print(f"Data validation report saved to {REPORT_FILE}")
    return pd_matrix, sites_matrix

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Batched pairwise dissimilarity matrix calculation with automatic memory optimization.\n\n'
                   'This script supports two processing modes:\n'
                   '1. BATCHED MODE (default): Preloads all reference data for maximum speed\n'
                   '2. STREAMING MODE: Loads references on-demand to minimize memory usage\n\n'
                   'The script automatically switches to streaming mode when:\n'
                   '- Reference data exceeds --max_ref_memory_gb limit\n'
                   '- Memory estimation indicates insufficient RAM\n'
                   '- --force_streaming flag is used\n\n'
                   'Streaming mode is slower but handles very large reference datasets that cannot fit in memory.',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--sc_files', type=str, nargs='+', required=True,
                       help='List of single-cell ALLC files or patterns (e.g., "*.allc.tsv.gz")')
    parser.add_argument('--ref_files', type=str, nargs='+', required=True,
                       help='List of reference ALLC files or patterns (e.g., "*.tsv.gz")')
    parser.add_argument('--output_dir', type=str, default='.',
                       help='Output directory for results')
    parser.add_argument('--min_reads', type=int, default=1,
                       help='Minimum number of reads for a methylation site (applies to both single cell and reference files)')
    parser.add_argument('--min_sites', type=int, default=100,
                       help='Minimum number of shared sites required')
    parser.add_argument('--n_processes', type=int, default=None,
                       help='Number of processes to use for parallel processing. If None, automatically determines optimal count based on system resources')
    parser.add_argument('--all_cytosines', action='store_true',
                       help='Use all methylation contexts (not just CpG sites). Default: CpG sites only')
    parser.add_argument('--force_streaming', action='store_true',
                       help='Force streaming mode (don\'t preload references). Use this for very large reference datasets.')
    parser.add_argument('--max_ref_memory_gb', type=float, default=8.0,
                       help='Maximum memory (GB) to use for preloading reference data. If references exceed this, streaming mode will be used.')
    
    args = parser.parse_args()
    
    # Record start time
    start_time = time.time()
    print(f"Starting analysis at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    calculate_pairwise_dissimilarity_matrix_batched(
        args.sc_files, args.ref_files, args.output_dir,
        args.min_reads, args.min_sites, args.n_processes,
        not args.all_cytosines,  # cpg_only is True when all_cytosines is False (default: CpG only)
        args.force_streaming, args.max_ref_memory_gb
    )
    
    # Record end time and calculate elapsed time
    end_time = time.time()
    elapsed_time = end_time - start_time
    
    # Format and print elapsed time
    hours = int(elapsed_time // 3600)
    minutes = int((elapsed_time % 3600) // 60)
    seconds = elapsed_time % 60
    
    print(f"\nAnalysis completed at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Total execution time: {hours:02d}:{minutes:02d}:{seconds:05.2f}")
    print(f"Total execution time: {elapsed_time:.2f} seconds") 