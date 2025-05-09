#!/usr/bin/env python3
import pandas as pd
import numpy as np
import os
import glob
import argparse
from tqdm import tqdm
import re
import gzip
import multiprocessing
from functools import partial

def read_single_cell_allc(file_path, cpg_only=True):
    """
    Reads a single-cell ALLC file (gzipped TSV, no header).
    Returns a dictionary with (chr, pos, strand) as keys and (mc_count, total_count, mc_frac) as values.
    
    Parameters:
    -----------
    file_path : str
        Path to the ALLC file
    cpg_only : bool
        If True, only consider CpG sites (methylated cytosines in CG context)
        If False, consider all methylation contexts
    """
    result = {}
    context_counts = {'CG': 0, 'CH': 0, 'Other': 0}
    
    try:
        # Define column names since the file has no header
        with gzip.open(file_path, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) < 6:
                    continue
                
                chrom = parts[0].replace("chr", "")
                pos = int(parts[1])
                strand = parts[2]
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
                result[(chrom, pos, strand)] = (mc_count, total_count, mc_frac)
        
        print(f"File: {os.path.basename(file_path)}")
        print(f"  Contexts found: CG={context_counts['CG']}, CH={context_counts['CH']}, Other={context_counts['Other']}")
        print(f"  Sites used: {len(result)}")
        
        return result
    except Exception as e:
        print(f"Error reading single-cell file {file_path}: {e}")
        return {}

def read_reference_allc(file_path, cpg_only=True):
    """
    Reads a reference ALLC file (gzipped TSV, with header).
    Returns a dictionary with (chr, pos, strand) as keys and (mc_count, total_count, mc_frac) as values.
    
    Parameters:
    -----------
    file_path : str
        Path to the reference ALLC file
    cpg_only : bool
        If True, only consider CpG sites (methylated cytosines in CG context)
        If False, consider all methylation contexts
    """
    result = {}
    context_counts = {'CG': 0, 'CH': 0, 'Other': 0}
    
    try:
        # Read the file and process it line by line for better control over filtering
        with gzip.open(file_path, 'rt') as f:
            # Skip header line
            header = f.readline()
            
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 7:  # Reference files have at least 7 columns
                    continue
                
                chrom = parts[0].replace("chr", "")
                pos = int(parts[1])
                strand = parts[2]
                mc_class = parts[3]
                mc_count = int(parts[4])
                total_count = int(parts[5])
                mc_frac = float(parts[6])
                
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
                
                result[(chrom, pos, strand)] = (mc_count, total_count, mc_frac)
        
        print(f"Reference: {os.path.basename(file_path)}")
        print(f"  Contexts found: CG={context_counts['CG']}, CH={context_counts['CH']}, Other={context_counts['Other']}")
        print(f"  Sites used: {len(result)}")
        
        return result
    except Exception as e:
        print(f"Error reading reference file {file_path}: {e}")
        return {}

def calculate_pairwise_dissimilarity(sc_data, ref_data, min_reads=1, min_sites=10):
    """
    Calculate pairwise dissimilarity between a single cell and a reference cell type.
    Uses the absolute difference between methylation fractions.
    Returns PD score, number of shared sites, and list of dissimilarity values.
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

def process_single_cell(sc_file, sc_name, ref_files, ref_names, min_reads=1, min_sites=10, cpg_only=True):
    """
    Process a single cell against all reference files.
    Returns a dictionary with results for this cell.
    """
    results = {
        'sc_name': sc_name,
        'pd_scores': {},
        'shared_sites': {},
        'detailed': {}
    }
    
    sc_data = read_single_cell_allc(sc_file, cpg_only)
    
    if not sc_data:
        print(f"Skipping empty single cell data from {sc_file}")
        return results
    
    # Process against each reference
    for j, ref_file in enumerate(ref_files):
        ref_name = ref_names[j]
        ref_data = read_reference_allc(ref_file, cpg_only)
        
        if not ref_data:
            print(f"Skipping empty reference data from {ref_file}")
            results['pd_scores'][ref_name] = np.nan
            results['shared_sites'][ref_name] = 0
            continue
        
        # Calculate pairwise dissimilarity
        pd_score, shared_sites, dissimilarities = calculate_pairwise_dissimilarity(
            sc_data, ref_data, min_reads, min_sites)
        
        # Store results
        results['pd_scores'][ref_name] = pd_score
        results['shared_sites'][ref_name] = shared_sites
        results['detailed'][ref_name] = dissimilarities
        
        # Print detailed results
        print(f"  {sc_name} vs {ref_name}: PD={pd_score:.2f}, Shared sites={shared_sites}")
    
    return results

def process_files(sc_files, ref_files, min_reads=1, min_sites=10, n_processes=None, cpg_only=True):
    """
    Process all single-cell files against all reference files using multiprocessing.
    Returns a DataFrame with dissimilarity scores and a dictionary with detailed results.
    """
    # Use up to 8 processes by default, or whatever the user specifies
    if n_processes is None:
        n_processes = min(8, max(1, multiprocessing.cpu_count() // 2))
    else:
        # Cap the maximum number of processes to avoid system overload
        n_processes = min(16, n_processes)
    
    # Extract cell and reference names
    sc_names = [os.path.basename(f).replace('allc_Methyl_', '').replace('.tsv.gz', '') for f in sc_files]
    ref_names = [os.path.basename(f).replace('_aggregated.tsv.gz', '') for f in ref_files]
    
    # Initialize results
    pd_matrix = pd.DataFrame(index=sc_names, columns=ref_names, dtype=float)
    sites_matrix = pd.DataFrame(index=sc_names, columns=ref_names, dtype=int)
    detailed_results = {}
    
    # Create a partial function with fixed parameters
    process_func = partial(process_single_cell, 
                          ref_files=ref_files, 
                          ref_names=ref_names,
                          min_reads=min_reads, 
                          min_sites=min_sites,
                          cpg_only=cpg_only)
    
    # Process each single cell in parallel
    print(f"Processing single cells using {n_processes} processes...")
    print(f"Using {'CpG sites only' if cpg_only else 'all methylation contexts'}")
    
    with multiprocessing.Pool(processes=n_processes) as pool:
        results = list(tqdm(
            pool.starmap(process_func, zip(sc_files, sc_names)), 
            total=len(sc_files),
            desc="Processing single cells"
        ))
    
    # Consolidate results
    for result in results:
        sc_name = result['sc_name']
        detailed_results[sc_name] = result['detailed']
        
        for ref_name in ref_names:
            pd_matrix.loc[sc_name, ref_name] = result['pd_scores'].get(ref_name, np.nan)
            sites_matrix.loc[sc_name, ref_name] = result['shared_sites'].get(ref_name, 0)
    
    return pd_matrix, sites_matrix, detailed_results

def calculate_pairwise_dissimilarity_matrix(sc_dir, ref_dir, output_dir='.', min_reads=1, min_sites=10, n_processes=None, cpg_only=True):
    """
    Calculate the pairwise dissimilarity matrix between single cells and reference cell types.
    
    Parameters:
    -----------
    sc_dir : str
        Directory containing single-cell files
    ref_dir : str
        Directory containing reference files
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
    
    # Prepare output file prefix based on methylation context
    context_label = "cpg" if cpg_only else "all"
    
    # Define output files
    pd_matrix_file = os.path.join(output_dir, f'{context_label}_pairwise_dissimilarity.csv')
    sites_matrix_file = os.path.join(output_dir, f'{context_label}_shared_sites.csv')
    
    # Find input files
    sc_files = sorted(glob.glob(os.path.join(sc_dir, "*tsv.gz")))
    ref_files = sorted(glob.glob(os.path.join(ref_dir, "*aggregated.tsv.gz")))
    
    print(f"Found {len(sc_files)} single-cell files")
    print(f"Found {len(ref_files)} reference files")
    
    if not sc_files or not ref_files:
        print("Error: No input files found")
        return None, None
    
    # Process files to calculate pairwise dissimilarity
    print("Calculating pairwise dissimilarity...")
    pd_matrix, sites_matrix, detailed_results = process_files(
        sc_files, ref_files, min_reads, min_sites, n_processes, cpg_only)
    
    # Save the raw PD matrix
    pd_matrix.to_csv(pd_matrix_file)
    sites_matrix.to_csv(sites_matrix_file)
    
    print(f"Pairwise dissimilarity matrix saved to {pd_matrix_file}")
    print(f"Shared sites matrix saved to {sites_matrix_file}")
    
    return pd_matrix, sites_matrix

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate pairwise dissimilarity matrix based on methylation data')
    parser.add_argument('--sc_dir', type=str, required=True, help='Directory containing single-cell ALLC files')
    parser.add_argument('--ref_dir', type=str, required=True, help='Directory containing reference ALLC files')
    parser.add_argument('--output_dir', type=str, default='.', help='Output directory for results')
    parser.add_argument('--min_reads', type=int, default=1, help='Minimum number of reads for a methylation site')
    parser.add_argument('--min_sites', type=int, default=10, help='Minimum number of shared sites required')
    parser.add_argument('--n_processes', type=int, default=None, 
                       help='Number of processes to use for parallel processing. Default: half of available CPU cores')
    parser.add_argument('--all_cytosines', action='store_true', 
                       help='Use all methylation contexts (not just CpG sites)')
    
    args = parser.parse_args()
    
    calculate_pairwise_dissimilarity_matrix(
        args.sc_dir, 
        args.ref_dir, 
        args.output_dir,
        args.min_reads,
        args.min_sites,
        args.n_processes,
        not args.all_cytosines  # cpg_only is True when all_cytosines is False
    ) 