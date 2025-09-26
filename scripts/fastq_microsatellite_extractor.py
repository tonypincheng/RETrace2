#!/usr/bin/env python3
"""
FASTQ Microsatellite Extractor for RETrace2 Paper

This script processes FASTQ sequencing files to extract UMI (Unique Molecular Identifier) 
sequences, library indices, and microsatellite (MS) sequences from sequencing reads.
It is specifically designed for the RETrace2 paper Figure 2 synthetic oligo library to test polymerases.

The script looks for poly-A repeats and AC dinucleotide repeats at expected positions
based on different library preparations used in the RETrace2 synthetic oligo experiments.

The script uses multiprocessing for efficient handling of large FASTQ files by splitting
them into 1,000,000 line chunks and processing them in parallel.

Usage:
    python fastq_microsatellite_extractor.py <input_fastq> <output_directory>

Arguments:
    input_fastq      Path to input FASTQ file (.fastq or .fastq.gz)
    output_directory Path to output directory where results will be saved

Example:
    python fastq_microsatellite_extractor.py data/sample.fastq.gz results/
    
Output:
    Creates a TSV file with columns: read, umi, lib, ms
    - read: Full sequencing read
    - umi: UMI sequence (positions 3-23)
    - lib: Library index sequence (positions 23-33)
    - ms: Microsatellite sequence (poly-A or AC repeats)

Supported Library Indices (RETrace2 Fig2 Synthetic Oligo Library):
    - CGCTCAGTTC: Looks for poly-A repeats at positions 97-122
    - TATCTGACCT: Looks for poly-A repeats at positions 77-107
    - ATATGAGACG: Looks for poly-A repeats at positions 83-123
    - CTTATGGAAT: Looks for AC repeats at positions 74-114
    - TAATCTCGTC: Looks for AC repeats at positions 84-134

Note: This script is specifically designed for analyzing the synthetic oligonucleotide 
library used in RETrace2 paper Figure 2 to test polymerases. The library indices and expected microsatellite 
positions are based on the experimental design described in that publication.
"""

import pandas as pd
import re
from tqdm import tqdm
import sys
import os
import subprocess
import multiprocessing
import glob

# Update 7/5/2022: add mutiprocessing by split into 1,000,000 line chunk

# Input validation
if len(sys.argv) != 3:
    print(f"Error: Expected 2 arguments, got {len(sys.argv)-1}")
    print("Usage: python fastq_microsatellite_extractor.py <input_fastq> <output_directory>")
    sys.exit(1)

in_path = sys.argv[1]   # 1st argument is input file path (.fastq or .fastq.gz)
out_dir = sys.argv[2]   # 2nd argument is output dir (return a dataframe .tsv)

# Validate input file exists
if not os.path.exists(in_path):
    print(f"Error: Input file '{in_path}' does not exist!")
    sys.exit(1)

# Validate input file format
if not (in_path.endswith('.fastq') or in_path.endswith('.fastq.gz')):
    print(f"Error: Input file must be .fastq or .fastq.gz format!")
    print(f"Provided file: {in_path}")
    sys.exit(1)

# Create output directory if it doesn't exist
if not os.path.exists(out_dir):
    print(f"Creating output directory: {out_dir}")
    os.makedirs(out_dir, exist_ok=True)

# Ensure output directory ends with slash
if not out_dir.endswith('/'):
    out_dir += '/'


def get_result_v2(path, out_dir):
    '''
    input take fastq file path and extract whole read sequence, umi, library index, and MS sequence --> output to df
    v2 also match position check within 5 repeats (+-5 or 10 bases)
    '''

    # Use list to collect results for better performance
    results_list = []
    
    n_lines = sum(1 for line in open(path))
    with open(path) as f:
        for i, line in enumerate(tqdm(f, total=n_lines)):

            # read every 2nd line for read sequence
            if i % 4 == 1:
                read = line.strip()  # CRITICAL: Remove newline character
                umi = read[3:23]
                lib = read[23:33]

                # filter library index
                # match with at least 8XA, first occurence, and +-5 base of known poistion
                if lib == 'CGCTCAGTTC':
                    ms_res = re.search('A{8}A+', read)  
                    if ms_res is not None and ms_res.span()[0] >= 97 and ms_res.span()[1] <= 122:
                        ms = ms_res.group()
                    else:
                        ms = float('NaN')

                elif lib == 'TATCTGACCT':
                    ms_res = re.search('A{8}A+', read)  
                    if ms_res is not None and ms_res.span()[0] >= 77 and ms_res.span()[1] <= 107:
                        ms = ms_res.group()
                    else:
                        ms = float('NaN')

                elif lib == 'ATATGAGACG':
                    ms_res = re.search('A{8}A+', read)  
                    if ms_res is not None and ms_res.span()[0] >= 83 and ms_res.span()[1] <= 123:
                        ms = ms_res.group()
                    else:
                        ms = float('NaN')

                # match with at least 6XAC, first occurence, and +-10 bases (5 repeat) of known poistion        
                elif lib == 'CTTATGGAAT':
                    ms_res = re.search('(?:AC){6}(?:AC)+', read)
                    if ms_res is not None and ms_res.span()[0] >= 74 and ms_res.span()[1] <= 114:
                        ms = ms_res.group()
                    else:
                        ms = float('NaN')

                elif lib == 'TAATCTCGTC':
                    ms_res = re.search('(?:AC){6}(?:AC)+', read)
                    if ms_res is not None and ms_res.span()[0] >= 84 and ms_res.span()[1] <= 134:
                        ms = ms_res.group()
                    else:
                        ms = float('NaN')

                else:
                    ms = float('NaN')  # if no match library index --> dont look for ms 

                # Append to list instead of DataFrame concatenation (much faster)
                results_list.append({'read': read, 'umi': umi, 'lib': lib, 'ms': ms})
    
    # Create DataFrame from list all at once (much more efficient)
    result = pd.DataFrame(results_list)
    
    # save
    out_name = os.path.basename(path).replace('.fastq','')
    result.to_csv(f'{out_dir}{out_name}_result.tsv', sep='\t', index=False)
    
    return(result)
    
    
print(in_path)

# create tmp dir
subprocess.call(f'mkdir {out_dir}tmp', shell=True)
tmp_dir = out_dir + "tmp/"

# split file into 1,000,000 line chunks
if in_path.endswith(".fastq.gz"):
    out_name = os.path.basename(in_path).replace('.fastq.gz','')
    subprocess.call(f'zcat {in_path} | split -l 1000000 - {tmp_dir}{out_name}_ --additional-suffix .fastq', shell=True)

elif in_path.endswith(".fastq"):
    out_name = os.path.basename(in_path).replace('.fastq','')
    subprocess.call(f'split -l 1000000 {in_path} {tmp_dir}{out_name}_ --additional-suffix .fastq', shell=True)
    
else:
    print('Error! wrong input file format!')
    
# run multiprocessing
jobs=[]
path_list = glob.glob(f'{tmp_dir}{out_name}*')
for path in path_list:
    p = multiprocessing.Process(target = get_result_v2, args = (path, tmp_dir))
    jobs.append(p)
    p.start()
        
#Join 
for p in jobs:
    p.join()

# combine all tsv output (handle headers properly)
tsv_files = glob.glob(f'{tmp_dir}{out_name}*_result.tsv')
if tsv_files:
    # Read first file with header
    combined_df = pd.read_csv(tsv_files[0], sep='\t')
    
    # Read remaining files without header and append
    for tsv_file in tsv_files[1:]:
        df = pd.read_csv(tsv_file, sep='\t')
        combined_df = pd.concat([combined_df, df], ignore_index=True)
    
    # Save final combined result
    combined_df.to_csv(f'{out_dir}{out_name}_result.tsv', sep='\t', index=False)
    print(f"Combined {len(tsv_files)} result files into final output")
else:
    print("Warning: No result files found to combine!")

# Clean up temporary directory
print("Cleaning up temporary files...")
subprocess.call(f'rm -rf {tmp_dir}', shell=True)
print(f"Processing complete! Results saved to: {out_dir}{out_name}_result.tsv")
