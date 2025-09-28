"""
RETrace2 Paper - Homopolymer Sequencing Accuracy Benchmarking Script

This script processes paired-end FASTQ files from AVITI sequencing runs to extract and analyze
microsatellite sequences for benchmarking homopolymer sequencing accuracy. Designed specifically
for the RETrace2 paper oligo model experiments (Figure 2).

The script extracts:
- UMI (Unique Molecular Identifier) sequences
- Library index sequences 
- Microsatellite sequences (homopolymer A repeats of 6+ bases)
- Full read sequences from both R1 and R2

Key features:
- Accounts for 3bp position shift due to ElementBio custom R1/R2 primers with TAA termini
- Filters for specific library indices used in the experiment
- Searches for homopolymer A sequences (minimum 6 consecutive A's) within defined position ranges
- Outputs results as TSV for downstream analysis

Usage: python script.py <R1.fastq> <R2.fastq> <output_directory>
"""

import pandas as pd
import re
from tqdm import tqdm
import sys
import os

# Input and output paths from command line arguments
in_path_r1 = sys.argv[1]   # 1st argument is input file path for read 1 (.fastq or .fastq.gz)
in_path_r2 = sys.argv[2]   # 2nd argument is input file path for read 2 (.fastq or .fastq.gz)
out_dir = sys.argv[3]      # 3rd argument is output directory (for .tsv dataframe)

def reverse_complement(dna_sequence):
    complement = dna_sequence.translate(str.maketrans('ATCG', 'TAGC'))
    return complement[::-1]

def get_result_v2(path_r1, path_r2, out_dir):
    '''
    Extracts whole read sequence, UMI, library index, and microsatellite (MS) sequence from Read1 and Read2.
    Adjusts position checks to account for initial 3bp trimming (this is special circumstance due to ElementBio use custom R1/R2 primer with TAA in the end).
    V2.1 update: make filter less stringent. Regex to find min 6xA. Poistion check allow +-10bp(repeats).
    '''
    result = []

    n_lines_r1 = sum(1 for line in open(path_r1, 'rt'))
    n_lines_r2 = sum(1 for line in open(path_r2, 'rt'))
    assert n_lines_r1 == n_lines_r2, "Read 1 and Read 2 files have different number of lines"

    with open(path_r1, 'rt') as f_r1, open(path_r2, 'rt') as f_r2:
        for i, (line_r1, line_r2) in enumerate(zip(tqdm(f_r1, total=n_lines_r1), tqdm(f_r2, total=n_lines_r2))):

            # Process every 2nd line for read sequence
            if i % 4 == 1:
                read_r1 = line_r1.rstrip()
                read_r2 = reverse_complement(line_r2.rstrip())
                
                # Adjusted UMI and lib extraction positions to account for initial 3bp trimming
                umi = read_r1[0:20]  # Adjusted: Original was [3:23], now starts from 0 due to 3bp trimming
                lib = read_r1[20:30]  # Adjusted: Original was [23:33]

                # Filter library index based on Read 1 and then extract microsatellite sequence from Read 1 & 2
                # Adjusted position checks by subtracting 3 from the start and end positions
                ms_r1 = ms_r2 = float('NaN')  # Default values
                if lib in ['CGCTCAGTTC', 'TATCTGACCT', 'ATATGAGACG', 'GCGCGATGTT', 'AGAGCACTAG']:
                    # Define the adjusted position ranges for each library
                    positions = {
                        'CGCTCAGTTC': (89, 124),
                        'TATCTGACCT': (69, 109),
                        'ATATGAGACG': (75, 125),
                        'GCGCGATGTT': (95, 125),
                        'AGAGCACTAG': (80, 125)
                    }

                    start_pos, end_pos = positions[lib]

                    # Process Read 1
                    ms_res = re.search('A{5}A+', read_r1)
                    if ms_res and start_pos <= ms_res.span()[0] <= end_pos:
                        ms_r1 = ms_res.group()

                    # Process Read 2
                    ms_res = re.search('A{5}A+', read_r2)
                    if ms_res and start_pos <= ms_res.span()[0] <= end_pos:
                        ms_r2 = ms_res.group()

                # Save results
                new_result = {'read_r1': read_r1, 'read_r2': read_r2, 'umi': umi, 'lib': lib, 'ms_r1': ms_r1, 'ms_r2': ms_r2}
                result.append(new_result)
    
    # Output the results to a DataFrame and save as TSV
    result_df = pd.DataFrame(result)
    out_name = os.path.basename(path_r1).replace('.fastq','').replace('.gz', '')
    result_df.to_csv(f'{out_dir}/{out_name}_R1R2result.tsv', sep='\t', index=False)
    
    return result_df

print(f"Processing files:\nRead 1: {in_path_r1}\nRead 2: {in_path_r2}\nOutput directory: {out_dir}")
get_result_v2(in_path_r1, in_path_r2, out_dir)



