#!/usr/bin/env python3

import argparse
import gzip
import sys


def count_cpgs(allc_file):
    """Count CpGs from allc file"""
    counts = {
        'total_cpgs': 0,
        'methylated_cpgs': 0,
        'unmethylated_cpgs': 0,
        'methylation_rate': 0.0
    }
    
    # Check if the file is gzipped
    open_func = gzip.open if allc_file.endswith('.gz') else open
    mode = 'rt' if allc_file.endswith('.gz') else 'r'
    
    try:
        with open_func(allc_file, mode) as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < 7:
                    continue
                    
                # Check if it's a CpG context (starts with 'CG')
                context = fields[3]
                if context.startswith('CG'):
                    counts['total_cpgs'] += 1
                    meth_level = int(fields[4])
                    unmeth_level = int(fields[5])
                    
                    if meth_level > 0:
                        counts['methylated_cpgs'] += 1
                    else:
                        counts['unmethylated_cpgs'] += 1
    except Exception as e:
        print(f"Error processing file {allc_file}: {e}", file=sys.stderr)
        return counts
    
    # Calculate methylation rate
    if counts['total_cpgs'] > 0:
        counts['methylation_rate'] = counts['methylated_cpgs'] / counts['total_cpgs']
    
    return counts

def write_output(counts, output_file):
    """Write counts to output file"""
    with open(output_file, 'w') as f:
        for key, value in counts.items():
            f.write(f"{key}\t{value}\n")

def main():
    parser = argparse.ArgumentParser(description='Count CpGs from allc files')
    parser.add_argument('--allc', required=True, help='allc file from methylpy')
    parser.add_argument('--output', required=True, help='Output file to write CpG counts')
    
    args = parser.parse_args()
    
    # Count CpGs
    print(f"Counting CpGs in {args.allc}...")
    counts = count_cpgs(args.allc)
    
    # Write output
    print(f"Writing CpG counts to {args.output}")
    write_output(counts, args.output)
    
    print("Done!")

if __name__ == "__main__":
    main() 