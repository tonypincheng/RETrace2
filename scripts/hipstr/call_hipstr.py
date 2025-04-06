#!/usr/bin/env python3

import argparse
import subprocess
import os
import sys
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description='Call microsatellites using HipSTR')
    parser.add_argument('--bam', required=True, help='Input BAM file')
    parser.add_argument('--reference', required=True, help='Reference genome FASTA')
    parser.add_argument('--regions', required=True, help='Regions BED file')
    parser.add_argument('--output', required=True, help='Output VCF file')
    parser.add_argument('--min-reads', type=int, default=10, help='Minimum number of reads required')
    parser.add_argument('--min-allele-freq', type=float, default=0.1, help='Minimum allele frequency')
    return parser.parse_args()

def call_hipstr(args):
    cmd = [
        'hipstr',
        '--bam', args.bam,
        '--fasta', args.reference,
        '--regions', args.regions,
        '--str-vcf', args.output,
        '--min-reads', str(args.min_reads),
        '--min-allele-freq', str(args.min_allele_freq),
        '--log', f'{args.output}.log'
    ]
    
    try:
        subprocess.run(cmd, check=True)
        print(f"Successfully generated {args.output}")
    except subprocess.CalledProcessError as e:
        print(f"Error running HipSTR: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    args = parse_args()
    call_hipstr(args)

if __name__ == '__main__':
    main() 