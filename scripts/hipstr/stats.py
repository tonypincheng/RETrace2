#!/usr/bin/env python3

import argparse
import sys
from collections import Counter
import os

def parse_args():
    parser = argparse.ArgumentParser(description='Generate statistics from HipSTR VCF file')
    parser.add_argument('vcf', help='Input VCF file')
    return parser.parse_args()

def parse_vcf(vcf_file):
    """Parse VCF file and extract metrics"""
    n_variants = 0
    n_samples = 0
    call_rate = 0
    allele_counts = Counter()
    sample_counts = {}
    motif_lengths = Counter()
    
    with open(vcf_file, 'r') as f:
        for line in f:
            # Skip header lines
            if line.startswith('##'):
                continue
            
            # Extract sample names from header
            if line.startswith('#CHROM'):
                fields = line.strip().split('\t')
                n_samples = len(fields) - 9  # Subtract standard VCF fields
                sample_counts = {sample: 0 for sample in fields[9:]}
                continue
            
            # Process variant lines
            n_variants += 1
            fields = line.strip().split('\t')
            
            # Extract variant info
            info = dict(item.split('=') if '=' in item else (item, True) 
                        for item in fields[7].split(';'))
            
            # Count by repeat motif length
            motif = info.get('MOTIF', 'N/A')
            motif_lengths[len(motif)] += 1
            
            # Process genotypes
            n_called = 0
            for i, sample in enumerate(sample_counts.keys()):
                gt = fields[9 + i].split(':')[0]
                
                # Skip missing genotypes
                if gt == './.':
                    continue
                
                n_called += 1
                sample_counts[sample] += 1
                
                # Count alleles
                alleles = gt.replace('|', '/').split('/')
                for allele in alleles:
                    if allele != '.':
                        allele_counts[allele] += 1
            
            # Calculate call rate for this variant
            call_rate += n_called / n_samples
    
    # Calculate average call rate
    if n_variants > 0:
        call_rate /= n_variants
    
    return {
        'n_variants': n_variants,
        'n_samples': n_samples,
        'call_rate': call_rate,
        'allele_counts': allele_counts,
        'sample_counts': sample_counts,
        'motif_lengths': motif_lengths
    }

def format_stats(stats):
    """Format statistics for output"""
    lines = []
    lines.append(f"Total variants: {stats['n_variants']}")
    lines.append(f"Total samples: {stats['n_samples']}")
    lines.append(f"Average call rate: {stats['call_rate']:.4f}")
    
    lines.append("\nVariants by motif length:")
    for length, count in sorted(stats['motif_lengths'].items()):
        lines.append(f"  {length} bp: {count} ({count/stats['n_variants']:.2%})")
    
    lines.append("\nTop 5 alleles:")
    for allele, count in stats['allele_counts'].most_common(5):
        lines.append(f"  Allele {allele}: {count}")
    
    lines.append("\nCall rate by sample:")
    for sample, count in stats['sample_counts'].items():
        rate = count / stats['n_variants'] if stats['n_variants'] > 0 else 0
        lines.append(f"  {sample}: {rate:.4f} ({count}/{stats['n_variants']})")
    
    return '\n'.join(lines)

def main():
    args = parse_args()
    
    if not os.path.exists(args.vcf):
        print(f"Error: VCF file {args.vcf} not found", file=sys.stderr)
        sys.exit(1)
    
    try:
        stats = parse_vcf(args.vcf)
        print(format_stats(stats))
    except Exception as e:
        print(f"Error processing VCF file: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main() 