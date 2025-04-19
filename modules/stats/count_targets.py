#!/usr/bin/env python3

import os
import argparse
import pysam
from collections import defaultdict



class MicrosatelliteCounter:
    """
    Process BAM files to count unique microsatellite targets.
    """
    
    def __init__(self, target_bed=None, min_depth=30):
        """Initialize counter with target regions and minimum read depth."""
        self.target_bed = target_bed
        self.min_depth = min_depth
        self.targets = self._load_targets() if target_bed else None
    
    def _load_targets(self):
        """Load target regions from BED file."""
        targets = {}
        with open(self.target_bed, 'r') as f:
            for line in f:
                fields = line.strip().split()
                if len(fields) >= 3:
                    chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                    target_id = f"{chrom}:{start}-{end}"
                    targets[target_id] = (chrom, start, end)
        return targets
    
    def count_targets(self, bam_file):
        """
        Count targets in a BAM file with coverage >= min_depth.
        Currently only support single-end reads counting.
        """
        depth_dict = defaultdict(int)
        
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            # For each target in the BED file
            for target_id, (chrom, start, end) in self.targets.items():
                # Get all reads that overlap this target region
                for read in bam.fetch(chrom, start, end):
                    # Check if read fully spans the target region
                    if read.reference_start <= start and read.reference_end >= end:
                        depth_dict[target_id] += 1
        
        return depth_dict
    
    def write_stats(self, targets_dict, output_file):
        """Write statistics to output file."""
        
        with open(output_file, 'w') as f_out:
            f_out.write("target_id\tdepth\tcovered_${self.min_depth}X\n")
            
            for target_id, depth in sorted(targets_dict.items()):
                covered = "True" if depth >= self.min_depth else "False"
                f_out.write(f"{target_id}\t{depth}\t{covered}\n")
        


def main():
    parser = argparse.ArgumentParser(description="Count unique microsatellite targets in a BAM file")
    parser.add_argument("--bam", required=True, help="Input BAM file (must be indexed)")
    parser.add_argument("--target_bed", required=True, help="BED file containing target microsatellite locations")
    parser.add_argument("--min_depth", type=int, default=30, help="Minimum read depth for a target to be counted")
    parser.add_argument("--output", required=True, help="Output filename for results")
    
    args = parser.parse_args()
    
    # Check if target_bed exists
    if not os.path.exists(args.target_bed):
        raise FileNotFoundError(f"Target BED file not found: {args.target_bed}")
        
    # Verify BAM file is indexed
    if not os.path.exists(args.bam + ".bai"):
        print(f"BAM index not found for {args.bam}. Creating index...")
        pysam.index(args.bam)
    
    counter = MicrosatelliteCounter(args.target_bed, args.min_depth)
    
    # Process BAM file
    print(f"Processing BAM file: {args.bam}")
    targets_dict = counter.count_targets(args.bam)
    
    # Write results
    counter.write_stats(targets_dict, args.output)
    print(f"Results written to {args.output}")
    print("Analysis complete.")


if __name__ == "__main__":
    main()
