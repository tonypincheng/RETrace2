#!/usr/bin/env python3
"""
parse_vcf.py - Parse HipSTR VCF output to extract microsatellite alleles

Output Dictionary Structure:
{
    "target_id": {                     # Target identifier from BED file
        "motif_len": int,              # Length of the repeat motif
        "num_copies": int,             # Number of repeat copies in reference
        "sample": {
            "sample_name": {           # Sample identifier from VCF
                "msCount": [int, ...], # List of observed microsatellite counts in base pairs diffrence from reference
                "allelotype": [int, int]  # Allele sizes in base pairs (not repeat units)
            },
            ...
        }
    },
    ...
}

This dictionary will be used for phylogenetic tree reconstruction.
"""

import pickle
import argparse

def load_target_bed(bed_file):
    """Load target information from BED file"""
    target_dict = {}
    with open(bed_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 6:
                continue
            
            chrom = fields[0]
            start = fields[1]
            end = fields[2]
            
            # Extract motif length, number of copies, and target ID directly from columns
            try:
                motif_len = int(fields[3])
                num_copies = int(fields[4])
                target_id = fields[5] 
                
                target_dict[target_id] = {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "motif_len": motif_len,    # Store motif length
                    "num_copies": num_copies   # Store number of copies
                }
            except (ValueError, IndexError):
                # Skip if we can't parse the fields correctly
                continue
                
    return target_dict

def parseVCF(HipSTR_vcf, target_bed, min_qual, min_reads, max_stutter):
    print("Extracting allelotype from HipSTR output")
    alleleDict = {}
    
    # Load target information from BED file
    print(f"Loading target information from {target_bed}")
    targetDict = load_target_bed(target_bed)
    print(f"Loaded {len(targetDict)} targets from BED file")
    
    with open(HipSTR_vcf) as vcf:
        for line in vcf:
            if not line.startswith('##'):
                if line.startswith('#'):
                    info_line = line.split()
                    sample_names = info_line[9:]
                else:
                    vcf_line = line.split()
                    # Extract target ID from VCF - use full target ID
                    target_id = vcf_line[2]
                    
                    if target_id not in targetDict:
                        print(f"Target {target_id} not found in target BED. Skipping...")
                        continue
                    
                    # Get reference information from the target dictionary
                    target_info = targetDict[target_id]
                    motif_len = target_info["motif_len"]
                    num_copies = target_info["num_copies"]
                    ref_MS_len = motif_len * num_copies
                    
                    alleleDict[target_id] = {
                        "motif_len": motif_len,
                        "num_copies": num_copies,
                        "sample": {}
                    }
                    
                    sample_format = vcf_line[9:]
                    for indx, val in enumerate(sample_format):
                        sample = sample_names[indx]
                        if val != '.' and ':' in val:  # Make sure entry isn't missing and has proper format
                            fields = val.split(':')
                            # Check if we have enough fields and none are '.'
                            if len(fields) > 7:  # Make sure we have all needed fields
                                allelotype = fields[1]
                                prob_genotype = fields[2]
                                num_reads = fields[4]
                                num_stutter = fields[6]
                                msCounts = fields[-2]
                                
                                # Check if any required fields are '.' before conversion
                                if all(x != '.' for x in [prob_genotype, num_reads, num_stutter]):
                                    try:
                                        prob_genotype_float = float(prob_genotype)
                                        num_reads_int = int(num_reads)
                                        num_stutter_int = int(num_stutter)
                                        
                                        if (prob_genotype_float >= min_qual and 
                                            num_reads_int >= min_reads and 
                                            (num_stutter_int/num_reads_int) <= max_stutter and 
                                            msCounts != '.'):
                                            
                                            alleleDict[target_id]["sample"][sample] = {}
                                            alleleDict[target_id]["sample"][sample]["msCount"] = []
                                            for msCount_info in msCounts.split(';'):
                                                [msCount, freq] = msCount_info.split('|')
                                                alleleDict[target_id]["sample"][sample]["msCount"].extend([int(msCount)] * int(freq))                                            
                                            # Add allelotype info which is shown as number of bases not number of repeats
                                            alleleDict[target_id]["sample"][sample]["allelotype"] = [int(allele) + ref_MS_len for allele in allelotype.split('|')]
                                    except (ValueError, ZeroDivisionError):
                                        continue

    return alleleDict, sample_names

def main():
    parser = argparse.ArgumentParser(description='Parse HipSTR VCF output to extract microsatellite alleles')
    parser.add_argument('--vcf', required=True, help='HipSTR VCF file')
    parser.add_argument('--target_bed', required=True, help='Target BED file')
    parser.add_argument('--output_pkl', required=True, help='Output pickle file for alleleDict')
    parser.add_argument('--output_samples', required=True, help='Output file for sample list')
    parser.add_argument('--min_qual', type=float, default=0.9, help='Minimum quality score (default: 0.9)')
    parser.add_argument('--min_reads', type=int, default=10, help='Minimum number of reads (default: 10)')
    parser.add_argument('--max_stutter', type=float, default=1, help='Maximum stutter ratio (default: 1)')
    
    args = parser.parse_args()
    
    # Parse VCF
    alleleDict, sample_names = parseVCF(
        args.vcf,
        args.target_bed,
        args.min_qual, 
        args.min_reads, 
        args.max_stutter
    )
    
    # Save alleleDict to pickle file
    with open(args.output_pkl, "wb") as f:
        pickle.dump(alleleDict, f)
    
    # Save sample list to a separate file
    with open(args.output_samples, "w") as f:
        f.write("\n".join(sample_names))
    
    print(f"Saved allele dictionary to {args.output_pkl}")
    print(f"Saved sample list to {args.output_samples}")

if __name__ == "__main__":
    main() 