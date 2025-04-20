#!/usr/bin/env python3

import os
import sys
import errno
import argparse
import csv

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Check and validate samplesheet format.')
    parser.add_argument('--samplesheet', required=True, help='Input samplesheet file.')
    parser.add_argument('--output', required=True, help='Output file for validated samplesheet.')
    
    return parser.parse_args()

def check_samplesheet(file_in, file_out):
    """
    Check that the samplesheet follows the required format.
    """
    required_headers = ["sample_id", "ms_fastq_1"]
    optional_headers = ["meth_fastq_1", "group", "color"]
    
    # Initialize sample tracking
    sample_ids = set()
    
    # Write validated samplesheet
    with open(file_in, "r") as fin, open(file_out, "w") as fout:
        reader = csv.DictReader(fin)
        header = reader.fieldnames
        
        # Check header
        if not header:
            print("ERROR: Samplesheet header is empty or malformed.", file=sys.stderr)
            sys.exit(1)
            
        # Check required headers
        missing_headers = set(required_headers) - set(header)
        if missing_headers:
            print(f"ERROR: Missing required column headers: {', '.join(missing_headers)}", file=sys.stderr)
            sys.exit(1)
        
        # Add optional headers that are present
        final_headers = required_headers.copy()
        for h in optional_headers:
            if h in header:
                final_headers.append(h)
        
        # Write header
        writer = csv.DictWriter(fout, fieldnames=final_headers)
        writer.writeheader()
        
        for row in reader:
            # Skip comment lines
            if row["sample_id"].startswith("#"):
                continue
                
            # Check for empty sample_id
            if not row["sample_id"]:
                print("ERROR: Sample ID is empty.", file=sys.stderr)
                sys.exit(1)
                
            # Check for duplicate sample_id
            if row["sample_id"] in sample_ids:
                print(f"ERROR: Duplicate sample ID found: {row['sample_id']}", file=sys.stderr)
                sys.exit(1)
            sample_ids.add(row["sample_id"])
            
            # Check ms_fastq_1 exists
            if not row["ms_fastq_1"]:
                print(f"ERROR: MS FASTQ file path is empty for sample: {row['sample_id']}", file=sys.stderr)
                sys.exit(1)
                
            # Create a new row with only the final headers
            new_row = {h: row.get(h, "") for h in final_headers}
            writer.writerow(new_row)
    
    # Check if samplesheet is empty
    if not sample_ids:
        print("ERROR: Samplesheet contains no valid data.", file=sys.stderr)
        sys.exit(1)
        
    print(f"INFO: Processed {len(sample_ids)} samples.")
    
def main():
    args = parse_args()
    
    # Check input file exists
    if not os.path.exists(args.samplesheet):
        print(f"ERROR: Samplesheet file does not exist: {args.samplesheet}", file=sys.stderr)
        sys.exit(1)
    
    # Make sure the output directory exists
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    
    # Check and validate samplesheet
    check_samplesheet(args.samplesheet, args.output)

if __name__ == "__main__":
    main() 