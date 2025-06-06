#!/bin/bash
set -euo pipefail

# Description: Process single-end FASTQ files matching a wildcard pattern (*R1*fastq.gz) to map to lambda DNA for bisulfite conversion rate calculation, saving summary statistics to summary.txt.
# Note: Ensure unmethylated lambda DNA is spiked into your library for accurate bisulfite conversion rate analysis.
# Usage: ./lambda_bisulfite.sh -i <input_pattern> -o <output_dir> -g <lambda_genome> -b <bismark_path> [-t <threads>]
# Example: ./lambda_bisulfite.sh -i "path/to/Methyl*R1*fastq.gz" -o path/to/output -g path/to/lambda -b /path/to/bismark -t 4

# Defaults
THREADS=4

# Usage function
usage() {
    echo "Usage: $0 -i <input_pattern> -o <output_dir> -g <lambda_genome> -b <bismark_path> [-t <threads>]"
    exit 1
}

# Parse arguments
while getopts "i:o:g:b:t:h" opt; do
    case $opt in
        i) INPUT_PATTERN="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        g) GENOME_PATH="$OPTARG" ;;
        b) BISMARK_PATH="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        h) usage ;;
    esac
done

# Validate inputs
[[ -z "${INPUT_PATTERN:-}" || -z "${OUTPUT_DIR:-}" || -z "${GENOME_PATH:-}" || -z "${BISMARK_PATH:-}" ]] && usage
[[ ! -d "$GENOME_PATH" ]] && { echo "Error: Genome path $GENOME_PATH not a directory"; exit 1; }
[[ ! -x "$BISMARK_PATH" ]] && { echo "Error: Bismark executable $BISMARK_PATH not found or not executable"; exit 1; }

# Create output directories
mkdir -p "$OUTPUT_DIR/trimmed" "$OUTPUT_DIR/bam" || { echo "Error: Failed to create output directories"; exit 1; }

# Summary file
SUMMARY_FILE="$OUTPUT_DIR/summary.txt"
> "$SUMMARY_FILE" # Clear or create summary file

# Find R1 files matching the pattern
mapfile -t R1_FILES < <(ls $INPUT_PATTERN 2>/dev/null | grep -E "R1.*fastq\.gz$")
if [[ ${#R1_FILES[@]} -eq 0 ]]; then
    echo "Error: No R1 FASTQ files found matching pattern '$INPUT_PATTERN'"
    ls $(dirname "$INPUT_PATTERN") 2>/dev/null || echo "Directory $(dirname "$INPUT_PATTERN") not accessible"
    exit 1
fi

echo "Found ${#R1_FILES[@]} R1 files:"
printf '%s\n' "${R1_FILES[@]}"

# Arrays to store stats for averaging
declare -a mapping_efficiencies
declare -a c_methylated_cpg

# Process FASTQ files
for r1 in "${R1_FILES[@]}"; do
    # Extract base name by removing 'fastq.gz'
    fbname=$(basename "$r1" .fastq.gz)
    # Sample name: take prefix before 'R1'
    sample=$(echo "$fbname" | sed 's/R1.*//')

    echo "Processing sample: $sample"

    # Trim reads (single-end)
    echo "Running Trim Galore for $sample"
    if ! trim_galore --quality 20 --stringency 3 --length 36 --rrbs --cores "$THREADS" \
        --output_dir "$OUTPUT_DIR/trimmed" "$r1"; then
        echo "Error: Trim Galore failed for $sample"
        continue
    fi

    # Verify trimmed file exists
    trimmed_r1="$OUTPUT_DIR/trimmed/${fbname}_trimmed.fq.gz"
    [[ ! -f "$trimmed_r1" ]] && { echo "Error: Trimmed file missing for $sample"; continue; }

    # Align with Bismark (single-end)
    echo "Running Bismark for $sample"
    if ! "$BISMARK_PATH" --genome "$GENOME_PATH" \
        "$trimmed_r1" \
        --output_dir "$OUTPUT_DIR/bam" -p "$THREADS"; then
        echo "Error: Bismark failed for $sample"
        continue
    fi

    # Extract stats from Bismark report
    report_file="$OUTPUT_DIR/bam/${fbname}_trimmed_bismark_bt2_SE_report.txt"
    if [[ -f "$report_file" ]]; then
        # Extract mapping efficiency
        mapping_eff=$(grep "Mapping efficiency:" "$report_file" | awk '{print $3}' | tr -d '%')
        # Extract C methylated in CpG context
        c_methyl_cpg=$(grep "C methylated in CpG context:" "$report_file" | awk '{print $6}' | tr -d '%')
        # Calculate bisulfite conversion rate
        bisulfite_conv=$(echo "100 - $c_methyl_cpg" | bc)

        # Store for averaging
        mapping_efficiencies+=("$mapping_eff")
        c_methylated_cpg+=("$c_methyl_cpg")

        # Save stats to summary file
        {
            echo "Summary for $sample:"
            echo "  Mapping efficiency: $mapping_eff% (percentage of reads mapping to lambda DNA)"
            echo "  C methylated in CpG context: $c_methyl_cpg%"
            echo "  Bisulfite conversion rate: $bisulfite_conv% (100 - C methylated in CpG context)"
            echo ""
        } >> "$SUMMARY_FILE"
    else
        echo "Error: Bismark report $report_file not found for $sample"
    fi

    # Clean up
    echo "Cleaning up for $sample"
    rm -f "$trimmed_r1" || echo "Warning: Failed to clean up trimmed file for $sample"
done

# Calculate and save average stats
if [[ ${#mapping_efficiencies[@]} -gt 0 ]]; then
    {
        echo "Summary Statistics Across All Files:"
        echo "Note: Mapping efficiency represents the percentage of reads mapping to lambda DNA, which should approximate your experimental spike-in rate."
        echo "      Since this is unmethylated lambda DNA, bisulfite conversion rate is 100 - C methylated in CpG context percentage."
        echo ""
    } >> "$SUMMARY_FILE"

    # Calculate averages
    total_mapping=0
    total_c_methyl_cpg=0
    count=${#mapping_efficiencies[@]}

    for eff in "${mapping_efficiencies[@]}"; do
        total_mapping=$(echo "$total_mapping + $eff" | bc)
    done
    for cpg in "${c_methylated_cpg[@]}"; do
        total_c_methyl_cpg=$(echo "$total_c_methyl_cpg + $cpg" | bc)
    done

    avg_mapping=$(echo "scale=2; $total_mapping / $count" | bc)
    avg_c_methyl_cpg=$(echo "scale=2; $total_c_methyl_cpg / $count" | bc)
    avg_bisulfite_conv=$(echo "scale=2; 100 - $avg_c_methyl_cpg" | bc)

    {
        echo "Average Mapping efficiency: $avg_mapping%"
        echo "Average C methylated in CpG context: $avg_c_methyl_cpg%"
        echo "Average Bisulfite conversion rate: $avg_bisulfite_conv%"
    } >> "$SUMMARY_FILE"
else
    echo "No files successfully processed for summary statistics"
    echo "No files successfully processed" >> "$SUMMARY_FILE"
fi

echo "Done"