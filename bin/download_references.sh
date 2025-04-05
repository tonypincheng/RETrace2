#!/usr/bin/env bash
# Script to download and prepare reference genomes for RETrace2

set -e

# Default settings
REF_DIR="references"
GENOME="mm39"
BUILD_BWA_INDEX=true
BUILD_SAMTOOLS_INDEX=true
BUILD_METHYLPL_INDEX=true

# Help message
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Download and prepare reference genomes for RETrace2."
    echo ""
    echo "Options:"
    echo "  -g, --genome GENOME     Specify genome to download (mm39, hg38)"
    echo "  -d, --directory DIR     Directory to store references (default: references)"
    echo "  --no-bwa                Don't build BWA index"
    echo "  --no-samtools           Don't build samtools index"
    echo "  --no-methylpl           Don't build methylpl index"
    echo "  -h, --help              Show this help message"
    exit 1
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -g|--genome)
            GENOME="$2"
            shift 2
            ;;
        -d|--directory)
            REF_DIR="$2"
            shift 2
            ;;
        --no-bwa)
            BUILD_BWA_INDEX=false
            shift
            ;;
        --no-samtools)
            BUILD_SAMTOOLS_INDEX=false
            shift
            ;;
        --no-methylpl)
            BUILD_METHYLPL_INDEX=false
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Create reference directory
mkdir -p "${REF_DIR}/${GENOME}/bwa" "${REF_DIR}/${GENOME}/methylpl"

echo "=== Downloading reference genome: ${GENOME} ==="

# Download the appropriate genome
case ${GENOME} in
    mm39)
        # UCSC mm39
        echo "Downloading mm39 genome..."
        wget -O "${REF_DIR}/${GENOME}/genome.fa.gz" "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz"
        gunzip -c "${REF_DIR}/${GENOME}/genome.fa.gz" > "${REF_DIR}/${GENOME}/genome.fa"
        ;;
        
    hg38)
        # UCSC hg38
        echo "Downloading hg38 genome..."
        wget -O "${REF_DIR}/${GENOME}/genome.fa.gz" "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
        gunzip -c "${REF_DIR}/${GENOME}/genome.fa.gz" > "${REF_DIR}/${GENOME}/genome.fa"
        ;;
        
    *)
        echo "Error: Unsupported genome ${GENOME}"
        usage
        ;;
esac

# Build BWA index if requested
if [ "$BUILD_BWA_INDEX" = true ]; then
    echo "Building BWA index..."
    bwa index -p "${REF_DIR}/${GENOME}/bwa/${GENOME}" "${REF_DIR}/${GENOME}/genome.fa"
fi

# Build samtools index if requested
if [ "$BUILD_SAMTOOLS_INDEX" = true ]; then
    echo "Building samtools index..."
    samtools faidx "${REF_DIR}/${GENOME}/genome.fa"
fi

# Build methylpl index if requested
if [ "$BUILD_METHYLPL_INDEX" = true ]; then
    echo "Building methylpl index..."
    if [ -x "$(command -v methylpl)" ]; then
        mkdir -p "${REF_DIR}/${GENOME}/methylpl"
        methylpl index -r "${REF_DIR}/${GENOME}/genome.fa" -o "${REF_DIR}/${GENOME}/methylpl"
    else
        echo "Warning: methylpl is not installed. Skipping methylpl index generation."
        echo "To build the methylpl index, install methylpl and run:"
        echo "  methylpl index -r ${REF_DIR}/${GENOME}/genome.fa -o ${REF_DIR}/${GENOME}/methylpl"
    fi
fi

echo "Reference preparation complete for ${GENOME}!"
echo "Reference location: ${REF_DIR}/${GENOME}"
echo ""
echo "Update your nextflow.config file with:"
echo "params.genomes_base = '$(realpath ${REF_DIR})'"

# Make script executable
chmod +x "$0" 