# Reference File Templates

This directory contains template files that are used by the RETrace2 pipeline.

## Reference Genome Management

For large reference genomes (such as mm39 or hg38), we recommend the following approach:

1. **Create a central reference repository:**
   - Store reference files in a central location accessible to all users
   - This avoids duplicating large files across multiple pipeline installations

2. **Use symbolic links:**
   - Create a `references/` directory in your pipeline installation
   - Use symbolic links to point to the central repository

3. **Automated download:**
   - The pipeline includes scripts to automatically download reference genomes
   - Run `bin/download_references.sh` to download and prepare reference genomes

## Recommended Reference Locations

### Option 1: Institution-wide Reference Repository
```
/path/to/shared/references/
├── mm39/               # Mouse genome (default)
│   ├── genome.fa
│   ├── genome.fa.fai
│   ├── bwa/
│   │   ├── mm39.amb
│   │   ├── mm39.ann
│   │   ├── mm39.bwt
│   │   ├── mm39.pac
│   │   └── mm39.sa
│   └── methylpl/
│       └── [methylpl index files]
├── hg38/
│   ├── genome.fa
│   ├── genome.fa.fai
│   ├── bwa/
│   │   ├── hg38.amb
│   │   ├── hg38.ann
│   │   ├── hg38.bwt
│   │   ├── hg38.pac
│   │   └── hg38.sa
│   └── methylpl/
│       └── [methylpl index files]
```

Update `nextflow.config` with:
```
params.genomes_base = "/path/to/shared/references"
```

### Option 2: Per-Project References
```
your-project/
├── RETrace2/
│   └── references/   # Symbolic links to actual reference files
└── reference-data/   # Actual reference files
```

## Directory Structure

Each reference genome should follow this structure:
```
genome_id/
├── genome.fa         # Reference genome FASTA
├── genome.fa.fai     # FASTA index
├── bwa/              # BWA indices
│   ├── genome_id.amb
│   ├── genome_id.ann
│   ├── genome_id.bwt
│   ├── genome_id.pac
│   └── genome_id.sa
└── methylpl/         # Methylpl indices
    └── [methylpl index files]
```

## Reference Download and Indexing

To download and prepare the mm39 reference genome:

```bash
# Download and index mm39 (default)
./bin/download_references.sh

# Download and index hg38
./bin/download_references.sh -g hg38
```

The script will:
1. Download the reference genome FASTA
2. Create BWA index
3. Create samtools index (faidx)
4. Create methylpl index 

# RETrace2

> ⚠️ **Under Active Development** ⚠️
> This pipeline is currently under development. Features and documentation may change. 