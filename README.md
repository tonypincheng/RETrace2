# RETrace2

> ⚠️ **Under Active Development** ⚠️

> This pipeline is currently under development. Features and documentation may change.

A Nextflow-based bioinformatics pipeline for RETrace2 analysis.

## Overview

RETrace2 is a comprehensive pipeline for processing and analyzing sequencing data to reconstruct single-cell phylogenetic trees using somatic microsatellite mutations. It integrates quality control, alignment, microsatellite genotype calling, and tree reconstruction into a streamlined workflow.
<br>


## Features

- Core pipeline:
  - FastQC and Adapter trimming
  - Alignment with BWA
  - Microsatellite calling with HipSTR
  - Phylogenetic tree reconstruction

- Optional analyses:
  - Tree bootstrapping
  - Tree accuracy evaluation (with ground truth)
  - Methylation analysis for cell type inference
<br>


## Pipeline Structure

```
RETrace2/
├── bin/                    # Utility scripts
├── modules/               # Nextflow modules
│   ├── mapping/          # Read mapping and alignment
│   ├── hipstr/          # Microsatellite calling
│   ├── phylo/           # Phylogenetic tree reconstruction
│   ├── bootstrap/       # Tree bootstrapping analysis
│   ├── methylation/     # Methylation analysis
│   └── evaluation/      # Tree accuracy evaluation
├── scripts/              # Supporting scripts
├── main.nf               # Main pipeline
└── nextflow.config      # Pipeline configuration
```
<br>


## Test Data

The repository includes small test datasets in the `data/` directory:

### HCT116 Human Cell Line Data
- 6 FASTQ files with ~100,000 reads each
- HCT116 cell line clones from a cell culture tree model
- 12k probe set targeting homopolymers (10-14bp repeat lengths)
- Includes ground truth data for validation
- Sufficient data to test the pipeline's basic functionality and ground truth validation

### MSH2 Mouse Data
- 6 FASTQ files with ~100,000 reads each
- Includes both microsatellite and methylation libraries
- Designed for testing dual-omic pipeline capabilities
- Sufficient data to test the pipeline's methylation analysis functionality

For details about the test data, see the [data/README.md](data/README.md) file.
<br><br>


## Getting Started

### Prerequisites
- nextflow=24.10.5
- fastqc=0.12.1
- multiqc=1.28
- python=3.9
- trim-galore=0.6.10
- methylpy=1.4.7
- bwa=0.7.19

> **Note:** You can also use Docker or Conda (see below) to handle dependencies automatically.

### Installation
```bash
git clone https://github.com/tonypincheng/retrace2.git
cd retrace2
```

### Example Usage
```bash
# Run core pipeline
nextflow run main.nf \
          --input_dir data/MSH2 \
          --fastq_pattern = "MS*.fastq.gz"
          --output_dir results/ \
          --genomes_base /path/to/reference/genomes \
          --target_bed resources/targets/mm39/RETrace2.mm39.1nt10-30bp.92460targets169818probes.bed
```

```bash
# To see all available option 
nextflow run main.nf --help
```
<br>


## Environments
RETrace2 supports multiple execution environments through Nextflow profiles.

### Standard (Default)
By default, the pipeline runs using your system's native tools without Docker or Conda.

```bash
# This runs using your locally installed packages
nextflow run main.nf
```

You can install the required packages from the environment.yml file:

```bash
# Install directly from the environment.yml file
conda env create -f environment.yml

# Activate the environment
conda activate retrace2
```

Or install packages individually according to the versions specified in environment.yml.

### Docker
To use Docker containers for all tools:

```bash
nextflow run main.nf -profile docker
```

If you don't have Docker installed, get it from the [official website](https://docs.docker.com/get-docker/).

### Conda
Alternatively, you can use Conda to automatically create environments with required dependencies:

```bash
nextflow run main.nf -profile conda
```

This will automatically create and manage Conda environments based on the requirements in `environment.yml`.
<br>


## Reference Genome Configuration

RETrace2 provides flexible options for specifying reference genomes:

### Option 1: Standard Directory Structure

The pipeline uses a standardized directory structure for organizing reference genomes:

```
genomes_base/
├── mm39/                  # Genome name
│   ├── raw_fasta/         
│   │   └── mm39.fa        # raw fasta files
│   │   └── mm39.fa.fai          
│   ├── bwa-index/         # BWA index files
│   │   └── mm39.fa        
│   └── methylpl-ref/      # Methylation references
│       └── mm39_f         # Forward methylation reference
│       └── mm39_r         # Reverse methylation reference
├── hg38/
│   └── ...
```

Use this structure with:

```bash
nextflow run main.nf --genomes_base /path/to/genomes_base --genome hg38
```

### Option 2: Direct Reference Path

Alternatively, specify reference files directly:

```bash
# For BWA alignment
nextflow run main.nf --bwa_index_path /path/to/specific/reference.fa

# For methylation analysis
nextflow run main.nf --run_methylation \
                    --methylpy_ref /path/to/methylation/reference_prefix \
                    --methylpy_ref_fasta /path/to/reference.fa
```

This option overrides the standard directory structure and is useful for:
- Custom or non-standard reference genomes
- References located in different directories
- Quick testing with specific reference files

For methylation analysis, the `--methylpy_ref` parameter specifies the prefix path for both forward and reverse methylation references. The pipeline will automatically append "_f" and "_r" to this prefix to locate the forward and reverse reference files, respectively.
<br><br>


## Contributing

This project is under active development. Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
<br>


## Contact

Pin-Chung (Tony) Cheng 
tonycheng521@gmail.com