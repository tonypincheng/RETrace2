# RETrace2

> ⚠️ **Under Active Development** ⚠️
> This pipeline is currently under development. Features and documentation may change.

A Nextflow-based bioinformatics pipeline for RETrace2 analysis.

## Overview

RETrace2 is a comprehensive pipeline for processing and analyzing sequencing data to reconstruct single-cell phylogenetic trees using somatic microsatellite mutations. It integrates quality control, alignment, microsatellite genotype calling, and tree reconstruction into a streamlined workflow.

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

## Test Data

The repository includes small test datasets (~60MB total) in the `data/` directory:

### HCT116 Human Cell Line Data
- 6 FASTQ files with ~100,000 reads each
- HCT116 cell line clones from a cell culture tree model
- 12k probe set targeting homopolymers (10-14bp repeat lengths)
- Includes ground truth data for validation
- Sufficient data to test the pipeline's basic functionality and ground truth validation

### MSH2 Mouse Data
- FASTQ files with ~100,000 reads each
- Includes both microsatellite (MS) and methylation (Methyl) libraries
- Designed for testing dual-omic pipeline capabilities
- Sufficient data to test the pipeline's methylation analysis functionality

For details about the test data, see the [data/README.md](data/README.md) file.

## Development Status

### Current Development
- [x] Basic pipeline structure
- [x] Module organization
- [x] Configuration setup
- [ ] HipSTR implementation
- [ ] Tree building implementation
- [ ] Optional modules implementation

### Planned Features
- [ ] Documentation
- [ ] Testing framework
- [x] Example datasets
- [ ] Container support

## Getting Started

### Prerequisites
- nextflow=24.10.5
- fastqc=0.12.1
- multiqc=1.28
- python=3.9
- trim-galore=0.6.10

> **Note:** You can also use Docker or Conda (see below) to handle dependencies automatically.

### Installation
```bash
git clone https://github.com/tonypincheng/retrace2.git
cd retrace2
```

### Usage
```bash
# Run core pipeline
nextflow run main.nf --input_dir data/ --output_dir results/
```
```bash
# Run with bootstrapping
nextflow run main.nf --input_dir data/ --output_dir results/ --run_bootstrap
```
```bash
# Run with evaluation
nextflow run main.nf --input_dir data/ --output_dir results/ --run_evaluation --ground_truth truth.nwk
```


## Execution Environments
RETrace2 supports multiple execution environments through Nextflow profiles.

### Standard (Default)
By default, the pipeline runs with the standard profile, which uses your system's native tools without Docker or Conda.

```bash
# This runs using your locally installed packages
nextflow run main.nf
```

Ensure all required dependencies are installed and in your PATH.

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
│   └── methylpy-ref/      # Methylation references
│       └── ...
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
nextflow run main.nf --bwa_index_path /path/to/specific/reference.fa
```

This option overrides the standard directory structure and is useful for:
- Custom or non-standard reference genomes
- References located in different directories
- Quick testing with specific reference files

## Contributing

This project is under active development. Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

Pin-Chung (Tony) Cheng 
tonycheng521@gmail.com