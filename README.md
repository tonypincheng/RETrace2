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

The repository includes a small test dataset (~60MB total) in the `data/` directory.
- 6 FASTQ files with ~100,000 reads each
- HCT116 cell line clones from a cell culture tree model
- 12k probe set targeting homopolymers (10-14bp repeat lengths)
- Sufficient data to test the pipeline's basic functionality

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
- Nextflow
- Python 3.x
- FastQC
- MultiQC
- BWA
- Samtools
- HipSTR


### Installation
```bash
git clone https://github.com/tonypincheng/retrace2.git
cd retrace2
```

### Usage
```bash
# Run core pipeline
nextflow run main.nf --input_dir data/ --output_dir results/

# Run with bootstrapping
nextflow run main.nf --input_dir data/ --output_dir results/ --run_bootstrap

# Run with evaluation
nextflow run main.nf --input_dir data/ --output_dir results/ --run_evaluation --ground_truth truth.nwk
```

## Reference Genome Configuration

RETrace2 provides flexible options for specifying reference genomes:

### Option 1: Standard Directory Structure (Recommended for Research Groups)

The pipeline uses a standardized directory structure for organizing reference genomes:

```
genomes_base/
├── hg38/                  # Genome name
│   ├── bwa-index/         # BWA index files
│   │   └── hg38.fa        # BWA index FASTA
│   └── methylpy/          # Methylation references
│       └── ...
├── mm10/
│   └── ...
```

Use this structure with:

```bash
nextflow run main.nf --genomes_base /path/to/genomes --genome hg38
```

### Option 2: Direct Reference Path (Recommended for Custom Setups)

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