# RETrace2

> ⚠️ **Under Active Development** ⚠️
> This pipeline is currently under development. Features and documentation may change.

A Nextflow-based bioinformatics pipeline for RETrace2 analysis.

## Overview

RETrace2 is a comprehensive pipeline for processing and analyzing sequencing data to reconstruct phylogenetic trees based on microsatellite mutations. The pipeline includes quality control, alignment, microsatellite calling, and phylogenetic tree reconstruction steps.

## Features

- Core pipeline:
  - FastQC for quality control
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
│   ├── core/             # Core phylogenetic tree pipeline
│   ├── bootstrap/        # Tree bootstrapping analysis
│   ├── methylation/      # Methylation analysis
│   └── evaluation/       # Tree accuracy evaluation
├── scripts/              # Supporting scripts
├── main.nf               # Main pipeline
└── nextflow.config      # Pipeline configuration
```

## Development Status

### Current Development
- [x] Basic pipeline structure
- [x] Core module organization
- [x] Configuration setup
- [ ] HipSTR implementation
- [ ] Tree building implementation
- [ ] Optional modules implementation

### Planned Features
- [ ] Documentation
- [ ] Testing framework
- [ ] Example datasets
- [ ] Container support

## Getting Started

### Prerequisites
- Nextflow
- BWA
- HipSTR
- Python (for tree building)

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
tonycheng521@gmal.com