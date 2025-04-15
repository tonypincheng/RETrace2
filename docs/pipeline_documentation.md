# Methylation Analysis Pipeline Documentation

## Overview
This pipeline performs methylation analysis using methylpy, starting from raw FASTQ files.

## Pipeline Steps
1. Quality Control (FastQC)
2. MultiQC Report Generation
3. Read Trimming (trim_galore)
4. Methylation Analysis (methylpy)

## Input Requirements
- Raw FASTQ files (single-end)
- Reference genome files:
  - Forward reference: `${genome}_f`
  - Reverse reference: `${genome}_r`
  - FASTA file: `${genome}.fa`

## Parameters
| Parameter | Description | Default |
|-----------|-------------|---------|
| `genomes_base` | Base directory for reference genomes | `/path/to/reference/genomes` |
| `genome` | Genome version (e.g., mm39) | `mm39` |
| `output_dir` | Output directory | `./results` |

## Dependencies
### Tools
- FastQC (v0.11.9)
- MultiQC (v1.14)
- trim_galore (v0.6.10)
- methylpy (v1.4.5)
- Picard (v3.0.0)

### Python Packages
- methylpy

## Process Details

### 1. Quality Control (FastQC)
- Input: Raw FASTQ files
- Output: FastQC reports (HTML and ZIP)
- Parameters:
  - Quality threshold: 20
  - Phred score: 33
  - Stringency: 3
  - Minimum length: 36

### 2. MultiQC Report
- Input: FastQC reports
- Output: MultiQC HTML report
- Aggregates all QC metrics

### 3. Read Trimming (trim_galore)
- Input: Raw FASTQ files
- Output: Trimmed FASTQ files
- Parameters:
  - Quality threshold: 20
  - Phred score: 33
  - Stringency: 3
  - Minimum length: 36
  - RRBS mode: enabled

### 4. Methylation Analysis (methylpy)
- Input: Trimmed FASTQ files
- Output: Methylation calls and statistics
- Parameters:
  - Minimum quality score: 30
  - Remove clonal reads: False
  - Trim reads: False
  - Number of processors: ${task.cpus}

## Output Structure
```
results/
├── fastqc/              # FastQC reports
├── multiqc_report.html  # MultiQC summary
├── trimmed/            # Trimmed FASTQ files
└── methylation/        # Methylation analysis results
    ├── *.log          # Log files
    └── */*            # Sample-specific results
```

## Troubleshooting
### Common Issues
1. **Memory Issues**
   - Solution: Adjust `task.memory` in nextflow.config

2. **Reference Genome Path**
   - Ensure correct path in `params.genomes_base`

3. **Quality Control Failures**
   - Check input FASTQ quality
   - Adjust trimming parameters if needed

## Development Log
| Date | Change | Description |
|------|--------|-------------|
| YYYY-MM-DD | Initial version | Basic pipeline structure |
| YYYY-MM-DD | Added trim_galore | Replaced trimmomatic with trim_galore |
| YYYY-MM-DD | Updated parameters | Adjusted quality thresholds |

## Future Improvements
- [ ] Add support for paired-end reads
- [ ] Implement batch processing
- [ ] Add more QC metrics
- [ ] Optimize resource usage 