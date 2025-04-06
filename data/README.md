# Test Datasets for RETrace2

This directory contains test data for the RETrace2 pipeline. These datasets are small subsets of real data, downsampled to allow for quick testing and demonstration of the pipeline functionality.

## Dataset Description

The test data consists of FASTQ files from HCT116 single cell samples with microsatellite targets. Each file contains approximately 100,000 reads (9-10MB per file).

Files are named according to the following convention:
`MS_1nt12k10-14bp_HCT116_SC_[clone-id]_[plate-info]_[reads-count].fastq.gz`

### Included Files:

- **Cell Clone 2-1-C4_3-7-G7 pair**:
  - MS_1nt12k10-14bp_HCT116_SC_2-1-C4_3-7-G7_20230110_Plate2-2G_100000reads.fastq.gz
  - MS_1nt12k10-14bp_HCT116_SC_2-1-C4_3-7-G7_20230110_Plate2-3G_100000reads.fastq.gz

- **Cell Clone 2-1-F7_3-3-B9 pair**:
  - MS_1nt12k10-14bp_HCT116_SC_2-1-F7_3-3-B9_20230110_Plate2-4B_100000reads.fastq.gz
  - MS_1nt12k10-14bp_HCT116_SC_2-1-F7_3-3-B9_20230110_Plate2-4E_100000reads.fastq.gz

- **Cell Clone 2-1-H9_3-4-D8 pair**:
  - MS_1nt12k10-14bp_HCT116_SC_2-1-H9_3-4-D8_20230110_Plate2-7E_100000reads.fastq.gz
  - MS_1nt12k10-14bp_HCT116_SC_2-1-H9_3-4-D8_20230110_Plate2-9E_100000reads.fastq.gz

## Data Content

These files contain microsatellite data with the following characteristics:
- Part of the HCT116 cell culture tree model
- 12k probe set (a smaller subset of the full probe set for testing)
- Targeting homopolymers between 10-14bp repeat lengths
- Single cell samples

## Usage

These test datasets are designed for quick validation of the RETrace2 pipeline. They can be used with the default pipeline configuration:

```bash
nextflow run main.nf --input_dir data/ --output_dir results/
```

## Data Size

Each file is approximately 10MB in size, with the total test dataset being around 60MB. This is intentionally kept small to:
1. Allow easy repository cloning
2. Enable quick pipeline testing
3. Demonstrate basic functionality without requiring extensive compute resources

## Source

These datasets are subsampled from full experimental data. The original data contains multiple million reads per sample and would require more computational resources to process.

## Limitations

Due to the reduced size of these test datasets:
- Coverage of microsatellite loci is lower than in production data
- Some analyses may show higher noise levels than with full datasets
- Results should be considered demonstrative rather than biologically meaningful

For production use, it's recommended to use the full datasets or your own experimental data. 