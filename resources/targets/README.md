# Target BED Files for RETrace2

This directory contains target BED files that define microsatellite loci for RETrace2 analysis. These files specify the genomic regions where microsatellites will be called using HipSTR.

## File Format

RETrace2 uses a **HipSTR-compatible BED format** with the following specifications:

### Column Structure
| Column | Name | Description | Example |
|--------|------|-------------|---------|
| 1 | Chromosome | Chromosome identifier | `chr6` |
| 2 | Start | Start position (**1-based**) | `8748168` |
| 3 | End | End position (**1-based, inclusive**) | `8748178` |
| 4 | Motif Length | Length of repeat unit in base pairs | `1` |
| 5 | Repeat Count | Number of repeat copies in reference genome | `10` |
| 6 | Target ID | Unique identifier for the target region | `chr6:8748168-8748178_10xT` |

### Example Lines
```
chr6	8748168	8748178	1	10	chr6:8748168-8748178_10xT
chr14	63837812	63837822	1	10	chr14:63837812-63837822_10xG
chr15	4945956	4945967	1	11	chr15:4945956-4945967_11xG
```

## Important Notes

### ⚠️ **1-Based Coordinates**
These BED files use **1-based coordinates** (start and end positions), which differs from the standard BED format that uses 0-based coordinates. This is required for HipSTR compatibility.

### Homopolymer Targets
RETrace2 focuses on **homopolymer microsatellites** (single nucleotide repeats):
- **Motif Length**: Always `1` for homopolymers (A, T, G, or C repeats)
- **Repeat Count**: Number of consecutive identical nucleotides in the reference genome
- **Target ID**: Descriptive format showing location and repeat type (e.g., `_10xT` = 10 consecutive T nucleotides)

### Target ID Convention
The target ID follows the pattern: `{chr}:{start}-{end}_{count}x{nucleotide}`
- `chr6:8748168-8748178_10xT` = 10 consecutive T nucleotides on chr6 from position 8748168-8748178

## Available Target Sets

### Mouse (mm39)
- **RETrace2.mm39.1nt10-30bp.92460targets169818probes.bed** (5.1MB)
  - 92,460 unique targets
  - 169,818 total probes
  - Homopolymer lengths: 10-30 base pairs
  - Comprehensive genome-wide coverage

### Human (hg38)
- **RETrace2.hg38.1nt10-25bp.68807targets119510probes.bed** (3.8MB)
  - 68,807 unique targets  
  - 119,510 total probes
  - Homopolymer lengths: 10-25 base pairs
  - Comprehensive genome-wide coverage

- **RETrace2.hg38.1nt10-14bp.12000probes.bed** (387KB)
  - ~12,000 probes
  - Homopolymer lengths: 10-14 base pairs  
  - Subset for testing and validation

## Usage in RETrace2

Specify the target BED file using the `--target_bed` parameter:

```bash
# Using provided targets
nextflow run main.nf --target_bed resources/targets/mm39/RETrace2.mm39.1nt10-30bp.92460targets169818probes.bed

# Using custom targets
nextflow run main.nf --target_bed /path/to/custom_targets.bed
```