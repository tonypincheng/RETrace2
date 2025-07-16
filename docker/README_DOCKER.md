# Docker Setup for RETrace2

This guide explains how to run RETrace2 using Docker containers.

## Quick Start

### 1. Build the Custom Docker Image

```bash
# Make the build script executable
chmod +x docker/build_docker.sh

# Build the image (this may take 5-10 minutes)
docker/build_docker.sh
```

Alternatively, build manually:
```bash
docker build -t retrace2/python:latest -f docker/Dockerfile .
```

### 2. Run the Pipeline with Docker

```bash
nextflow run main.nf -profile docker \
  --samplesheet /path/to/samplesheet.csv \
  --genome_base /path/to/genome_base \
  --genome mm39 \
  --output_dir results
```

## How It Works

### Architecture
- **Public Images**: Most processes use pre-built images from biocontainers (e.g., FastQC, BWA, Trim Galore)
- **Custom Image**: Python-based processes use our custom `retrace2/python:latest` image
- **Automatic**: Nextflow handles pulling and running containers automatically

### Container Mapping
| Process | Container |
|---------|-----------|
| FASTQC, MULTIQC | quay.io/biocontainers/fastqc, multiqc |
| TRIM_GALORE | quay.io/biocontainers/trim-galore |
| BWA_MEM | quay.io/biocontainers/bwa |
| COUNT_TARGETS | quay.io/biocontainers/pysam |
| METHYLPY | retrace2/python:latest |
| Custom Python processes | retrace2/python:latest |

### What's in the Custom Image
The `retrace2/python:latest` image includes:
- Python 3.9 + scientific packages (matplotlib, seaborn, pandas, etc.)
- Bioinformatics tools (pysam, biopython, ete3, scikit-bio)
- System tools (bcftools, samtools, tabix)
- **HipSTR** (compiled from source)

## Running Options

### Local Environment
```bash
# Activate your conda environment
conda activate retrace2

# Run without any profile (uses local tools)
nextflow run main.nf --samplesheet samplesheet.csv [other params]
```

### Docker Profile
```bash
# No need to activate conda - containers provide all tools
nextflow run main.nf -profile docker --samplesheet samplesheet.csv [other params]
```


## Volume Mounting

When using Docker, ensure your data is accessible to containers:

```bash
# If your data is outside the project directory, you may need to mount volumes
nextflow run main.nf -profile docker \
  --samplesheet /data/samplesheet.csv \
  --genome_base /references/genome_base \
  -v /data:/data \
  -v /references:/references
```

## Troubleshooting

### Permission Issues
If you encounter permission errors:
```bash
# The docker profile includes user mapping to fix ownership
# This is already configured in nextflow.config
```

### Custom Image Updates
To update the custom image with new dependencies:
1. Modify `docker/Dockerfile`
2. Rebuild: `docker/build_docker.sh`
3. Optionally push to registry for sharing

### Registry Alternative
If you can't build locally, you can pull a pre-built image:
```bash
docker pull retrace2/python:latest
```

## Benefits of Docker Approach

✅ **Reproducible**: Same environment across machines  
✅ **Portable**: Works on any system with Docker  
✅ **Isolated**: No conflicts with system packages  
✅ **Scalable**: Ready for cloud/cluster deployment  
✅ **Maintainable**: Individual containers for different tools  

## Comparison with Local Setup

| Aspect | Local Conda | Docker |
|--------|-------------|--------|
| Setup | Manual conda env | Automatic containers |
| Reproducibility | Environment-dependent | Fully reproducible |
| Resource Usage | Lower overhead | Slight container overhead |
| Debugging | Direct access | Container access needed |
| Sharing | Requires env setup | Just needs Docker | 