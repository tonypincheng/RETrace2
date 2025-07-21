# RETrace2 Docker Setup

## Quick Start

1. **Build the custom container:**
   ```bash
   cd docker
   ./build_docker.sh
   ```

2. **Run with Docker profile:**
   ```bash
   nextflow run main.nf -profile docker [other parameters]
   ```

## Important Limitations

### S3/FUSE Filesystem Compatibility Issue

**Docker has known compatibility issues with FUSE-based filesystems (including S3 mounts via mountpoint-s3).**

If your data is stored on S3 mounts, you have two options:

1. **Copy data to local storage** (recommended for Docker):
   ```bash
   # Copy data from S3 mount to local directory
   cp -r /mnt/your/data/path /local/data/path
   
   # Update your samplesheet to use local paths instead of S3 paths
   # Example: change /mnt/data/sample1.fastq.gz to /local/data/path/sample1.fastq.gz
   
   # Then run pipeline with updated samplesheet
   nextflow run main.nf -profile docker --samplesheet /path/to/updated_samplesheet.csv [other parameters]
   ```

2. **Use Singularity instead of Docker** (better S3 compatibility):
   ```bash
   # Note: Singularity profile not yet implemented
   # Contact maintainers if Singularity support is needed
   ```

## Container Details

- **Custom container**: `retrace2/python:latest` 
- **Base image**: python:3.9-slim
- **Key tools**: Python dependencies, HipSTR, Picard, PyQt5
- **Public containers**: BWA, samtools, FastQC, etc. from biocontainers

## Volume Mounting

The Docker profile automatically mounts common directories:
- `/home` - Home directories
- `/opt` - Optional software
- `/data` - Data directories  
- `/shared` - Shared filesystems

For custom mount points, use `--docker-runOptions`:
```bash
nextflow run main.nf -profile docker \
  --docker-runOptions "-v /custom/path:/custom/path" \
  [other parameters]
```

## File Structure

```
docker/
├── Dockerfile              # Custom container definition
├── build_docker.sh         # Build script
└── README_DOCKER.md        # This documentation
```

## Building Custom Container

The custom container includes tools not available in standard biocontainers:

### Included Tools
- **HipSTR**: Microsatellite analysis
- **Picard**: Java-based genomics tools  
- **PyQt5**: GUI support for ete3 tree visualization
- **methylpy**: Methylation analysis
- **Python packages**: pysam, ete3, biopython, pandas, etc.

### Build Process
```bash
cd docker
docker build -t retrace2/python:latest .
```

### Build Options
```bash
# Build with different tag
docker build -t my-retrace2:v1.0 .

# Build with no cache (clean build)
docker build --no-cache -t retrace2/python:latest .
```

## Container Registry

Currently using local containers only. For production deployment:

1. **Tag for registry:**
   ```bash
   docker tag retrace2/python:latest your-registry.com/retrace2/python:latest
   ```

2. **Push to registry:**
   ```bash
   docker push your-registry.com/retrace2/python:latest
   ```

3. **Update nextflow.config** to use registry images

## Docker Hub Registry (Recommended)

**✅ Public Image Available:** `tonypincheng/retrace2-python:latest`

The pipeline is already configured to use a public Docker image. No setup required!

### For AWS Instances
- ✅ **No rebuilding** needed after restarts
- ✅ **Automatic pulling** from Docker Hub
- ✅ **Zero configuration** - just run with `-profile docker`

### Custom Image (Optional)
If you need to modify the container:
```bash
cd docker
./build_docker.sh           # Build the image locally
./setup_registry.sh         # Push to your Docker Hub
```

### Benefits
- ✅ No rebuilding after AWS restarts
- ✅ Automatic image pulling
- ✅ Works across multiple machines
- ✅ Version control for images

## Troubleshooting

### Common Issues

1. **Permission errors**: Check user mapping in docker.runOptions
2. **Mount failures**: Verify paths exist and are readable  
3. **S3 mount errors**: Copy data locally or use Singularity
4. **GUI errors**: PyQt5 should be properly configured with offscreen display
5. **Image not found**: Use Docker Hub registry (see above) to avoid rebuilding

### Debug Commands

```bash
# Test container functionality
docker run --rm retrace2/python:latest python -c "import ete3; print('Success')"

# Test file access
docker run --rm -v /your/data/path:/data retrace2/python:latest ls /data

# Interactive debugging
docker run -it --rm retrace2/python:latest bash
```

### Log Analysis

Check Nextflow logs for Docker-specific errors:
```bash
tail -f .nextflow.log | grep -i docker
``` 