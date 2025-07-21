#!/bin/bash

# Build script for RETrace2 Docker image

set -e

IMAGE_NAME="retrace2/python"
TAG="latest"

echo "Building Docker image: ${IMAGE_NAME}:${TAG}"

# Build the image
docker build -t ${IMAGE_NAME}:${TAG} -f Dockerfile ..

echo "Build completed successfully!"
echo "Image: ${IMAGE_NAME}:${TAG}"

echo ""
echo "To run with Docker profile, use:"
echo "nextflow run main.nf -profile docker [your parameters]" 