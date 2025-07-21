#!/bin/bash

# Script to set up Docker Hub registry for RETrace2

set -e

IMAGE_NAME="retrace2/python"
LOCAL_TAG="latest"

echo "=== RETrace2 Docker Hub Setup ==="
echo ""

# Check if Docker image exists locally
if ! docker images ${IMAGE_NAME}:${LOCAL_TAG} | grep -q retrace2; then
    echo "Error: Local image ${IMAGE_NAME}:${LOCAL_TAG} not found!"
    echo "Please build it first: ./build_docker.sh"
    exit 1
fi

echo "Setting up Docker Hub registry for RETrace2..."
echo ""
echo "ℹ️  A public Docker image is already available at: tonypincheng/retrace2-python:latest"
echo "   You can use this image directly, or create your own custom version."
echo ""
echo "Choose an option:"
echo "1. Use existing public image (tonypincheng/retrace2-python:latest)"
echo "2. Create your own custom image"
echo ""

read -p "Enter choice (1 or 2): " choice

if [ "$choice" = "1" ]; then
    echo ""
    echo "✅ Using public image: tonypincheng/retrace2-python:latest"
    echo "   Your pipeline is already configured to use this image!"
    echo "   No additional setup needed."
    exit 0
fi

echo ""
echo "Prerequisites for custom image:"
echo "1. Create Docker Hub account at: https://hub.docker.com"
echo "2. Create public repository at: https://hub.docker.com/repository/create"
echo "   - Repository name: retrace2-python"
echo "   - Visibility: Public"
echo ""

read -p "Enter your Docker Hub username: " username

if [ -z "$username" ]; then
    echo "Error: Username cannot be empty"
    exit 1
fi

REGISTRY_IMAGE="${username}/retrace2-python:${LOCAL_TAG}"

echo ""
echo "=== Commands to Run ==="
echo ""
echo "1. Login to Docker Hub:"
echo "   docker login"
echo ""
echo "2. Tag your image:"
echo "   docker tag ${IMAGE_NAME}:${LOCAL_TAG} ${REGISTRY_IMAGE}"
echo ""
echo "3. Push to Docker Hub:"
echo "   docker push ${REGISTRY_IMAGE}"
echo ""
echo "4. Update nextflow.config:"
echo "   Change 'retrace2/python:latest' to '${REGISTRY_IMAGE}'"
echo ""

read -p "Do you want to run these commands now? (y/n): " confirm

if [ "$confirm" = "y" ] || [ "$confirm" = "Y" ]; then
    echo ""
    echo "Running Docker Hub setup..."
    
    echo "Step 1: Logging in to Docker Hub..."
    docker login
    
    echo "Step 2: Tagging image..."
    docker tag ${IMAGE_NAME}:${LOCAL_TAG} ${REGISTRY_IMAGE}
    
    echo "Step 3: Pushing to Docker Hub..."
    docker push ${REGISTRY_IMAGE}
    
    echo ""
    echo "✅ Successfully pushed to Docker Hub!"
    echo ""
    echo "🔧 IMPORTANT: Update your nextflow.config file:"
    echo "   Replace 'retrace2/python:latest' with '${REGISTRY_IMAGE}'"
    echo ""
    echo "🧪 Test the setup:"
    echo "   docker pull ${REGISTRY_IMAGE}"
    echo ""
    echo "🚀 After updating nextflow.config, your pipeline will automatically"
    echo "   pull the image from Docker Hub on any machine - no more rebuilding!"
    
else
    echo ""
    echo "Setup instructions saved. Run the commands above when ready."
fi

echo ""
echo "Registry Image: ${REGISTRY_IMAGE}" 