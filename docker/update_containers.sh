#!/bin/bash

# Script to update all container references to use Docker Hub image

set -e

echo "=== RETrace2 Container Update Script ==="
echo ""

# Get the Docker Hub username and image
read -p "Enter your Docker Hub username: " username

if [ -z "$username" ]; then
    echo "Error: Username cannot be empty"
    exit 1
fi

OLD_IMAGE="retrace2/python:latest"
NEW_IMAGE="${username}/retrace2-python:latest"

echo ""
echo "Updating container references:"
echo "  FROM: ${OLD_IMAGE}"
echo "  TO:   ${NEW_IMAGE}"
echo ""

# List of files that contain the container references
FILES=(
    "modules/hipstr/hipstr.nf"
    "modules/phylo/phylo.nf"
    "modules/methylation/methylation.nf"
    "modules/infer_celltype/infer_celltype.nf"
    "modules/stats/stats.nf"
    "modules/evaluation/evaluation.nf"
)

# Count total occurrences
total_count=0
for file in "${FILES[@]}"; do
    if [ -f "$file" ]; then
        count=$(grep -c "container \"${OLD_IMAGE}\"" "$file" 2>/dev/null || echo 0)
        total_count=$((total_count + count))
        if [ $count -gt 0 ]; then
            echo "  📁 $file: $count occurrence(s)"
        fi
    fi
done

echo ""
echo "Total occurrences to update: $total_count"

if [ $total_count -eq 0 ]; then
    echo "No occurrences found. Nothing to update."
    exit 0
fi

echo ""
read -p "Proceed with the update? (y/n): " confirm

if [ "$confirm" != "y" ] && [ "$confirm" != "Y" ]; then
    echo "Update cancelled."
    exit 0
fi

echo ""
echo "Updating files..."

# Update each file
updated_count=0
for file in "${FILES[@]}"; do
    if [ -f "$file" ]; then
        # Create backup
        cp "$file" "${file}.backup"
        
        # Perform replacement
        if sed -i "s|container \"${OLD_IMAGE}\"|container \"${NEW_IMAGE}\"|g" "$file"; then
            count=$(grep -c "container \"${NEW_IMAGE}\"" "$file" 2>/dev/null || echo 0)
            if [ $count -gt 0 ]; then
                echo "  ✅ Updated $file ($count occurrences)"
                updated_count=$((updated_count + count))
            fi
        else
            echo "  ❌ Failed to update $file"
            # Restore backup on failure
            mv "${file}.backup" "$file"
        fi
    else
        echo "  ⚠️  File not found: $file"
    fi
done

echo ""
if [ $updated_count -gt 0 ]; then
    echo "✅ Successfully updated $updated_count container references!"
    echo ""
    echo "🔧 Next steps:"
    echo "1. Test your pipeline: nextflow run main.nf -profile docker [your parameters]"
    echo "2. Docker will automatically pull ${NEW_IMAGE} when needed"
    echo ""
    echo "📋 Backup files created (*.backup) - remove them if everything works:"
    echo "   find modules/ -name '*.backup' -delete"
else
    echo "❌ No files were updated"
fi

echo ""
echo "🚀 Your pipeline is now configured to use Docker Hub!"
echo "   No more rebuilding needed after AWS restarts!" 