#!/bin/bash
# Loop through all .fastq.gz files in the current directory
for file in *.fastq.gz; do
    # Extract the beginning of the filename until the first underscore
    folder_name=$(echo "$file" | cut -d'_' -f1)
    
    # Create a folder with the extracted name if it doesn't exist
    if [ ! -d "$folder_name" ]; then
        mkdir "$folder_name"
        echo "Created folder: $folder_name"
    fi

    echo "Ungzipping $file"
    gunzip $file

    mv $(basename "$file" .gz) $folder_name
done

echo "Completed successfully."

