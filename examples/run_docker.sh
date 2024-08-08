#!/bin/bash

# List all subfolders in the current directory
for dir in */; do
    echo "Checking folder: $dir"
    cd "$dir"
    if [ -f run_docker.sh ]; then
        echo "Running run_docker.sh in $dir"
        ./run_docker.sh
    else
        echo "No run script in $dir"
    fi
    cd ..
done
