#!/bin/bash

# Original metadata file
original_metadata="/workdir/joo29/metadata_cpp.txt"
# New file to be created, without carriage return characters
cleaned_metadata="/workdir/joo29/metadata.txt"

# Remove carriage return characters (\r) and write to new file
tr -d '\r' < "$original_metadata" > "$cleaned_metadata"

echo "Carriage return characters removed. Cleaned file created at: $cleaned_metadata"
