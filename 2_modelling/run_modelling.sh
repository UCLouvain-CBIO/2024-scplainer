#!/bin/bash

process_file() {
  file_name=$(basename "$1")
  echo "Processing file: $file_name"
  R CMD BATCH --no-save --no-restore "$1"
}

# Specify the directory you want to process
directory="~/PhD/asca-scp/scripts/2_modelling"
# Expand the tilde using eval
directory=$(eval echo "$directory")

# Use a for loop to iterate through the files in the directory
for file in "$directory"/*; do
  if [ -f "$file" ] && [[ "$file" == *.R ]]; then
    process_file "$file"
  fi
done
