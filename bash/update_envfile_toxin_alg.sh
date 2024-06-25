#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <ENVFILE_PATH> <BLAST_PATH> <BLAST_DB_PATH> <MERCI_PATH> <MERCI_MOTIF_PATH>"
    exit 1
fi

# Assign the provided arguments to variables
ENVFILE_PATH="$1"
NEW_BLAST_PATH="$2"
NEW_BLAST_DB_PATH="$3"
NEW_MERCI_PATH="$4"
NEW_MERCI_MOTIF_PATH="$5"

# Check if the "envfile" exists
if [ -e "$ENVFILE_PATH" ]; then
    # Use sed to replace the paths in the file, preserving the "$" symbol
    sed -i "s#^BLAST:.*#BLAST: \$${NEW_BLAST_PATH}#" "$ENVFILE_PATH"
    sed -i "s#^BLAST database:.*#BLAST database: \$${NEW_BLAST_DB_PATH}#" "$ENVFILE_PATH"
    sed -i "s#^MERCI:.*#MERCI: \$${NEW_MERCI_PATH}#" "$ENVFILE_PATH"
    sed -i "s#^MERCI motif file:.*#MERCI motif file: \$${NEW_MERCI_MOTIF_PATH}#" "$ENVFILE_PATH"
    
    echo "Paths in $ENVFILE_PATH have been updated."
else
    echo "Error: $ENVFILE_PATH does not exist."
fi
