#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <NETMHCpan_SCRIPT_PATH> <NEW_NMHOME_PATH>"
    exit 1
fi

# Assign the provided arguments to variables
NETMHC_SCRIPT_PATH="$1"
NEW_NMHOME_PATH="$2"

# Check if the "netMHCpan" script exists
if [ -e "$NETMHC_SCRIPT_PATH" ]; then
    # Use sed to replace the NMHOME path in the file
    sed -i "s#^setenv\tNMHOME.*#setenv\tNMHOME\t$NEW_NMHOME_PATH#" "$NETMHC_SCRIPT_PATH"
    
    # Create a "tmp" directory in the new NMHOME path if it doesn't exist
    TMP_DIR="$NEW_NMHOME_PATH/tmp"
    if [ ! -d "$TMP_DIR" ]; then
        mkdir -p "$TMP_DIR"
        echo "Created a 'tmp' directory in $NEW_NMHOME_PATH"
    fi
    
    echo "NMHOME path in $NETMHC_SCRIPT_PATH has been updated."
else
    echo "Error: $NETMHC_SCRIPT_PATH does not exist."
fi
