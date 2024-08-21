#!/bin/bash

# Check if the correct number of arguments is passed
if [ $# -ne 1 ]; then
    echo "Usage: $0 {on|off}"
    exit 1
fi

# Assign the parameter to a variable
action=$1

# Define the file to operate on
file="_config.yml"

# Perform the action based on the parameter
if [ "$action" == "on" ]; then
    sed -i '' '8s/^#//' "$file"    # Uncomment line 8 (remove # at the start)
    sed -i '' '9s/^[^#]/#&/' "$file"  # Comment line 9 (add # at the start if not present)
elif [ "$action" == "off" ]; then
    sed -i '' '8s/^[^#]/#&/' "$file"  # Comment line 8 (add # at the start if not present)
    sed -i '' '9s/^#//' "$file"    # Uncomment line 9 (remove # at the start)
else
    echo "Invalid parameter. Use 'on' or 'off'."
    exit 1
fi
