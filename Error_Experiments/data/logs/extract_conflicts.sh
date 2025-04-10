#!/bin/bash

# Output file
output="conflict_summary.txt"
> "$output"  # Clear previous content

for i in $(seq 1 20); do
    file="miller_rabin_k$i.log"

    if [[ -f "$file" ]]; then
        # Extract the conflict number
        line=$(grep "Processed up to 100000000, conflicts found:" "$file")

        if [[ -n "$line" ]]; then
            # Use sed to extract the number after the last colon
            conflict=$(echo "$line" | sed -E 's/.*conflicts found: ([0-9]+)/\1/')
            echo "k$i: $conflict" >> "$output"
        else
            echo "Line not found for i = $i"
        fi
    else
        echo "File not found: $file"
    fi
done

# Sort the output by k-number (just in case order is off)
sort -n -t 'k' -k2 "$output" -o "$output"

echo "Extraction complete. Results in $output."
