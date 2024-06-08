#!/bin/bash
# Enter your path to FASTQ files.
DIRECTORY="/path/to/fastq/files"

analyze_fastq() {
    for file in "$DIRECTORY"/*.fastq; do
        if [ -f "$file" ]; then
            echo "Analyzing file: $(basename "$file")"
            analyze_file "$file"
        fi
    done
}

analyze_file() {
    local file="$1"
    local total_reads=0
    local min_length=999999
    local max_length=0

    while read -r line; do
        if [[ "$((++total_reads % 4))" == "2" ]]; then
            local length=${#line}
            if [[ "$length" -lt "$min_length" ]]; then
                min_length="$length"
            fi
            if [[ "$length" -gt "$max_length" ]]; then
                max_length="$length"
            fi
        fi
    done < "$file"

    echo "Total Reads: $total_reads"
    echo "Read Length: $min_length - $max_length"
}

analyze_fastq
