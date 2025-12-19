#!/bin/bash

set -eou pipefail

# Default values
gff=""
dfam=""
output=""
threads=4

# Function to display usage help
usage() {
    echo "Usage: $0 -g /path/to/gff -d /path/to/dfam -o /path/to/output -t num_threads"
    echo "  -g, --gff         Input GFF file"
    echo "  -d, --dfam        Input DFAM file"
    echo "  -o, --out         Output GFF file"
    echo "  -t, --threads     Number of threads for parallel processing (default: 4)"
    exit 1
}

# Parse command-line options
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -g|--gff) gff="$2"; shift 2 ;;
        -d|--dfam) dfam="$2"; shift 2 ;;
        -o|--out) output="$2"; shift 2 ;;
        -t|--threads) threads="$2"; shift 2 ;;
        *) usage ;;
    esac
done

# Validate required parameters
if [[ -z "$gff" || -z "$dfam" || -z "$output" ]]; then
    echo "Error: Missing required parameters."
    usage
fi

# Check if gff and dfam files exist
if [[ ! -f "$gff" ]]; then
    echo "Error: GFF file does not exist or is not a file."
    exit 1
fi

if [[ ! -f "$dfam" ]]; then
    echo "Error: DFAM file does not exist or is not a file."
    exit 1
fi

# Create directories for chunking and processing
chunk_dir="chunks"
processed_dir="processed_chunks"
mkdir -p "$chunk_dir"
mkdir -p "$processed_dir"

# Step 2: Split the GFF file into chunks
split -l 1000 "$gff" "$chunk_dir/chunk_"  # Split the GFF file into chunks of 10,000 lines each

# Step 3: Process each chunk in parallel
process_chunk() {
    local chunk_file="$1"
    local chunk_output="$2"
    local dfam="$3"

    # Process the chunk
    while IFS= read -r line; do
        if [[ "$line" =~ ^# ]]; then
            echo "$line" >> "$chunk_output"
        else
            # Extract repeat_class_family using shell built-in tools
            class_family=$(echo "$line" | sed -E 's|.*repeat_class_family=([^;]*);.*|\1|' |sed 's/?//g')

            classification=$(awk -F ',' -v class_family="$class_family" '
                    BEGIN { class_family = tolower(class_family) }
                    { if (tolower($8) ~ ("(^|;)" class_family "(;|$)")) print $0 }' "$dfam")

            # Lookup classification in the dfam_map
            if [[ -n "\$classification" ]]; then
                IFS=',' read -r type class order superfamily family subfamily subfamily2 rm_codes <<< "$classification"

                # Replace repeat_class_family with detailed classification (fast inline substitution)
                reformated_line=$(echo "${line}" | sed -E "s|repeat_class_family=${class_family}|repeat_type=${type};repeat_class=${class};repeat_order=${order};repeat_superfamily=${superfamily};repeat_family=${family};repeat_subfamily=${subfamily};repeat_subfamily2=${subfamily2};repeatmasker_code=${class_family}|")
                reformated_line=$(echo "$reformated_line" | sed -E 's/=(;|$)/=NA\1/g')
            else
                # Replace repeat_class_family with detailed classification (fast inline substitution)
                reformated_line=$(echo "${line}" | sed -E "s|repeat_class_family=${class_family}|repeat_type=NA;repeat_class=NA;repeat_order=NA;repeat_superfamily=NA;repeat_family=NA;repeat_subfamily=NA;repeat_subfamily2=NA;repeatmasker_code=${class_family}|")
            fi

            # update the column 3 if it's a tandem repeat
            reformated_line=$(echo "${reformated_line}" |awk -F'\t' '{if ($9 ~ /repeat_type=Tandem Repeat/) $3="tandem_repeat"; print $0}')

            echo "$reformated_line" >> "$chunk_output"
        fi
    done < "$chunk_file"
}

# Step 4: Process chunks in parallel
export -f process_chunk  # Export the function to be used by parallel

# Process each chunk in parallel and save to corresponding processed file
ls "$chunk_dir" | parallel -j "$threads" --bar "process_chunk $chunk_dir/{} $processed_dir/processed_{}" "$dfam"

# Step 5: Combine all processed chunks into the final output
cat "$processed_dir"/processed_* > "$output"

# Clean up chunk and processed directories
rm -rf "$chunk_dir" "$processed_dir"

echo "Processing complete. Output written to $output." 