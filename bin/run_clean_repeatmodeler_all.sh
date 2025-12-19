#!/usr/bin/env bash
set -euo pipefail

# Simple runner: apply bin/clean_repeat_families.py to all RepeatModeler FASTA files
# Input pattern: cache/repeats_annotation/repeatmodeler/*/*fa
# For each file, determine the species directory name and find the mnemonic by
# reading the first sequence header in the corresponding
# cache/repeats_annotation/clean_repeat_families/<species>/*fa file.
# Outputs from clean_repeat_families.py are moved into the corresponding
# cache/repeats_annotation/clean_repeat_families/<species>/ directory.

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd "$SCRIPT_DIR/.." && pwd)

CLASSIF_CSV="$ROOT_DIR/data/repeats_classification.csv"

if [ ! -f "$CLASSIF_CSV" ]; then
    echo "Classification CSV not found: $CLASSIF_CSV" >&2
    exit 1
fi

# basic checks: python script exists, python3 available, input dir exists
clean_script="$ROOT_DIR/bin/clean_repeat_families.py"
if [ ! -f "$clean_script" ]; then
    echo "ERROR: cleaner script not found: $clean_script" >&2
    exit 1
fi
if ! command -v python3 >/dev/null 2>&1; then
    echo "ERROR: python3 not found in PATH" >&2
    exit 1
fi

repeatmodeler_dir="$ROOT_DIR/cache/repeats_annotation/repeatmodeler"
if [ ! -d "$repeatmodeler_dir" ]; then
    echo "ERROR: repeatmodeler input dir not found: $repeatmodeler_dir" >&2
    exit 1
fi

shopt -s nullglob
files=("$repeatmodeler_dir"/*/*fa)
if [ ${#files[@]} -eq 0 ]; then
    echo "ERROR: no FASTA files found under $repeatmodeler_dir" >&2
    exit 1
fi

shopt -s nullglob

for fasta in "$ROOT_DIR"/cache/repeats_annotation/repeatmodeler/*/*fa; do
    echo "\n==== Processing: $fasta ===="
    species_dir=$(basename "$(dirname "$fasta")")
    clean_dir="$ROOT_DIR/cache/repeats_annotation/clean_repeat_families/$species_dir"

    # determine out_base for produced files
    out_base="$ROOT_DIR/$(basename "$fasta" .fa)"

    # try to find an existing cleaned file to extract mnemonic; prefer any .cleaned_for_mchelper.fa
    mnemonic=''
    if [ -d "$clean_dir" ] && compgen -G "$clean_dir"/*cleaned_for_mchelper.fa > /dev/null 2>&1; then
        first_header=$(head -n 1 "$clean_dir"/*cleaned_for_mchelper.fa | sed -n '1p' || true)
        if [ -n "$first_header" ]; then
            # remove leading '>' and take up to first '_' (e.g. AGAAC_TE1#... -> AGAAC)
            mnemonic=$(echo "$first_header" | sed 's/^>//' | cut -d'_' -f1)
        fi
    else
        # fallback: use first 5 letters of species_dir uppercased
        mnemonic=$(echo "${species_dir:0:5}" | tr '[:lower:]' '[:upper:]')
        echo "Mnemonic not found in $clean_dir; falling back to $mnemonic"
    fi
    echo "Mnemonic: $mnemonic"

    # run the cleaner script
    echo "Running: python3 $ROOT_DIR/bin/clean_repeat_families.py -f $fasta -c $CLASSIF_CSV -m $mnemonic"
    if python3 "$ROOT_DIR/bin/clean_repeat_families.py" -f "$fasta" -c "$CLASSIF_CSV" -m "$mnemonic"; then
        echo "Cleaning succeeded for $fasta"
    else
        echo "Cleaning FAILED for $fasta" >&2
        continue
    fi
    # move outputs into clean_dir (create dir if missing)
    mkdir -p "$clean_dir"
    mchelper_out="${out_base}.cleaned_for_mchelper.fa"
    table_out="${out_base}.families_table.tsv"
    if [ -f "$mchelper_out" ]; then
        mv -f "$mchelper_out" "$clean_dir/"
        echo "Moved: $(basename "$mchelper_out") -> $clean_dir/"
    else
        echo "Warning: expected output not found: $mchelper_out" >&2
    fi
    if [ -f "$table_out" ]; then
        mv -f "$table_out" "$clean_dir/"
        echo "Moved: $(basename "$table_out") -> $clean_dir/"
    else
        echo "Warning: expected table not found: $table_out" >&2
    fi

done
