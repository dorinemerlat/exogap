#!/usr/bin/env bash
# Simple: skip first line (header), extract first CSV column
# Usage: bin/get_sample_ids.sh samplesheet.csv

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 samplesheet.csv" >&2
  exit 1
fi

sed '1d' "$1" | cut -d',' -f1 | sed -e 's/^ *"//' -e 's/" *$//' -e 's/^ *//' -e 's/ *$//' -e '/^$/d'
