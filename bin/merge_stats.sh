#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="/enadisk/gstock/metainvert/results/myriapods_annotation/exogap/intermediate_results/repeats_annotation/repeatmasker"
STATS_BASE="/tempor/merlat/exogap/cache/preprocessing/seqkit_stats"
OUT="/tempor/merlat/exogap/cache/preprocessing/merge_stats/all/all.stats"

# (Re)create output
: > "$OUT"

# Get subdirectories (ids) in BASE_DIR
mapfile -t ids < <(find "$BASE_DIR" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort)

if [[ ${#ids[@]} -eq 0 ]]; then
  echo "Error: no subdirectories (ids) found in: $BASE_DIR" >&2
  exit 1
fi

first_id="${ids[0]}"
first_stats="${STATS_BASE}/${first_id}/${first_id}.stats"

if [[ ! -f "$first_stats" ]]; then
  echo "Error: stats file not found for first id '$first_id': $first_stats" >&2
  exit 1
fi

# 1) Header: take first line of the first id stats file, replace first field by "genome_id"
awk -F'\t' 'BEGIN{OFS="\t"} NR==1 { $1="genome_id"; print; exit }' "$first_stats" >> "$OUT"

# 2) For every id: take second line, replace first field by id, append to OUT
for id in "${ids[@]}"; do
  stats_file="${STATS_BASE}/${id}/${id}.stats"

  if [[ ! -f "$stats_file" ]]; then
    echo "Warning: missing stats file for id '$id': $stats_file (skipping)" >&2
    continue
  fi

  awk -F'\t' -v id="$id" 'BEGIN{OFS="\t"} NR==2 { $1=id; print; found=1; exit } END{ if(!found) exit 3 }' \
    "$stats_file" >> "$OUT" \
    || echo "Warning: could not read 2nd line in '$stats_file' (skipping)" >&2
done

echo "Wrote: $OUT"
