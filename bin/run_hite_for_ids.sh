#!/usr/bin/env bash
set -euo pipefail

# Usage: run_hite_for_ids.sh --ids id1,id2,... [--singularity IMAGE] [--genome-dir DIR]
# Example:
#   ./bin/run_hite_for_ids.sh --ids cylindroiulus-punctatus,other-species
# The script will:
# - create /data/merlat/sandbox/hite/<id> for each id
# - expect a genome file named genome_<id>.fa in the current directory or in --genome-dir
# - bind the current directory and the work_dir into the singularity container

IMAGE="/enadisk/tempor/merlat/metainvert/.singularity/dorinemerlat-exogap_hite-1.0.img"
GENOME_DIR="."
IDS=""
THREADS=30
OUT_BASE="out"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --ids)
      IDS="$2"
      shift 2
      ;;
    --singularity)
      IMAGE="$2"
      shift 2
      ;;
    --genome-dir)
      GENOME_DIR="$2"
      shift 2
      ;;
    --threads)
      THREADS="$2"
      shift 2
      ;;
    --out-base)
      OUT_BASE="$2"
      shift 2
      ;;
    -h|--help)
      sed -n '1,120p' "$0"
      exit 0
      ;;
    *)
      echo "Unknown arg: $1" >&2
      exit 1
      ;;
  esac
done

if [[ -z "$IDS" ]]; then
  echo "Please pass --ids id1,id2,..." >&2
  exit 1
fi

IFS=',' read -r -a ID_ARR <<< "$IDS"

PWD_DIR="$(pwd)"

for id in "${ID_ARR[@]}"; do
  id=$(echo "$id" | xargs)
  if [[ -z "$id" ]]; then
    continue
  fi

  work_dir="/data/merlat/sandbox/hite/${id}"
  mkdir -p "$work_dir"
  echo "Created work dir: $work_dir"

  # genome file name - try different common patterns
  candidates=("genome_${id}.fa" "genome_${id}.fasta" "${id}.fa" "${id}.fasta")
  genome_file=""
  for c in "${candidates[@]}"; do
    if [[ -f "$GENOME_DIR/$c" ]]; then
      genome_file="$GENOME_DIR/$c"
      break
    fi
  done

  if [[ -z "$genome_file" ]]; then
    echo "Genome file for id '$id' not found in $GENOME_DIR. Tried: ${candidates[*]}" >&2
    continue
  fi

  echo "Using genome file: $genome_file"

  # prepare output dir under current directory to keep results per id
  out_dir="$PWD_DIR/${OUT_BASE}_${id}"
  mkdir -p "$out_dir"

  # Run singularity: bind current dir and work dir. Also bind genome dir to /data/genomes inside container
  # The HiTE command expects python /HiTE/main.py ... --work_dir <work_dir>
  echo "Running HiTE for $id..."

  singularity run \
    --bind "$PWD_DIR:$PWD_DIR" \
    --bind "$work_dir:$work_dir" \
    --bind "$GENOME_DIR:/data/genomes" \
    "$IMAGE" \
    python /HiTE/main.py \
      --genome "$genome_file" \
      --thread "$THREADS" \
      --out_dir "$out_dir" \
      --domain 1 \
      --plant 0 \
      --recover 1 \
      --search_struct 1 \
      --work_dir "$work_dir"

  echo "Finished $id. Output in $out_dir, work dir $work_dir"
done
