#!/usr/bin/env bash
set -euo pipefail

# run_barrnap.sh
# Args:
#   1: id (sample identifier)
#   2: genome fasta path
#   3: type ("mito" or "nucl" or other)
#   4: threads (optional; default 1)

usage() {
  echo "Usage: $0 <id> <genome_fasta> <type> [threads]" >&2
  exit 1
}

if [[ ${1-} == "" || ${2-} == "" || ${3-} == "" ]]; then
  usage
fi

id="$1"
genome="$2"
type="$3"
threads="${4:-1}"

if [[ ! -s "$genome" ]]; then
  echo "Error: genome file not found or empty: $genome" >&2
  exit 2
fi

# Map type -> kingdom as requested: mito -> mito; nucl -> euk; otherwise keep type
if [[ "$type" == "mito" ]]; then
  kingdom="mito"
elif [[ "$type" == "nucl" ]]; then
  kingdom="euk"
else
  kingdom="$type"
fi

# Run barrnap
/biolo/barrnap/bin/barrnap --quiet --kingdom "$kingdom" --threads "$threads" \
  --outseq "${id}_barrnap_${type}.fa" "$genome" > "${id}_barrnap_${type}.gff.tmp"

# Normalize IDs and expand to gene + rRNA feature rows
grep -v '^#' "${id}_barrnap_${type}.gff.tmp" \
  | sed -E 's/barrnap:[0-9]+(\.[0-9]+)?/barrnap/g' \
  | awk -v OFS='\t' '{
      nb++;
      # gene feature
      print $1,$2,"gene",$4,$5,$6,$7,$8,"ID=barrnap_" nb ";" $9 ";family=Gene,rRNA";
      # child rRNA feature
      print $1,$2,$3,$4,$5,$6,$7,$8,"ID=rnammer_" nb "-rRNA;Parent=barrnap_" nb ";" $9 ";family=Gene,rNA";
    }' \
  > "${id}_barrnap_${type}.gff"

echo "Done. Outputs: ${id}_barrnap_${type}.gff ${id}_barrnap_${type}.fa" >&2
