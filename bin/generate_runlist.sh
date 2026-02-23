#!/bin/bash
# Generate VARUS Runlist.txt from SRA run IDs
# Usage:
#   ./generate_runlist.sh SRR123456 SRR789012 ... > Runlist.txt

set -euo pipefail

if [[ $# -eq 0 ]]; then
  echo "Usage: $0 SRR123456 SRR789012 ... > Runlist.txt" >&2
  exit 1
fi

get_metadata() {
  local sra_id="$1"

  # Search
  local search_url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=${sra_id}&usehistory=y&retmode=json"
  local search_result webenv query_key
  search_result=$(curl -s "$search_url")

  webenv=$(echo "$search_result" | jq -r '.esearchresult.webenv')
  query_key=$(echo "$search_result" | jq -r '.esearchresult.querykey')

  [[ -z "$webenv" || "$webenv" == "null" ]] && return 1

  # Fetch metadata
  local fetch_url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=sra&WebEnv=${webenv}&query_key=${query_key}&retmode=json"
  local metadata uid expxml wrapped
  metadata=$(curl -s "$fetch_url")

  uid=$(echo "$metadata" | jq -r '.result.uids[0]')
  expxml=$(echo "$metadata" | jq -r --arg uid "$uid" '.result[$uid].expxml')
  wrapped="<ROOT>$expxml</ROOT>"

  # Parse XML values
  local total_spots total_bases library_layout
  total_spots=$(echo "$wrapped" | xq -r '.ROOT.Summary.Statistics["@total_spots"]')
  total_bases=$(echo "$wrapped" | xq -r '.ROOT.Summary.Statistics["@total_bases"]')
  library_layout=$(echo "$wrapped" | xq -r '.ROOT.Library_descriptor.LIBRARY_LAYOUT | keys[0]')

  # Average read length
  local avg_len
  avg_len=$(echo "scale=2; $total_bases / $total_spots" | bc)

  # Paired flag
  local paired=0
  local paired=0
  [[ "$library_layout" == "PAIRED" ]] && paired=1

  echo -e "${sra_id}\t${total_spots}\t${total_bases}\t${avg_len}\t${paired}\t0"
}

# VARUS header
echo -e "@Run_acc\ttotal_spots\ttotal_bases\tavg_len\tbool:paired\tcolor_space"

# Loop over arguments
for sra_id in "$@"; do
  sra_id=$(echo "$sra_id" | tr -d '[:space:]')
  [[ -z "$sra_id" ]] && continue

  if out="$(get_metadata "$sra_id")"; then
    echo "$out"
  else
    echo -e "${sra_id}\t0\t0\t0\t0\t0"
  fi

  sleep 0.5
done
