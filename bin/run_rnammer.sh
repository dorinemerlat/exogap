#!/usr/bin/env bash
set -euo pipefail

# run_rnammer.sh
# Args:
#   1: id (sample identifier)
#   2: genome fasta path

usage() {
  echo "Usage: $0 <id> <genome_fasta>" >&2
  exit 1
}

if [[ ${1-} == "" || ${2-} == "" ]]; then
  usage
fi

id="$1"
genome="$2"

if [[ ! -s "$genome" ]]; then
  echo "Error: genome file not found or empty: $genome" >&2
  exit 2
fi

# Run rnammer (Eukaryote; 5S/SSU/LSU)
rnammer -S euk -m tsu,ssu,lsu -multi \
  -gff "${id}_rnammer.gff.tmp" \
  -f   "${id}_rnammer.fa.tmp" \
  -h   "${id}_rnammer.report" < "$genome"

# Fix GFF: rename 8s_rRNA -> 5s_rRNA, then emit gene + rRNA features
grep -v '^#' "${id}_rnammer.gff.tmp" \
  | awk -F'\t' -v OFS='\t' '$9 == "8s_rRNA" { $9 = "5s_rRNA" } { print }' \
  | awk -v OFS='\t' '{
      nb++;
      print $1,$2,"gene",$4,$5,$6,$7,$8,"ID=rnammer_" nb ";Name=" $9 ";family=Gene,rRNA";
      print $1,$2,$3,$4,$5,$6,$7,$8,"ID=rnammer_" nb "-rRNA;Name=" $9 ";Parent=rnammer_" nb ";family=Gene,rRNA";
    }' \
  > "${id}_rnammer.gff"

# Fix FASTA header molecule tag
sed "s|molecule=8s_rRNA|molecule=5s_rRNA|g" "${id}_rnammer.fa.tmp" > "${id}_rnammer.fa"

# Normalize case in all produced files
sed -i 's/s_rRNA/S_rRNA/g' ${id}_rnammer.*

echo "Done. Outputs: ${id}_rnammer.gff ${id}_rnammer.fa ${id}_rnammer.report" >&2
