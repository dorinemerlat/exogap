#!/usr/bin/env bash
# Script to search SRA database for RNA-seq data
# Usage: search_sra.sh -t <taxid> -n <name> [-k <api_key>]

#SBATCH --job-name=dwnld
#SBATCH --output=/home/merlat/log/sra/%j.out
#SBATCH --error=/home/merlat/log/sra/%j.err
#SBATCH --mail-type ALL
#SBATCH --mail-user dorine.merlat@etu.unistra.fr
#SBATCH --partition lab
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks-per-node=1         # Number of tasks (Nextflow runs as a single task)

TAXID=""
NAME=""
API_KEY=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -t)
            TAXID="$2"
            shift 2
            ;;
        -n)
            NAME="$2"
            shift 2
            ;;
        -k|--key)
            API_KEY="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

echo "Searching SRA for taxid: ${TAXID}, name: ${NAME}"

# Check mandatory arguments
if [[ -z "${TAXID}" ]]; then
    echo "Error: -t (taxid) is required" >&2
    exit 1
fi

if [[ -z "${NAME}" ]]; then
    echo "Error: -n (name) is required" >&2
    exit 1
fi

# Export API key only if provided
if [[ -n "${API_KEY}" ]]; then
    export NCBI_API_KEY="${API_KEY}"
fi

{ # try a first time
    esearch -db sra -query "((((txid${TAXID}[Organism:exp]) AND \"paired\"[Layout]) AND \"illumina\"[Platform]) AND \"transcriptomic\"[Source]) AND \"filetype fastq\"[Properties]" \
        > esearch.out

} || { # try a second time
    sleep $(shuf -i 5-30 -n 1)
    esearch -db sra -query "((((txid${TAXID}[Organism:exp]) AND \"paired\"[Layout]) AND \"illumina\"[Platform]) AND \"rna data\"[Filter]) AND \"filetype fastq\"[Properties]" \
        > esearch.out
}

# remove header
count=$(grep "<Count>" esearch.out | cut -d '>' -f 2 | cut -d '<' -f 1)
echo "${TAXID},${NAME},${count}" > ${TAXID}_sra.count

if [[ ${count} != 0 ]] ; then
    timer=$((3 * (1 + RANDOM % count)))
    sleep ${timer}
    efetch -format runinfo < esearch.out | cut -f 1 -d ',' | sed 1d > efetch.out

    while read -r rna; do
        echo "${rna}" > ${rna}.list
    done < efetch.out
else
    touch 0000.list
fi
