#!/bin/bash
#SBATCH --job-name=repeatmasker
#SBATCH --output=log_rm/%j.out
#SBATCH --error=log_rm/%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dorine.merlat@etu.unistra.fr
#SBATCH --partition=lab
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=55

set -euo pipefail

# Arguments
ID=$1

GENOME="/enadisk/tempor/merlat/exogaptwo/cache/repeats_annotation/unmask_genome/${ID}/genome_${ID}.fa"
LIB="/enadisk/tempor/merlat/exogaptwo/cache/repeats_annotation/cd_hit/repeats_families/repeats_families.nr.fa"
OUTDIR="/enadisk/tempor/merlat/exogaptwo/cache/repeats_annotation/repeatmasker/${ID}"

echo "Genome   : $GENOME"
echo "ID       : $ID"
echo "Outdir   : $OUTDIR"
echo "Temp dir : $TMPDIR"

# Create output directory
mkdir -p "$OUTDIR"

# Move to node /tmp (or $TMPDIR if defined)
WORKDIR=${TMPDIR:-/tmp}/repeatmasker_${ID}_${SLURM_JOB_ID}
mkdir -p "$WORKDIR"
cd "$WORKDIR"

echo "Running RepeatMasker in $WORKDIR"

/enadisk/tempor/merlat/metainvert/.singularity/quay.io-biocontainers-repeatmasker-4.1.5--pl5321hdfd78af_0.img RepeatMasker \
    -pa 55 \
    -e ncbi \
    -lib "$LIB" \
    -a -gccalc -norna -excln -s -xsmall -gff \
    "$GENOME"

# Unzip output if needed
gunzip $GENOME*.gz

# Move results back to final directory

mv ${GENOME}.* $OUTDIR/

echo "Job finished. Results stored in $OUTDIR"
