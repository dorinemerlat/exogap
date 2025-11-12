#!/usr/bin/env python3
"""
Clean header of a repeat families file.
Usage: clean_repeat_families.py -f <fasta> -c <classification.csv> -s <source> -m <mnemonic>
"""

import sys
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import csv

def load_classification(path):
    classification = {}

    # read CSV without index column
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Get the codes, split by ';' and remove extra spaces
            codes = [c.strip() for c in row["Code"].split(";")]

            # Prepare the value: all columns except "Code"
            value = {k: v for k, v in row.items() if k != "Code"}

            # Add one entry per code
            for code in codes:
                classification[code.lower()] = value

    return classification


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('-f', '--fasta', required=True, help='Input FASTA file')
    p.add_argument('-c', '--classification', required=True, help='Classification CSV file')
    p.add_argument('-m', '--mnemonic', required=True, help='Mnemonic of species')
    p.add_argument('-s', '--source', required=True, choices=['RM2', 'HC', 'HLC'], help='Source tag to add to headers')
    args = p.parse_args()

    classification = load_classification(args.classification)

    # determine output paths
    base = args.fasta.rsplit('.', 1)[0]
    main_out = f"{base}.cleaned.fa"
    mchelper_out = f"{base}.cleaned_for_mchelper.fa"

    # iterate sequences and build cleaned records
    out_records = []
    out_records_mchelper = []

    counter = 1  # start numbering sequences at 1

    for rec in SeqIO.parse(args.fasta, 'fasta'):
        seq_id = rec.id
        header = rec.description

        # try to extract a classification token from the seq_id or header
        parts = seq_id.split('#', 1)
        repeat_id = parts[0].strip()

        if len(parts) > 1:
            repeat_classification = parts[1].strip()
        else:
            repeat_classification = ''

        # If the classification is already 'Unknown', keep it
        if repeat_classification.lower() == 'unknown':
            clean_classification = 'Unknown'
            clean_classification_short = 'Unknown'
            matched = None
        else:
            matched = classification.get(repeat_classification.lower()) if repeat_classification else None
            if matched is not None:
                cls = (matched.get('Class') or '').strip()
                order = (matched.get('Order') or '').strip()
                superfamily = (matched.get('SuperFamily') or '').strip()
                family = (matched.get('Family') or '').strip()
                subfamily = (matched.get('SubFamily') or '').strip()
                subfamily2 = (matched.get('Subfamily2') or '').strip()
                clean_classification = '/'.join(filter(None, [cls, order, superfamily, family, subfamily, subfamily2]))
                clean_classification_short = '/'.join(filter(None, [cls, order, superfamily]))
            else:
                clean_classification = 'Unknown'
                clean_classification_short = 'Unknown'

        # New sequence ID format: MNEMONIC_TE1, MNEMONIC_TE2, etc.
        new_id = f"{args.mnemonic}_TE{counter}"
        counter += 1

        # prepare SeqRecord for general output
        repeat_header = f"{args.source}_{new_id}#{clean_classification}"
        main_record = SeqRecord(rec.seq, id=repeat_header, description='')
        out_records.append(main_record)

        # prepare SeqRecord for mchelper (short)
        if matched and 'CMC' in clean_classification:
            clean_classification_short = 'ClassII/TIR/CMC'
        if matched and 'Mutator' in clean_classification:
            clean_classification_short = 'ClassII/TIR/MULE'
        repeat_header = f"{args.source}_{new_id}#{clean_classification_short}"
        mchelper_record = SeqRecord(rec.seq, id=repeat_header, description='')
        out_records_mchelper.append(mchelper_record)

    # write outputs in append mode so we don't overwrite existing files
    with open(main_out, 'a') as out:
        SeqIO.write(out_records, out, 'fasta')
    with open(mchelper_out, 'a') as out:
        SeqIO.write(out_records_mchelper, out, 'fasta')


if __name__ == '__main__':
    main()
