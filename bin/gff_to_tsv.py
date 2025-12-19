#!/usr/bin/env python3

import argparse
import pandas as pd

def parse_attributes(attribute_string):
    """Parses the attributes column in a GFF file into a dictionary."""
    attributes = {}
    for attr in attribute_string.split(";"):
        key_value = attr.split("=")
        if len(key_value) == 2:
            key, value = key_value
            attributes[key] = value
    return attributes

def parse_gff_line(line):
    """Parse a GFF line into a dictionary."""
    fields = line.strip().split("\t")

    first_fields = {
        "seqid": fields[0],
        "source": fields[1],
        "type": fields[2],
        "start": fields[3],
        "end": fields[4],
        "score": fields[5],
        "strand": fields[6],
        "phase": fields[7],
    }

    attribute_fields = parse_attributes(fields[8])

    return {**first_fields, **attribute_fields}

def parse_gff(gff_file):
    """Parses the GFF file and returns a DataFrame."""
    entries = []

    with open(gff_file, 'r') as file:
        for line in file:
            if not line.startswith("#"):  # Ignore header lines
                entry = parse_gff_line(line)
                entries.append(entry)

    # Get all column names from the first entry to ensure all attributes are included
    all_columns = set()
    for entry in entries:
        all_columns.update(entry.keys())

    # Classic GFF columns (always present)
    classic_columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase']
    
    # Get the attribute columns (everything except the classic columns)
    attribute_columns = sorted(all_columns - set(classic_columns))
    
    # Final column order: first the classic columns, then the sorted attribute columns
    final_columns = classic_columns + attribute_columns

    # Convert the list of dictionaries into a DataFrame with the correct column order
    df = pd.DataFrame(entries, columns=final_columns)

    return df

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Parse a GFF file and extract attributes.")
    parser.add_argument('-i', '--input', required=True, help="Input GFF file")
    parser.add_argument('-o', '--output', required=True, help="Output TSV file")

    # Parse the command-line arguments
    args = parser.parse_args()

    # Parse the GFF file
    df = parse_gff(args.input)
    
    # Save the result to a TSV file
    df.to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    main()
