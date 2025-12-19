#!/usr/bin/env python3

import argparse
import pandas as pd
import pybedtools

def process_tsv(input_file, genome_size):
    """Process the TSV file to summarize repeats"""
    # Define the grouping columns
    grouping_columns = [
        "genome", "repeat_type", "repeat_class", "repeat_order", "repeat_superfamily", 
        "repeat_family", "repeat_subfamily", "repeat_subfamily2"
    ]

    # Columns to keep include grouping columns + 'total_length'
    columns_to_keep = grouping_columns + ["total_length"] 

    # Read the TSV file into a DataFrame
    df = pd.read_csv(input_file, sep="\t", low_memory=False, dtype=str)


    # Ensure required columns exist
    missing_columns = [col for col in columns_to_keep if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing columns in input file: {', '.join(missing_columns)}")

    # Select only the necessary columns
    df = df[columns_to_keep]

    # Fill NaN values with a placeholder ('NA')
    df[grouping_columns] = df[grouping_columns].fillna("NA")

    # Convert total_lenght to int (
    df["total_length"] = df["total_length"].astype(int)

    # Group by the specified columns and aggregate
    result = df.groupby(grouping_columns).aggregate(
        coverage=("total_length", "sum"),
        load=("total_length", "count")
    ).reset_index() 

    # Calculate coverage percentage
    result["coverage"] = (result["coverage"] * 100) / genome_size

    # Rename columns to remove 'repeat_' prefix
    result.rename(columns={col: col.replace("repeat_", "") for col in result.columns}, inplace=True)

    return result

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Summarize repeats TSV")
    parser.add_argument('-i', '--input', required=True, help="Input TSV file")
    parser.add_argument('-o', '--output', required=True, help="Output TSV file")
    parser.add_argument('-s', '--size', required=True, type=int, help="Genome size")

    # Parse the command-line arguments
    args = parser.parse_args()

    # Process the TSV file
    result = process_tsv(args.input, args.size)

    # Save the result to a TSV file
    result.to_csv(args.output, sep="\t", index=False)

    print(f"Summarized file saved to {args.output}")

if __name__ == "__main__":
    main()
