#!/usr/bin/env python3

import argparse
import gffpandas.gffpandas as gffpd
import pandas as pd
import os

def parse_attributes(s):
    if pd.isna(s) or s == "":
        return {}
    attributes = {}
    for item in str(s).split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            attributes[k] = v
    return attributes

def parse_gff(path):
    """Parse a GFF file into a DataFrame and expand attributes into columns."""
    gff_obj = gffpd.read_gff3(path)
    name = os.path.basename(path).split(".")[0]

    gff = gff_obj.df.copy()
    gff.insert(0, "genome_id", name)

    # filter to keep only lines with type "Transposable_elements"
    gff = gff[gff["type"] == "Transposable_elements"].copy()

    # Parse attributes into dict
    attr_dicts = gff["attributes"].apply(parse_attributes)

    # Expand dicts into columns
    attr_df = pd.json_normalize(attr_dicts)

    # Merge expanded attributes with the main dataframe
    gff = gff.reset_index(drop=True)

    attr_dicts = gff["attributes"].apply(parse_attributes)
    attr_df = pd.json_normalize(attr_dicts).reset_index(drop=True)

    gff = pd.concat([gff.drop(columns=["attributes"]), attr_df], axis=1)

    return gff

def summarize_gff(df):
    """Summarize the GFF DataFrame by counting entries per genome_id and type."""
    summary = df.groupby(["genome_id", "total_length", "repeat_type", "repeat_class", "repeat_order", "repeat_superfamily", "repeat_family"]).size().reset_index(name="count")

    # if total_length is a string
    df["total_length"] = pd.to_numeric(df["total_length"], errors="coerce")

    summary = (
        df.groupby(
            ["genome_id", "repeat_type", "repeat_class", "repeat_order",
             "repeat_superfamily", "repeat_family"],
            dropna=False
        )
        .agg(
            count=("genome_id", "size"),
            total_length_sum=("total_length", "sum")
        )
        .reset_index()
    )

    return summary

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Parse a GFF file and extract attributes.")
    parser.add_argument('-i', '--input', required=True, help="Input GFF file")
    parser.add_argument('-o', '--output', required=True, help="Output TSV file")

    # Parse the command-line arguments
    args = parser.parse_args()

    # Parse the GFF file
    df = parse_gff(args.input)
    summary = summarize_gff(df)

    # Save the result to a TSV file
    summary.to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    main()
