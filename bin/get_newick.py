#!/usr/bin/env python3

#************************************************#
#            get_taxonomic_lineage.py            #
#            written by Dorine Merlat            #
#                 April 01, 2022                 #
#                                                #
#    Sorting classified and unclassified from    #
#       a library  of repetitive elements.       #
#************************************************#

import sys
import os
import pandas as pd
import requests
import re

def main():
    # Check if there is exactly one argument provided
    if len(sys.argv) != 2:
        print("Error: The info file is missing.")
        sys.exit(1)
    info_file = sys.argv[1] 

    # Check if the provided argument is a valid path
    if not os.path.isfile(info_file):
        print(f"Error: {info_file} is not a valid path.")
        sys.exit(1)

    # Check if the provided file is a valid CSV file
    try:
        df = pd.read_csv(info_file, dtype=str)
    except pd.errors.EmptyDataError:
        print("Error: The info file is not valid.")
        sys.exit(1)
    
    # Check if the required columns are present in the DataFrame
    required_columns = ["id", "taxid"]
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        print(f"Error: The CSV file is missing the required columns: {', '.join(missing_columns)}")
        sys.exit(1)
    
    # Get the phylogeny
    taxid_list = '%2C'.join(df["taxid"].tolist())
    url = "https://lbgi.fr/api/taxonomy/tree/" + taxid_list
    r = requests.get(url, headers={ "Accept" : "application/json"})
    newick = r.json()['data']['tree']

    # remove the taxid id 'root'
    if newick.endswith(")1;"):
        newick = newick[:-3] + ';'  # Remove the ")1"
        newick = newick[1:]   # Remove the first "("
    
    # Group by 'taxid' and aggregate 'id' values into lists
    grouped = df.groupby('taxid')['id'].apply(list).reset_index()
    for index, row in grouped.iterrows(): 
        newick = newick.replace(row['taxid'], ','.join(map(str, row['id']))) # it will join only if they are several ids for the same taxid

    remaining_taxids = re.findall(r'\)(\d+)', newick)
    for taxid in remaining_taxids:
        url = "https://lbgi.fr/api/taxonomy/description/{}".format(str(taxid))
        r = requests.get(url, headers={ "Accept" : "application/json"})
        name = r.json()['data']['name']
        newick = newick.replace(taxid, name)
    
    # Get the root
    # put the line below in a try, if doesn't work then the root is the value (composed of letter, number and '-' character) before the ';' character
    try:
        root = re.findall(r'\)(\w+);', newick)[0]
    except:
        root = re.findall(r'([\w-]+);$', newick)[0]

    print(root)
    with open(root + ".nwk", "w") as file:
        file.write(newick)

if __name__ == "__main__":
    main()
