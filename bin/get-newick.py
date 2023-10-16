#!/usr/bin/env python3

#************************************************#
#            get_taxonomic_lineage.py            #
#            written by Dorine Merlat            #
#                 April 01, 2022                 #
#                                                #
#    Sorting classified and unclassified from    #
#       a library  of repetitive elements.       #
#************************************************#

import requests
import argparse
import json
import re
import newick

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--ids", help="list of genome's taxids", type=str, required=True)
    args = parser.parse_args()

    IDs = dict()
    for i in re.findall("(\d+):([\w\s]+)", args.ids):
        IDs[i[1]] = i[0]

    url = "https://lbgi.fr/api/taxonomy/tree/{}".format('%2C'.join(IDs.values()))
    r = requests.get(url, headers={ "Accept" : "application/json"})
    newick_file = r.json()['data']['tree']

    # add all ID for genomes with same taxID (present only once in the newick)
    duplicates_ids = list( {x for x in IDs.values() if list(IDs.values()).count(x) > 1} )
    duplicate_names = dict()
    for duplicate_id in duplicates_ids:
        for name, taxid in IDs.items():
            if taxid == duplicate_id:
                if taxid not in duplicate_names:
                    duplicate_names[taxid] = [name]
                else:
                    duplicate_names[taxid] = duplicate_names[taxid] + [name]

    for taxid, names in duplicate_names.items():
        names = ','.join(names)
        newick_file = newick_file.replace(taxid, names)

    # add other IDs
    for name, taxid in IDs.items():
        if taxid in newick_file:
            newick_file = newick_file.replace(taxid, name)

    # Add other names
    taxids = re.findall(r'\d+', newick_file)
    for taxid in taxids:
        url = "https://lbgi.fr/api/taxonomy/description/{}".format(str(taxid))
        r = requests.get(url, headers={ "Accept" : "application/json"})
        name = r.json()['data']['name']

        newick_file = newick_file.replace(taxid, name)

    # Use newick to reformat the file
    newick_file = newick_file.replace(' ', '-')
    tree = newick.loads(newick_file)[0]
    tree.remove_redundant_nodes(keep_leaf_name=True)

    with open(f'{tree.name}.tree.ascii_art'.lower(), 'w') as file:
        file.write(tree.ascii_art().replace('-', ' '))

    with open(f'{tree.name}.tree'.lower(), 'w') as file:
        file.write(newick.dumps(tree).replace('-', ' '))

if __name__ == "__main__":
    main()
