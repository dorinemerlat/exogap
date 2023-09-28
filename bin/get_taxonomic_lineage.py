#!/usr/bin/env python

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

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--taxid", help="genome's taxid", type=str, required=True)
    parser.add_argument("-f", "--fasta", help="fasta", type=str, required=True)
    args = parser.parse_args()

    taxonomy = {'file': args.fasta, 'superkingdom': { }, 'kingdom': { }, 'phylum': { },
    'subphylum': { },'class': { }, 'subclass': { }, 'order': { }, 'family': { },
    'genus': { }, 'species': { }}

    r = requests.get("https://lbgi.fr/api/taxonomy/lineage/{}".format(args.taxid),
                    headers={ "Accept" : "application/json"})

    r = r.json()['data']

    for i in r:
        if 'rank' in i:
            taxonomy[i['rank']] = {'name': i['name'], 'taxid': i['id']}

    taxonomy = str(taxonomy).replace("'",'"')
    print(taxonomy)

if __name__ == "__main__":
    main()
