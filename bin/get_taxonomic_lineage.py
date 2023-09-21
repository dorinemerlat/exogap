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
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--taxid", help="genome's taxid", type=str, required=True)
    parser.add_argument("-n", "--name", help="genome's name", type=str, required=True)
    parser.add_argument("-o", "--output", help="output file", type=str, required=True)
    args = parser.parse_args()

    df = pd.DataFrame(index = ['name', 'taxid'],
                    columns = ['genome', 'superkingdom', 'kingdom', 'phylum', 'subphylum',
                                'class', 'subclass', 'order', 'family', 'genus', 'species'])

    r = requests.get("https://lbgi.fr/api/taxonomy/lineage/{}".format(args.taxid),
                    headers={ "Accept" : "application/json"})

    r = r.json()['data']

    name = {'genome': args.name}
    taxid = {'genome' : args.taxid}
    for i in r:
        if 'rank' in i:
            name[i['rank']] = i['name']
            taxid[i['rank']] = str(i['id'])

    df.loc['name'] = name
    df.loc['taxid'] = taxid

    df.to_csv(args.output, index = False, sep = '\t')

if __name__ == "__main__":
    main()
