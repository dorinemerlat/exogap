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
import newick
import json

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--taxids", help="list of genome's taxids", type=json.loads, required=True)
    parser.add_argument("-i", "--infos", help="info on genomes", type=json.loads, required=True)
    args = parser.parse_args()

    print(args.infos)

    taxids = [str(x) for x in args.taxids]
    taxids = '%2C'.join(taxids)

    r = requests.get("https://lbgi.fr/api/taxonomy/tree/{}".format(taxids),
                    headers={ "Accept" : "application/json"})
    tree = r.json()['data']["tree"]
    tree = newick.loads(tree)[0]
    tree.remove_redundant_nodes(keep_leaf_name=True)

    print(tree)
    # for genome in args.infos:
    #     print(genome)


if __name__ == "__main__":
    main()
