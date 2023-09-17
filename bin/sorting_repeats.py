#!/usr/bin/env python

#************************************************#
#               sorting_repeats.py               #
#            written by Dorine Merlat            #
#                 April 01, 2022                 #
#                                                #
#    Sorting classified and unclassified from    #
#       a library  of repetitive elements.       #
#************************************************#

from Bio import SeqIO
import re
import argparse

def sorting(input, unclassified_file, classified_file):
    unclassified = open(unclassified_file, 'w')
    classified = open(classified_file, 'w')

    for seq_record in SeqIO.parse(input, 'fasta'):
        if seq_record.id.split('_')[-1] == "Unknown":
            unclassified.write(seq_record.format('fasta'))
        else :
            classified.write(seq_record.format('fasta'))

    unclassified.close()
    classified.close()


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", help='FASTA file to sort (required)', type=str, required=True, metavar="FASTA")
    parser.add_argument("-c", "--classified", help='FASTA file output with the classified repetitive elements (required)', type=str, required=True, metavar="FASTA")
    parser.add_argument("-u", "--unclassified", help='FASTA file output with the unclassified repetitive elements (required)', type=str, required=True, metavar="FASTA")

    args = parser.parse_args()

    sorting(args.input, args.unclassified, args.classified)


if __name__ == "__main__":
    main()
