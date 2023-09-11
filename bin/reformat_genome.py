#!/usr/bin/env python

from Bio import SeqIO
import re
import json
import argparse
import os

def change_seq_names(fasta_in, outfiles, prefix):

    # Names of out files
    if outfiles is None:
        extension = re.search(r'(?P<name>.fa|fasta)', fasta_in).group('name')
        fasta_out = fasta_in.replace(extension, '_rename.fa')
        json_out = fasta_in.replace(extension, '_rename.json')
    else:
        fasta_out = outfiles + '.fa'
        json_out = outfiles + '.json'

    dictionnary = {}
    num_seq = 1

    with open(fasta_out, 'w') as file :
        for record in SeqIO.parse(fasta_in, 'fasta') :
            if record.id not in dictionnary :
                new_name = prefix + str(num_seq)
                num_seq += 1
                dictionnary[record.id] = new_name
            record.id = dictionnary[record.id]
            record.description = ""
            file.write(record.format('fasta'))

    with open(json_out, 'w') as file:
        json.dump(dictionnary, file)

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help='FASTA file to reformat', type=str, required=True, metavar="FILE")
    parser.add_argument("-o", "--output", help='Output name (without .extension)', type=str, default = None, metavar="FILE")
    parser.add_argument("-p", "--prefix", help='Prefix (by default: SEQ)', type=str, default = 'SEQ', metavar="FILE")
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        print("Input file doesn't exist.")
        exit()


    change_seq_names(args.input, args.output, args.prefix)

    print('\nReformating done\n')

if __name__ == "__main__":
    main()
