# Argument 1 : le gff à reformater
# Argument 2 : le gff reformaté

import re
import sys

file_in_name = sys.argv[1]
file_out_name = sys.argv[2]

file_in = open(file_in_name, "r")
lines = file_in.readlines()

file_out = open(file_out_name, "w")

gene_dict={}
mrna_dict={}
count=0

for line in lines :
    if re.search(r'\tgene\t', line) :
        count += 1

        # On récupère la dernière colonne du GFF
        end_line_regex = re.compile (".[^\t]+\n")
        end_line_old = end_line_regex.findall(line)
        end_line_old = end_line_old[0]

        # On modifie la ligne
        end_line_new = '\tg' + str(count) + '\n'
        line = line.replace(end_line_old, end_line_new)
        file_out.write(line)

    elif re.search(r'\tmRNA\t', line) :
        # On récupère la dernière colonne du GFF
        end_line_regex = re.compile (".[^\t]+\n")
        end_line_old = end_line_regex.findall(line)
        end_line_old = end_line_old[0]

        # On modifie la ligne
        end_line_new = '\tg' + str(count) + '-t1\n'
        line = line.replace(end_line_old, end_line_new)
        file_out.write(line)

    elif re.search(r'\t(three_prime_UTR|three_prime_UTR|CDS|exon)\t', line) :
        # On récupère la dernière colonne du GFF
        end_line_regex = re.compile (".[^\t]+\n")
        end_line_old = end_line_regex.findall(line)
        end_line_old = end_line_old[0]

        # On modifie la ligne
        end_line_new = '\ttranscript_id   g' + str(count) + '-t1;  gene_id g' + str(count) +'\n'
        line = line.replace(end_line_old, end_line_new)
        file_out.write(line)

file_in.close()
file_out.close()

