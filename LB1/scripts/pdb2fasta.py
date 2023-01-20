#!/usr/bin/python3

import sys
from Bio import SeqIO

PDBFile = sys.argv[1]
with open(PDBFile, 'r') as pdb_file:
    for record in SeqIO.parse(pdb_file, 'pdb-atom'):
        print('>' + record.id)
        print(record.seq)
