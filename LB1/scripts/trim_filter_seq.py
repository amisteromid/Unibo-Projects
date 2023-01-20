#!/usr/bin/python3

#python3 scripts/trim_filter_seq.py pfam00014.csv.2 data/5PTI.fasta 

from sys import argv
import pandas as pd
import seaborn as sns
sns.set_style("whitegrid")
import matplotlib.pyplot as plt



def trim_fasta(fasta,pfam):
	# This function gets a fasta file and a table for domain information, Then it returns the fasta file containing only the domain. (filter domain: length <40)
	pdb_id = fasta[5:9] #data/1234.fasta
	start = int(pfam.loc[pdb_id,'PDB residues'].split('-')[0])-1
	end = int(pfam.loc[pdb_id,'PDB residues'].split('-')[1])-1
	chain_id = pfam.loc[pdb_id,'PDB chain ID']
	
	with open(fasta) as fh:
		key = False
		for line in fh:
			if line[0] == '>' and line.split('|')[1][1:].find(chain_id) > -1:
				key = True
				info = line
				seq=''
			elif key == True:
				seq += line.rstrip()
			elif line[0] == '>' and line.split('|')[1][1:].find(chain_id) == -1:
				key = False
		seq=seq[start:end+1]
	if len(seq) < 40:
		return
	else:
		return info + seq



if __name__ == "__main__":
	pfam = argv[1]
	pfam = pd.read_csv(pfam, index_col = 0)
	
	# Length of domain
	'''
	lengths = pfam.loc[:,pfam.columns[-1]].tolist()
	lengths = list(map(lambda x: int(x.split('-')[1])-int(x.split('-')[0]),lengths))
	sns.displot(data=lengths)
	plt.show()
	'''
	
	# Trim fasta files
	fasta = argv[2]
	new_fasta = trim_fasta(fasta,pfam)
	
	print(new_fasta)

			
