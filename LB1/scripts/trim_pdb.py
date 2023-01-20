#!/usr/bin/python3

from sys import argv
import prody
import pandas as pd




if __name__ == "__main__":
	pfam = argv[1]
	pdbfile = argv[2]
	pfam = pd.read_csv(pfam, index_col = 3)
	pfam.drop(pfam.columns[0],axis=1,inplace= True)
	pdb_id = pdbfile[5:9] #data/1234.pdb
	start = int(pfam.loc[pdb_id,'PDB residues'].split('-')[0])
	end = int(pfam.loc[pdb_id,'PDB residues'].split('-')[1])
	chain_id = pfam.loc[pdb_id,'PDB chain ID']
	
	model = prody.parsePDB(pdbfile, model=1)
	extract = model.select('chain %s resnum %dto%d' % (chain_id,start,end))
	prody.writePDB(pdbfile+'.2', extract)
	
