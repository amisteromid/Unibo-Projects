#!/usr/bin/python3

from sys import argv
import pandas as pd



if __name__ == "__main__":
	pfam = argv[1]
	pfam = pd.read_csv(pfam, index_col = 3)
	pfam.drop(pfam.columns[0],axis=1,inplace= True)
	pdb = argv[2]
	pdb = pd.read_csv(pdb,index_col = 0).index.tolist()
	#print (len(pdb),len(pfam.index))
	for i in pfam.index:
		if i not in pdb:
			pfam = pfam.drop(i,axis=0)
	#print (len(pfam.index))
	pfam.to_csv("pfam00014.csv.3")
