#!/usr/bin/python3

from sys import argv

def get_best_evalue(filename,pos=1,e=100):
	d={}
	with open(filename) as f:
		for line in f:
			if line.find('#')==0: continue
			v=line.split()
			d[v[0]] = min(d.get(v[0],e),float(v[pos]))
	return d
			
			
			
if __name__=="__main__":
	filename = argv[1]
	d = get_best_evalue(filename)
	for k in d.keys():
		print (k,d[k])
