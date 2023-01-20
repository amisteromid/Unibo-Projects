#!/usr/bin/python3

from sys import argv


def parse_fasta(fasta):
	d = {}
	with open(fasta) as f:
		for line in f:
			if line.find(">")==0:
				pid = line.split("|")[1]
			else:
				d[pid]=d.get(pid,'')+line.rstrip()
	return d
	
	


if __name__ == "__main__":
	ids = argv[1]
	fasta = argv[2]
	
	d = parse_fasta(fasta)
	lid = open(ids).read().rstrip().split('\n')
	for pid in lid:
		if pid in d:
			print (">"+pid)
			print (d[pid])
