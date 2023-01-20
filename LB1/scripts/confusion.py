#!/usr/bin/python3

from sys import argv
import numpy as np
import math

def get_cm(filename,th=1):
	dic={'TP':0, 'TN':0, 'FN':0, 'FP':0}
	with open(filename) as fh:
		for line in fh:
			line=line.rstrip().split()
			line[1],line[2]=float(line[1]),float(line[2])
			if line[-1] == 1:
				if line[1]<th:
					dic['TP']+=1
				else:
					dic['FP']+=1
			else:
				if line[1]<th:
					dic['FN']+=1
				else:
					dic['TN']+=1
	return dic

def get_accuracy(dic):
	t = dic['TP']+dic['TN']
	a = dic['TP']+dic['TN'] + dic['FP']+dic['FN']
	return t/a

def get_mcc(dic):
	a = dic['TP']*dic['TN'] - dic['FP']*dic['FN']
	b = (dic['TP']+dic['FP']) * (dic['TP']+dic['FN']) * (dic['TN']+dic['FP']) * (dic['TN']+dic['FN'])
	return a/ ( math.sqrt(b))

def main():
	filename=argv[1]
	th = float(argv[2])
	dic=get_cm(filename,th)
	for k,v in dic.items():
		print (k,v)
	accuracy = get_accuracy(dic)
	mcc = get_mcc(dic)
	#specificity = dic['TN']/(dic['TN']+dic['FP'])
	#sensitivity = dic['TP']/(dic['TP']+dic['FN'])
	print ('TH: ',th,'Q2: ',accuracy,'MCC: ',mcc)
	#return 1-specificity, sensitivity
	#return sensitivity
main()
'''
th = [10**(-i) for i in range(29,1,-1)]
#print (th)
for i in list(th):
	print(main(i))
'''
