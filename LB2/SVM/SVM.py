#!/usr/bin/python3

# python3 SVM.py ../LB2/training_set.tsv ../LB2/benchmark_set.tsv


from sys import argv
import pandas as pd
import numpy as np
from sklearn import svm
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics import precision_recall_fscore_support
from sklearn.metrics import accuracy_score
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.ProtParamData import kd

def encode(seq,k):
	# Takes as input K (first kth residues of seq) and the sequence
	aas = "ARNDCQEGHILKMFPSTWYV"
	matrice = np.zeros([20,k])
	for i,res in enumerate(seq[:k]):
		#print (aas.find(res),i)
		matrice[aas.find(res),i]=1
	matrice = np.reshape(matrice,20*k)
	return matrice


def encode1(seq,k):
	# Takes as input K (first kth residues of seq) and the sequence        need correction
	aas = "ARNDCQEGHILKMFPSTWYV"
	matrice = np.zeros([20,k])
	for i,res in enumerate(seq[:k]):
		#print (aas.find(res),i)
		matrice[aas.find(res),i]=1
	matrice = np.reshape(matrice,20*k)
	return matrice

	
def encode2(seq,k):
	aas = "ARNDCQEGHILKMFPSTWYV"
	matrice = np.zeros([20])
	seq = seq[:k]
	for i in range(len(aas)):
		matrice[i] = seq.count(aas[i])/k
	return matrice

	
def encode4(seq,k):
	aas = "ARNDCQEGHILKMFPSTWYV"
	matrice = np.zeros([20+k])
	hydrophobicity = np.zeros([k])
	seq = seq[:k]
	for i in range(20):
		matrice[i] = seq.count(aas[i])/k
	startandend= seq[0]+seq[1]+seq[-2]+seq[-1]
	seq= ProteinAnalysis(seq)
	matrice[20]=kd[startandend[0]]
	matrice[21]=kd[startandend[1]]
	matrice[22:20+k-2]=seq.protein_scale(kd,5)
	matrice[20+k-2]=kd[startandend[-2]]
	matrice[20+k-1]=kd[startandend[-1]]
	return matrice

def scores(true,pred):
	TP,TN,FP,FN=0,0,0,0
	for i in range(len(true)):
		if true[i]==pred[i] and true[i]==1:
			TP += 1
		elif true[i]==pred[i] and true[i]==0:
			TN += 1
		elif true[i]!=pred[i] and true[i]==1:
			FN+=1
		elif true[i]!=pred[i] and true[i]==0:
			FP+=1
	accuracy = (TP+TN)/(TP+TN+FN+FP)
	precision = TP/(TP+FP)
	recall = TP/(TP+FN)
	MCC = (TP*TN - FP*FN) / np.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
	F1 = 2*((precision*recall)/(precision+recall))
	return [accuracy,precision,recall,F1,MCC]

if __name__ == "__main__":
	dataset = argv[1]
	df = pd.read_csv(dataset, sep='\t')
	
	# Vector of labels
	labels = np.where(df['Class']=='SP',1,0)
	
	# Vector of folds
	#print (df.columns)
	folds = df['Cross-validation fold']
	
	
	grids = {}
	all_scores=[]
	for k in [20,22,24]:
		# We add each encoded seq to a new list
		encoded_seqs = []
		for i in range(len(df)):
			# first encoding is positional based and second one is compact
			encoded_seqs.append(encode4(df.loc[i,'Sequence (first 50 N-terminal residues)'],k)) #change here for encoding type
		encoded_seqs = np.array(encoded_seqs)
		
		for g in [0.5,1,"scale"]:
			for c in [1,2,4]:
				mcc = []
				for i in range(5):
					mysvm = svm.SVC(C = c, kernel = 'rbf', gamma = g)
					mysvm.fit(encoded_seqs[folds != i],labels[folds != i])
					y_pred = mysvm.predict(encoded_seqs[folds == i])
					mcc.append(matthews_corrcoef(labels[folds == i],y_pred))
					if k==22 and g=='scale' and c==1:                                        # after running one time, change the grid to appropriate one
						all_scores.append(scores(labels[folds == i],y_pred))
				final_mcc = np.mean(mcc)
				grids[(k,g,c)]=final_mcc
	print (max(grids, key=grids.get), max(grids.values()))
	print (np.mean(all_scores,axis=0))
	arguments = max(grids, key=grids.get)
		
	
	
	# train again with the best hyperparameters
	encoded_seqs = []
	for i in range(len(df)):
		# first encoding is positional based and second one is compact
		encoded_seqs.append(encode4(df.loc[i,'Sequence (first 50 N-terminal residues)'],arguments[0])) #change here for encoding type
	encoded_seqs = np.array(encoded_seqs)
	mysvm = svm.SVC(C = arguments[2], kernel = 'rbf', gamma = arguments[1])
	mysvm.fit(encoded_seqs,labels)
	
	
	
	# Test model with benchmarking
	benchmarking = argv[2]
	df = pd.read_csv(benchmarking, sep='\t')
	labels = np.where(df['Class']=='SP',1,0)
	
	encoded_seqs = []
	for i in range(len(df)):
		encoded_seqs.append(encode4(df.loc[i,'Sequence (first 50 N-terminal residues)'], arguments[0])) #change here for encoding type
	encoded_seqs = np.array(encoded_seqs)
	
	y_pred = mysvm.predict(encoded_seqs)
	print(matthews_corrcoef(labels,y_pred))
	print (scores(labels,y_pred))
	
	df['svm4_prediction'] = y_pred
	df.to_csv('predicted-dataframe.csv', index=None)
	#print (df)
	#all_negatives = df.loc[df['Class']=='NO_SP','UniProtKB accession']
	#all_negatives.to_csv('all_negatives.txt', header=None, index=None, sep=' ', mode='a')
	
	
