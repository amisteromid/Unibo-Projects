#!/usr/bin/python3

# python3 vH-train.py LB2/training_set.tsv LB2/benchmark_set.tsv 

from sys import argv
import pandas as pd
import numpy as np


def make_profile(df_sp):
	# Order of Aminoacids are the one below
	aas = "ARNDCQEGHILKMFPSTWYV"
	#profile = np.ones([20,6])
	profile = np.ones([20,15])
	for i in range(len(df_sp)):
		seq = df_sp.iloc[i,]['Sequence (first 50 N-terminal residues)']
		length = df_sp.iloc[i,]['SP cleavage-site annotation'].count('S')
		seq = seq[length-13:length+2]
		for ii in range(15):
			aa_pos = aas.find(seq[ii])
			profile[aa_pos,ii] += 1
	profile = profile/ (20+len(df_sp))
	
	random_compo='''
	Ala (A) 8.25   Gln (Q) 3.93   Leu (L) 9.65   Ser (S) 6.64
	Arg (R) 5.53   Glu (E) 6.72   Lys (K) 5.80   Thr (T) 5.35
	Asn (N) 4.06   Gly (G) 7.07   Met (M) 2.41   Trp (W) 1.10
	Asp (D) 5.46   His (H) 2.27   Phe (F) 3.86   Tyr (Y) 2.92
	Cys (C) 1.38   Ile (I) 5.91   Pro (P) 4.74   Val (V) 6.86
	'''
	aa_freq_rand = {}
	random_compo = random_compo.replace('\n','').replace('\t','   ').split('   ')[1:-1]
	for each in random_compo:
		aa = each[5]
		freq = float(each[8:]) /100
		aa_freq_rand[aa] = freq
	
	for row in range(len(profile)):
		profile[row,:] = profile[row,:]/ aa_freq_rand[aas[row]]
	profile = np.around(np.log2(profile),2)
	#print (profile)
	return profile
	
	
def predict(dff, profile):
	'''
	This function add two columns to the dataframe
	1) best position that signal peptide could potentially start
	2) best score associated to that potential signal peptide
	'''
	for i in dff.index:
		best_pos_score = [-9999,-99999999999]
		seq = dff.loc[i,]['Sequence (first 50 N-terminal residues)']
		for ii in range(len(seq)-15+1):
			tmp = 0
			for iii in range(ii,15+ii):
				#print(iii, iii-ii,len(seq))
				score = profile[aas.find(seq[iii]), iii-ii]
				tmp += score
			if best_pos_score[1] < tmp:
				best_pos_score = [ii,round(tmp,2)]
		dff.loc[i,"best_pos"]=best_pos_score[0]
		dff.loc[i,"best_score"]=best_pos_score[1]
	
	#print (dff[dff['Class']=="NO_SP"]["best_score"].mean())
	#print (dff[dff['Class']=="SP"]["best_score"].mean())
	
	return dff
	
	
	
def F1(df,score):
	'''
	This function takes 
	1) Dataframe with with scores and classes
	2) A range of numbers 
	to assess scores as a threshold for seperating the classes
	and finally shows us the F1 score assiciated to each threshold
	
	'''
	TP = (df['Class']=='SP') & (df['best_score']>=score)
	TP = len(df.loc[TP,])
	TN = (df['Class']=='NO_SP') & (df['best_score']<=score)
	TN = len(df.loc[TN,])
	FP = (df['Class']=='NO_SP') & (df['best_score']>=score)
	FP = len(df.loc[FP,])
	FN = (df['Class']=='SP') & (df['best_score']<=score)
	FN = len(df.loc[FN,])
	accuracy = (TP+TN)/(TP+TN+FN+FP)
	precision = TP/(TP+FP)
	recall = TP/(TP+FN)
	MCC = (TP*TN - FP*FN) / np.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
	F1 = 2*((precision*recall)/(precision+recall))
	print ("for score of",score,accuracy,precision,recall,F1,MCC)
	#print ('\n\nwith threshold',score,'\tF1: ',F1)
	return score, F1



if __name__ == "__main__":
	aas = "ARNDCQEGHILKMFPSTWYV"
	
	dataset = argv[1]
	df = pd.read_csv(dataset, sep='\t')
	df_sp = df[df['Class']=='SP']
	
	# Initialize two new columns and split the folds
	df['best_pos']=np.NaN
	df['best_score']=-1000
	
	# loop through all folds to come up with an average threshold score
	list_of_thresholds=[]
	for i in range(5):
		df_tmp_train = df_sp[df_sp["Cross-validation fold"] != i]
		profile = make_profile(df_tmp_train)
		
		df_tmp_test = df[df["Cross-validation fold"] == i]
		predict(df_tmp_test,profile)
		
		score=1
		thresholds=['-',0]
		for i in range(1,20):
			score += 0.5
			thre_and_F1 = F1(df_tmp_test,score) 
			if thre_and_F1[1]> thresholds[1]:
				thresholds = thre_and_F1
		print (thresholds)
		list_of_thresholds.append(thresholds[0])
	final_threshold = np.mean(list_of_thresholds)
	print (final_threshold)
		
	# make profile with the whole training data again!
	profile = make_profile(df_sp)
	predict(df,profile)
	print ('Traininging F1 score: ',F1(df,final_threshold)[1])
	
	# And evaluate with benchmarking data
	dataset = argv[2]
	df = pd.read_csv(dataset, sep='\t')
	predict(df,profile)
	print ('Benchmarking F1 score: ',F1(df,final_threshold)[1])
	
	
	# to a predicted dataframe
	for i in df.index:
		if df.loc[i,'best_score'] > final_threshold:
			df.loc[i,'svm4_prediction'] = 1
		else:
			df.loc[i,'svm4_prediction'] = 0
	#print (df)
	df.to_csv('predicted-dataframe-vh.csv', index=None)
