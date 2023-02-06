#!/usr/bin/python3

#python3 vHr-predict.py LB2/training_set.tsv results/profile.npy 

from sys import argv
import pandas as pd
import numpy as np
import seaborn as sns
from sklearn.metrics import roc_curve
from sklearn.metrics import RocCurveDisplay
import matplotlib.pyplot as plt


def predict(df, profile):
	'''
	This function add two columns to the dataframe
	1) best position that signal peptide could potentially start
	2) best score associated to that potential signal peptide
	'''
	for i in range(len(df)):
		best_pos_score = [-9999,-99999999999]
		seq = df.iloc[i,]['Sequence (first 50 N-terminal residues)']
		for ii in range(len(seq)-15+1):
			tmp = 0
			for iii in range(ii,15+ii):
				#print(iii, iii-ii,len(seq))
				score = profile[aas.find(seq[iii]), iii-ii]
				tmp += score
			if best_pos_score[1] < tmp:
				best_pos_score = [ii,round(tmp,2)]
		df.loc[i,"best_pos"]=best_pos_score[0]
		df.loc[i,"best_score"]=best_pos_score[1]
	
	
	#print (df[df['Class']=="NO_SP"]["best_score"].mean())
	#print (df[df['Class']=="SP"]["best_score"].mean())
	#print (df.columns)
	#print (df)
	
	return df

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
	F1 = 2*((precision*recall)/(precision+recall))
	#print ('\n\nwith threshold',score,'\nF1: ',F1)
	return score, F1
	
def plot(df, best_threshold):
	plot = sns.stripplot(y =df['best_score'], x=df.Class, hue = df.Class, dodge=True)
	plot.axhline(best_threshold, color = 'dodgerblue', linestyle = '--', label = 'threshold')
	plt.legend()
	
	
	#fpr, tpr, thresholds = roc_curve(df.Class, df['best_score'], pos_label='SP')
	#plt.plot(fpr,tpr)
	#RocCurveDisplay(fpr=fpr, tpr=tpr).plot()
	plt.show()	

if __name__ == "__main__":
	# Read the training data
	dataset = argv[1]
	df = pd.read_csv(dataset, sep='\t')
	
	# Initialize two new columns and split the folds
	df['best_pos']=np.NaN
	df['best_score']=-1
	'''
	df0 = df[df["Cross-validation fold"] == 0]
	df1 = df[df["Cross-validation fold"] == 1]
	df2 = df[df["Cross-validation fold"] == 2]
	df3 = df[df["Cross-validation fold"] == 3]
	df4 = df[df["Cross-validation fold"] == 4]
	dfs = [df0,df1,df2,df3,df4]
	'''
	# Read the profileS
	aas = "ARNDCQEGHILKMFPSTWYV"
	profile = argv[2]
	'''
	threshold_of_folds = []
	profile0 = np.load(profile.replace('.','0.'))
	profile1 = np.load(profile.replace('.','1.'))
	profile2 = np.load(profile.replace('.','2.'))
	profile3 = np.load(profile.replace('.','3.'))
	profile4 = np.load(profile.replace('.','4.'))
	profiles = [profile0,profile1,profile2,profile3,profile4]
	'''
	profile = np.load(profile)
	
	'''
	for each in range(len(profiles)):
		# score the rows in traning data
		predict(dfs[each],profiles[each])
	
		# assess different thresholds
		score=1
		thresholds=['-',0]
		for i in range(1,20):
			score += 0.5
			thre_and_F1 = F1(dfs[each],score)
			if thre_and_F1[1]> thresholds[1]:
				thresholds = thre_and_F1
			
		# Finish
		plot(dfs[each], thresholds[0])
		threshold_of_folds.append(thresholds[0])
		print ("best threshold is:",thresholds)
	print (threshold_of_folds)
	'''
	predict(df,profile)
	score=1
	thresholds=['-',0]
	for i in range(1,20):
		score += 0.5
		thre_and_F1 = F1(df,score)
		if thre_and_F1[1]> thresholds[1]:
			thresholds = thre_and_F1
		
	# Finish
	plot(df, thresholds[0])
	print ("best threshold is:",thresholds)
	
	
	
	
	
