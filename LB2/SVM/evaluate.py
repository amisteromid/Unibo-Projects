#!/usr/bin/python3

# python3 evaluate.py annotation-negatives.tsv predicted-dataframe.csv
# python3 evaluate.py annotation-negatives.tsv ../predicted-dataframe-vh.csv 

from sys import argv
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def confusion(df):
	TP,FP,TN,FN = 0,0,0,0
	for i in range(len(df)):
		if df.loc[i,'Class'] == 'SP' and df.loc[i,'svm4_prediction'] == 1:
			TP += 1
		elif df.loc[i,'Class'] == 'SP' and df.loc[i,'svm4_prediction'] == 0:
			FN += 1
		elif df.loc[i,'Class'] != 'SP' and df.loc[i,'svm4_prediction'] == 1:
			FP += 1
		elif df.loc[i,'Class'] != 'SP' and df.loc[i,'svm4_prediction'] == 0:
			TN += 1
	return TP,FP,TN,FN


def TM(annotation):
	auto_annot = ['ECO:0000256','ECO:0000313','ECO:0007829']
	for i in range(len(annotation)):
		#print (annotation.columns)
		if annotation.notnull().iloc[i,1]:
			tms = annotation.iloc[i,1].split('TRANSMEM ')
			for tm in tms[1:]:
				#print (tm)
				if "Helical" in tm and 'evidence' in tm and \
				auto_annot[0] not in tm and auto_annot[1] not in tm and\
				auto_annot[2] not in tm and int(tm.split()[0].split("..")[0]) <= 35:
					annotation.iloc[i,1]=1
				else:
					annotation.iloc[i,1]=0
		else:
			annotation.iloc[i,1]=0
			
def transit(annotation):
	auto_annot = ['ECO:0000256','ECO:0000313','ECO:0007829']
	for i in range(len(annotation)):
		#print (annotation.columns)
		if annotation.notnull().iloc[i,2]:
			tms = annotation.iloc[i,2].split('TRANSIT ')
			for tm in tms[1:]:
				#print (tm)
				if 'ECO' in tm and auto_annot[0] not in tm and auto_annot[1] not in tm and\
				auto_annot[2] not in tm:
					if 'Mitochondrion' in tm:
						annotation.iloc[i,2]='m'
					elif 'Peroxisome' in tm:
						annotation.iloc[i,2]='p'
					elif 'Chloroplast' in tm:
						annotation.iloc[i,2]='c'
					else:
						annotation.iloc[i,2]=1
				else:
					annotation.iloc[i,2]=0
		else:
			annotation.iloc[i,2]=0
	
def length_distribution(df_sp):
	s_dist = []
	for i in range(len(df_sp)):
		if df_sp.iloc[i,]['Class']=='SP' :
			print (df_sp.iloc[i,])
			s_dist.append(df_sp.iloc[i,]['SP cleavage-site annotation'].count('S'))
	
	sns.histplot(data=s_dist).set(xlabel='Length')
	plt.show()
	
	
if __name__ == "__main__":
	# the dataframe with their svm prediction
	df = argv[2]
	df = pd.read_csv(df)
	#print(df)
	TP,FP,TN,FN = confusion(df)
	FPR = FP/(FP+TN)
	#print (FPR)
	length_distribution(df)
	
	'''
	# the annotation dataframe and extracting annotations
	annotation = argv[1]
	annotation = pd.read_csv(annotation,index_col=0,sep='\t')
	#print(annotation['Transit peptide'].describe())
	TM(annotation)
	transit(annotation)
	#print(annotation['Transit peptide'].unique())
	#print (annotation)
	
	
	# New dataframe
	new_df = df[['UniProtKB accession','Class','svm4_prediction']]
	new_df.set_index('UniProtKB accession',inplace=True)
	new_df = new_df.join(annotation[['Transit peptide','Transmembrane']])
	print (new_df)
	
	FP = (new_df['Class']=='NO_SP') & (new_df['svm4_prediction']==1)
	all_N = (new_df['Class']=='NO_SP')
	print ("FPR: ",len(new_df[FP])/len(new_df[all_N]))
	
	
	fptm = (new_df['Class']=='NO_SP') & (new_df['svm4_prediction']==1) & (new_df['Transmembrane']==1)
	tm = (new_df['Transmembrane']==1)
	print ("TM: ",len(new_df[fptm])/len(new_df[tm]))

	fptransit = (new_df['Class']=='NO_SP') & (new_df['svm4_prediction']==1) & (new_df['Transit peptide']!=0)
	transit = (new_df['Transit peptide']!=0)
	print ("TP: ",len(new_df[fptransit])/len(new_df[transit]))
	
	'''
