#!/usr/bin/python3

from sys import argv
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def length_distribution(df_sp):
	s_dist = []
	for i in range(len(df_sp)):
		s_dist.append(df_sp.iloc[i,]['SP cleavage-site annotation'].count('S'))
	
	sns.histplot(data=s_dist).set(xlabel='Length')
	plt.show()

	

def compo(df_sp):
	aa_freq= {}
	for i in range(len(df_sp)):
		length = df_sp.iloc[i,]['SP cleavage-site annotation'].count('S')
		for ii in range(length):
			aa = df_sp.iloc[i,]['Sequence (first 50 N-terminal residues)'][ii] 
			if aa in aa_freq:
				aa_freq[aa] += 1
			else:
				aa_freq[aa] = 1
	summation = sum(aa_freq.values())
	for key in aa_freq:
		aa_freq[key] = round(aa_freq[key]/summation * 100, 2)
		
	
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
		freq = float(each[8:]) 
		aa_freq_rand[aa] = freq
	
	freqs = {}
	for aa in aa_freq:
		freqs[aa] = [aa_freq[aa],aa_freq_rand[aa]]
	df = pd.DataFrame(freqs).T
	df.columns = ['signal sequence', 'random']
	df['aa'] = df.index
	df = pd.melt(df, id_vars=['aa'], value_vars=['signal sequence', 'random'])
	order = [i for i in 'GAVPLIMFWYSTCNQHDEKR']
	sns.catplot(data=df, kind="bar",\
    x="aa", y="value", hue="variable",\
    palette="dark", alpha=.6, height=6,\
    order = order).set(xlabel='Residue type', ylabel='Frequency(%)')
	plt.show()


def species1(df_sp):
	plt.pie(df_sp['Kingdom'].value_counts(),  autopct='%.0f%%')
	plt.legend(labels=list(set(df_sp['Kingdom'])))
	plt.show()
def species2(df_sp):
	al = df_sp['Taxa'].value_counts().sum()
	top_10 = df_sp['Taxa'].value_counts().nlargest(10).sum()
	new_row = pd.Series({ 'other': al - top_10 })
	concatanation = pd.concat([df_sp['Taxa'].value_counts().nlargest(10), new_row])
	#print (concatanation)
	clr = sns.color_palette("Paired",11)
	plt.pie(concatanation,  autopct='%.0f%%', colors = clr)
	plt.legend(labels=concatanation.index, loc='upper left', bbox_to_anchor=(0.9, 1))
	plt.show()
	
def logo(df_sp):
	seqs = []
	for i in range(len(df_sp)):
		seq = df_sp.iloc[i,]['Sequence (first 50 N-terminal residues)']
		length = df_sp.iloc[i,]['SP cleavage-site annotation'].count('S')
		seqs.append(seq[length-13:length+2])
		print (seqs[-1])
	

def merged_dist(dfs): #merge both benchmark and training set plot of distribution of siganl lengths
	dfs['Length']=0
	for i in dfs.index:
		dfs.loc[i,'Length'] = dfs.loc[i,'SP cleavage-site annotation'].count('S')
	#print (dfs)
	sns.kdeplot(dfs["Length"],hue=dfs["Dataset"], shade = True)
	plt.show()
	
	
if __name__ == "__main__":
	dataset = argv[1]
	df = pd.read_csv(dataset, sep='\t')
	df_sp = df[df['Class']=='SP']
	#length_distribution(df_sp)
	#compo(df_sp)
	#species1(df)
	#species2(df)
	#logo(df_sp)
	
	'''
	dataset2 = argv[2]
	df2 = pd.read_csv(dataset2, sep='\t')
	df2_sp = df2[df2['Class']=='SP']
	df2_sp['Dataset']="benchmark"
	df_sp['Dataset']='training'
	dfs = pd.concat([df_sp,df2_sp])
	#print (dfs)
	merged_dist(dfs)
	'''
	
	
	
	
