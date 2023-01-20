#!/usr/bin/python3

from sys import argv
import pandas as pd
import seaborn as sns
sns.set_style("whitegrid")
import matplotlib.pyplot as plt

'''
df = pd.read_csv(argv[1], sep=' ',names = ['ID','e-value','label'] ,header = None)
df.loc[df['e-value']==10, 'e-value'] = 1
plot = sns.stripplot(y =df['e-value'], x=df.label, hue = df.label, dodge=True)
plot.set(yscale='log', ylim=(10e-30,10))
plot.axhline(1.75e-9, color = 'dodgerblue')
'''
df2 = pd.read_csv('overall_th.txt.2',header=None, names = ['1-specificity','sensitivity'])
sns.scatterplot(x=df2['1-specificity'],y=df2['sensitivity'])


plt.show()



