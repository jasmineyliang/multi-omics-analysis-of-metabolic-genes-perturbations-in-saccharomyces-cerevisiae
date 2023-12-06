
import numpy as np
import pandas as pd

# this script filters out the genes with expression value = 0 in all samples
# read data from txt
df = pd.read_csv("transcriptome.txt", sep="\t")
print(df)
df2 = df.loc[df.max(axis=1)>50]
# print(df.iloc[:,5:])
# df2 = df.loc[(df.iloc[:,5:]!=0).any(axis=1)]
print(df2)

df2.to_csv('filtered_transcriptome.txt', sep='\t', mode='a', index=None)

# drop genes with all zeros expression
# df = df.loc[(df.iloc[:,5:-1]!=0).any(axis=1)]






