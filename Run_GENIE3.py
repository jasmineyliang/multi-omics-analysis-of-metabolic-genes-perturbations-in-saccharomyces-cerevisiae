
import numpy as np
import pandas as pd
from GENIE3 import *

# read data from txt
df = pd.read_csv("proteome.txt", sep="\t")
# drop genes with all zeros expression
df = df.iloc[:,0:-3]
df = df.loc[(df.iloc[:,7:]!=0).any(axis=1)]
# df = df.loc[df.max(axis=1)>50]
#transform into array
dfArray = df.to_numpy()
# transpose row and col (GENIE3 row = condition, col = gene)
Transcrip = np.transpose(dfArray)

# extract data only
TranscripData = Transcrip[7:,:]
print (TranscripData)

# extract genename
TranscripName = Transcrip[1].tolist()
print (TranscripName)

#GENIE3
VIM = GENIE3(TranscripData, K='all', nthreads=1000)
get_link_list(VIM, gene_names=TranscripName, file_name = 'ProteomicGraph.tsv')
