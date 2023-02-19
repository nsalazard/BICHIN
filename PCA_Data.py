import pandas as pd 
import numpy as np
import random as rd
from sklearn.decomposition import PCA
from sklearn import preprocessing
import matplotlib.pylab as plt


genes=['gene' + str(i) for i in range(0,100)]

wt = ['wt' + str(i) for i in range(0,5)]
ko = ['ko' + str(i) for i in range(0,5)]


data = pd.DataFrame(columns=[*wt,*ko] index=genes )


for gene in data.index:
	data.loc[gene,'wt0':'wt4'] = np.random.poisson( lam= rd.randrange(10,1000), size=5)
	data.loc[gene,'ko0':'ko4'] = np.random.poisson( lam= rd.randrange(10,1000), size=5)









