import pandas as pd 
import numpy as np
from sklearn.decomposition import PCA
from sklearn import preprocessing
import matplotlib.pylab as plt
import sys 

	
TMAX=int(sys.argv[1])



pca=PCA(n_components=3) #objeto PCA
data2=pd.read_csv(f'Genes_datos/genes_{TMAX-100}.csv', header = None)
data3=data2.T
scaled_data = preprocessing.scale(data3.T)
pca.fit(scaled_data)


print('corre')
tick=100
times=int(TMAX/tick)
for i in range(times):
	time=i*tick
	data2=pd.read_csv(f'Genes_datos/genes_{time}.csv', header = None)
	data3=data2.T
	scaled_data = preprocessing.scale(data3.T)

	pca_data=pca.transform(scaled_data)

	per_var = np.round(pca.explained_variance_ratio_* 100, decimals=1)
	labels = ['PC' + str(x) for x in range(1,len(per_var)+1)]

	#plt.show()
	plt.close()

	pca_df = pd.DataFrame(pca_data, columns=labels)
	pca_df.to_csv(f'PCA_datos/PCA_eden_{time}.csv',index=False)

	plt.scatter(pca_df.PC1,pca_df.PC2)

	
	'''
	plt.clf()

	fig = plt.figure()
	ax = fig.add_subplot(projection='3d')
	ax.scatter(pca_df.PC1,pca_df.PC2,pca_df.PC3)
	ax.set_xlim(-4,4)
	ax.set_ylim(-4,4)
	ax.set_zlim(-4,4)

	ax.set_xlabel('PC1')
	ax.set_ylabel('PC2')
	ax.set_zlabel('PC3')
	plt.title(f"PCA espacio genetico tiempo: {time}")
	plt.savefig(f'per_var/im{time}')
	plt.clf()
	plt.close()
	'''

