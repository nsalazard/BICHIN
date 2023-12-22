import pandas as pd 
import numpy as np
from sklearn.decomposition import PCA
from sklearn import preprocessing
import matplotlib.pylab as plt
import sys 

df = pd.read_csv('config.csv')
TMAX=df['TMAX'][0]
tick=df['data_tick'][0]

pca=PCA(n_components=3) #objeto PCA
data2=pd.read_csv(f'Genes_datos/genes_{TMAX-tick}.csv', header = None)
data3=data2.T
scaled_data = preprocessing.scale(data3.T) 
pca.fit(scaled_data)

pci=[]
pci_np_i=[]
pci_np_val=[]
for i in range(3):
	carga_genes= pd.Series(pca.components_[i])
	organizado_carga_gen=carga_genes 
	pci.append(organizado_carga_gen.to_frame()) #toda la info se va a esta lista, un espacio con index y otro con el valor
	pci_np_i.append(pci[i].index.to_numpy(copy=True).reshape(-1,1)) #extrae los indices
	pci_np_val.append(pci[i].to_numpy()) #extrae los valores

pci_np_datos=[]
for i in range(3):
	q=np.concatenate((pci_np_i[i],pci_np_val[i]),axis=1)
	pci_np_datos.append(q)


datos=np.stack((pci_np_datos[0],pci_np_datos[1],pci_np_datos[2]))

#Los 3 arrays organizados en una lista
genes_pca=[]
for i in range(8):
	p=(np.absolute(datos[:,i,1]).tolist())
	genes_pca.append(p) #valor de cada gen en cada componente de PCA

barWidth = 0.25
fig = plt.subplots(figsize =(12, 8),dpi=300)


br1 = np.arange(len(genes_pca[0]))

colores=['darkviolet',
		'lightcoral','firebrick','darkred',
		'darkorange',
		'darkgreen','limegreen','lightgreen'
]
separacion=4

x=1


for i in range(8):
	x0=x+i*barWidth
	br=[x0,x0+1*separacion,x0+2*separacion]
	plt.bar(br,genes_pca[i] ,color=colores[i], width = barWidth,edgecolor ='grey', label =f'Gene{i}')


plt.xticks([1.75,5.75,9.75],['PC1','PC2','PC3'])
plt.ylabel('Component Load Scores')
plt.title('Absolute value of the linear decomposition of PCA')
plt.legend()
plt.savefig('dist_PCA.png',bbox_inches='tight')
plt.clf()

print('corre')
times=int(TMAX/tick)
for i in range(times):
	time=i*tick
	data2=pd.read_csv(f'Genes_datos/genes_{time}.csv', header = None)
	data3=data2.T
	scaled_data = preprocessing.scale(data3.T)
	labels2=["gene " + str(x) for x in range(8)]
	
	genes = pd.read_csv(f'Genes_datos/genes_{time}.csv', names=labels2)
	
	pca_data=pca.transform(scaled_data) #Transforma los datos al espacio PCA
	
	per_var = np.round(pca.explained_variance_ratio_* 100, decimals=1) #explained variance ratio da un vector con la varianza de cada dimension
	#print(np.sum(pca.explained_variance_ratio_))
	
	labels = ['PC' + str(x) for x in range(1,len(per_var)+1)]
	
	#plt.show()
	plt.close()
	pca_df = pd.DataFrame(pca_data, columns=labels)
	joined= pd.concat([pca_df,genes],axis=1)
	pca_df.to_csv(f'PCA_datos/PCA_eden_{time}.csv',index=False)
	joined.to_csv(f'PCA_datos_y_genes/PCA_eden_{time}.csv',index=False)
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

