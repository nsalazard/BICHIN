import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances_argmin_min 
from sklearn.metrics import silhouette_score
import sys 

def encontraridx(array, value):# esta funcion empareja un valor al valor mas cercano en una lista dada y retorna el indice de la pareja encon
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


TMAX=int(sys.argv[1])

ticks=100
times=int(TMAX/ticks)

for i in range(times):
	time=i*ticks
	print(time,f'{(time/TMAX)*100}%')
	dataframe = pd.read_csv(f"PCA_datos/PCA_eden_{time}.csv")


	X = np.array(dataframe[["PC1","PC2","PC3"]])

	# Silhouette Method
	range_n_clusters = range(2, 10)
	silhouette_avg = []
	for num_clusters in range_n_clusters:
	# initialise kmeans
		kmeans = KMeans(n_clusters=num_clusters, n_init= 60)
		kmeans.fit(X)
		cluster_labels = kmeans.labels_
	# silhouette score
		silhouette_avg.append(silhouette_score(X, cluster_labels))

	plt.plot(range_n_clusters,silhouette_avg)
	plt.xlabel("Values of K") 
	plt.ylabel("Silhouette score") 
	plt.title(f"Silhouette analysis For Optimal k {time}")
	plt.clf()
	plt.close()

	aux= max(silhouette_avg)
	index = silhouette_avg.index(aux)
	#print(index)

	if(index >= 0):
		n_clus = 2+index

	else:
		n_clus = 1

	kmeans = KMeans(n_clusters=n_clus, n_init = 6).fit(X)
	centroids = kmeans.cluster_centers_
	#print(centroids)

	# Predicting the clusters
	labels = kmeans.predict(X)
	# Getting the cluster centers
	C = kmeans.cluster_centers_
	
	colores=['purple','cyan', 'red','blue']
	asignar=[]
	norm=np.double(4)
	for row in labels:

		color=(np.absolute(C[row]+norm)/np.max(np.absolute(C[row]+norm)))
		max_idx=encontraridx(color,np.max(color))
		color[max_idx]=color[max_idx]*2
		color=color/np.max(color)
		if np.array_equal(color,np.ones((3))):
			color=color/5
		asignar.append(color)
		
	#print('||||||||||||||||||')
	#print(C+norm)
	fig2 = plt.figure()
	ax = fig2.add_subplot(projection='3d')
	#print('X:',type(X),'asignar:',len(asignar))
	ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=asignar,s=5)
	ax.scatter(C[:, 0], C[:, 1], C[:, 2], marker='*', c="pink", s=20)

	ax.set_title(f"Data Analysis tiempo: {time}")
	ax.set_xlabel('PC1')
	ax.set_ylabel('PC2')
	ax.set_zlabel('PC3')
	ax.set_xlim3d(-4, 4)
	ax.set_ylim3d(-4, 4)
	ax.set_zlim3d(-4, 4)
	plt.savefig(f'IMG/Clusters_{time}.png')
	plt.clf()
	plt.close()
	#plt.show()
