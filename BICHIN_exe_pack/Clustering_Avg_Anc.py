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

def Distance(M_old, M_actual, n_old, n_clus):
	#create comparison matrix -> #Num Filas - n_clus, #Num columnas - n_old
	#Recordar la menor distancia en todo el sistema
	M_Dis = np.zeros([n_clus, n_old])
	ancestors =  np.zeros(n_clus)
	for ii in range(n_old):
		M_oaux = M_old[ii,0:8]
		for jj in range(n_clus):
			M_acaux = M_actual[jj,0:8]
			M_Dis[ii][jj] = np.sum(np.abs(M_acaux- M_oaux))
	for kk in range(n_clus*n_old):
		mini = np.argmin(M_Dis)
		r = mini//n_clus
		c = mini-(n_clus*row)
		if(ancestors[r] == 0 and M_Dis[r][c] != 1000 and c not in ancestors): 
			ancestors[r] = c
		M_Dis[r][c] = 1000
	return ancestors

n_old=0
for i in range(times):
	time=i*ticks
	print(time,f'{(time/TMAX)*100}%')
	DF_PCA = pd.read_csv(f"PCA_datos/PCA_eden_{time}.csv")
	DF_BICHINES = pd.read_csv(f"Genes_datos/genes_{time}.csv", names=["F","L1","L2","L3","B","R3", "R2", "R1"]) # Ni idea de tus carpetas


	PCA = np.array(DF_PCA[["PC1","PC2","PC3"]])
	BICHINES = np.array(DF_BICHINES[["F","L1","L2","L3","B","R3", "R2", "R1"]])

	# Silhouette Method
	range_n_clusters = range(2, 10)
	silhouette_avg = []
	for num_clusters in range_n_clusters:
	# initialise kmeans
		kmeans = KMeans(n_clusters=num_clusters, n_init= 60)
		kmeans.fit(PCA)
		cluster_labels = kmeans.labels_
	# silhouette score
		silhouette_avg.append(silhouette_score(PCA, cluster_labels))

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

	kmeans = KMeans(n_clusters=n_clus, n_init = 6).fit(PCA)
	centroids = kmeans.cluster_centers_
	#print(centroids)

	# Predicting the clusters
	labels = kmeans.predict(PCA)
	# Getting the cluster centers
	C = kmeans.cluster_centers_
	
	colores=['purple','cyan', 'red','blue']
	asignar=[]
	norm=np.double(4)
	Gen_Avg = np.zeros([n_clus,10]) # El último elemento corresponde a su antepasado
	num = 0

	for row in labels:

		color=(np.absolute(C[row]+norm)/np.max(np.absolute(C[row]+norm)))
		max_idx=encontraridx(color,np.max(color))
		color[max_idx]=color[max_idx]*2
		color=color/np.max(color)
		if np.array_equal(color,np.ones((3))):
			color=color/5
		asignar.append(color)
		# Genetic Average
		for ii in range(8):
			Gen_Avg[row][ii] += BICHINES[num, ii]
		Gen_Avg[row][9] += 1 # Guarda el numero de bichines en ese cluster
	print(Gen_Avg)
	print("holi")
	print(BICHINES[1,1])
	#Divide en el numero total de bichines por cluster
	for jj in range(n_clus):
		for kk in range(8):
			Gen_Avg[jj][kk] /= Gen_Avg[jj][9]

	# La funcion Distance va a devolver una matriz donde los indices representan los clusters actuales
	# El numero en ese indice es el numero del cluster ancestro
	#if time > 0:
	#	Distance(Gen_Avg_OLD, Gen_Avg, n_old, n_clus)
	#if(i>0):
		#Distance(Gen_Avg_OLD,Gen_Avg, n_old, n_clus)
	Gen_Avg_OLD = np.copy(Gen_Avg)
	n_old = n_clus
	#Distance(Gen_Avg_OLD, Gen_Avg, n_old, n_clus)
	

		
		
	#print('||||||||||||||||||')
	#print(C+norm)
	fig2 = plt.figure()
	ax = fig2.add_subplot(projection='3d')
	#print('PCA:',type(PCA),'asignar:',len(asignar))
	ax.scatter(PCA[:, 0], PCA[:, 1], PCA[:, 2], c=asignar,s=5)
	ax.scatter(C[:, 0], C[:, 1], C[:, 2], marker='*', c="pink", s=20)

	ax.set_title(f"Data Analysis time: {time}")
	ax.set_xlabel('PC1')
	ax.set_ylabel('PC2')
	ax.set_zlabel('PC3')
	ax.set_xlim3d(-4, 4)
	ax.set_ylim3d(-4, 4)
	ax.set_zlim3d(-4, 4)
	plt.savefig(f'IMG/Clusters_{time}.png',bbox_inches='tight')
	plt.clf()
	plt.close()
	#plt.show()
	
