import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import os

def check_folder(path):
	if not os.path.exists(path):
		os.makedirs(path)

def find_silhouette_score(PCA_position):
	range_n_clusters = range(1, 10)
	silhouette_avg = []
	centroids_list=[]
	for n_clusters in range_n_clusters: #This for loop runs k-means for 2-10 clusters, the better segmentation is the one with the highest silhouette score

		kmeans = KMeans(n_clusters=n_clusters, n_init= 100)
		kmeans.fit(PCA_position)
		cluster_labels = kmeans.labels_ 
		print(cluster_labels)
		centroids = kmeans.cluster_centers_
		print(centroids)
		if len(set(cluster_labels)) == 1:
			silhouette_avg.append(1)
			centroids_list.append(centroids[0])
			return silhouette_avg,centroids_list
		else:
			silhouette_avg.append(silhouette_score(PCA_position, cluster_labels))
			centroids_list.append(centroids)

	return silhouette_avg,centroids_list

range_n_clusters = range(2, 10)

results_folder='PCA_results/'
clustering_folder='Clustering_results/'



config = pd.read_csv('config.csv')
TMAX=config['TMAX'][0]
tick=config['data_tick'][0]
distribution=str(config['food_dis'][0])
times=int(TMAX/tick)

for i in range(times):
	time=i*tick
	print(time)
	current_data=pd.read_csv(results_folder+f"PCA_datos/PCA_{distribution}_{time}.csv")
	PCA_position = np.array(current_data[["PC1","PC2","PC3"]])

	silhouette_avg,centroids_list = find_silhouette_score(PCA_position)
	
	index = silhouette_avg.index(max(silhouette_avg)) #This is the clustering with the best silhouette score
	
	
	if(index >= 0): #If the best silhouette score is positive, then the clustering is better than random
		n_clus = 2+index
		centroids=centroids_list[index]

	else:
		n_clus = 1
		centroids=centroids_list[0]


	kmeans = KMeans(n_clusters=n_clus, n_init = 100).fit(PCA_position)

	labels = kmeans.predict(PCA_position)
	labels_dataframe=pd.DataFrame(labels, columns=["Cluster"])
	gene_data=pd.read_csv(results_folder+f"PCA_datos_y_genes/PCA_{distribution}_{time}.csv")
	complete_data=pd.concat([labels_dataframe,gene_data],axis=1)
	complete_data=complete_data.sort_values(by=['Cluster'])
	complete_data.to_csv(clustering_folder+f"Clustering_{distribution}_{time}.csv", index=False)


	np.savetxt(clustering_folder+f"centroids_{distribution}_{time}.csv", centroids, delimiter=",")