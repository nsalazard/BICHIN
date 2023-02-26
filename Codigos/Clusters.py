import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances_argmin_min 
from sklearn.metrics import silhouette_score


dataframe = pd.read_csv(r"C:\Users\LENOVO\Documents\2023-1\BICHIN\DATA\PCA_eden_300.csv")
dataframe.head()

X = np.array(dataframe[["PC1","PC2","PC3"]])

# Silhouette Method
range_n_clusters = range(2, 10)
silhouette_avg = []
for num_clusters in range_n_clusters:
 # initialise kmeans
    kmeans = KMeans(n_clusters=num_clusters, n_init= 'auto')
    kmeans.fit(X)
    cluster_labels = kmeans.labels_
 # silhouette score
    silhouette_avg.append(silhouette_score(X, cluster_labels))

plt.plot(range_n_clusters,silhouette_avg)
plt.xlabel("Values of K") 
plt.ylabel("Silhouette score") 
plt.title("Silhouette analysis For Optimal k")
plt.show()

aux= max(silhouette_avg)
index = silhouette_avg.index(aux)
print(index)
if(index == 0):
    n_clus = 2
elif(index == 1):
    n_clus = 3
else:
    n_clus = 1

kmeans = KMeans(n_clusters=n_clus, n_init = 6).fit(X)
centroids = kmeans.cluster_centers_
print(centroids)

# Predicting the clusters
labels = kmeans.predict(X)
# Getting the cluster centers
C = kmeans.cluster_centers_
colores=['purple','cyan', 'red','blue']
asignar=[]
for row in labels:
    asignar.append(colores[row])
 
fig2 = plt.figure()
ax = fig2.add_subplot(projection='3d')
ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=asignar,s=5)
ax.scatter(C[:, 0], C[:, 1], C[:, 2], marker='*', c="pink", s=20)

ax.set_title("Data Analysis")
ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.set_zlabel('PC3')
plt.savefig('IMG/C00.png')
plt.show()
