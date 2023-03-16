import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import math
TMAX=int(sys.argv[1])
ticks=100
times=int(TMAX/ticks)
for i in range(times):
    time=i*ticks
    data = pd.read_csv(f'Datos_Completos/Datos_{time}.csv') #csv with genetic and pca data for each bichin at time t
    genes_clusters = pd.DataFrame(columns=data.columns[[4,5,6,7,8,9,10,11]]) #get genetic data
    n_clust_max=data["Cluster"].max() #number of clusters
    clusters=[]
    for i in range(n_clust_max):
        dfaux=data.loc[data['Cluster'] == i] #get bichines that belong to cluster i 
        avg=[]
        clusters.append(i)
        for i in range(8):
            avg.append(dfaux.iloc[:,4+i].sum()) #sum values for each gene of bichines in cluster i
        bichines_cluster = len(dfaux) 
        avg=[x/bichines_cluster for x in avg]
        genes_clusters.loc[len(genes_clusters)] = avg
        clustersdf=pd.DataFrame(clusters, columns=["Cluster"])
        complete_data=pd.concat([clustersdf,genes_clusters],axis=1)
        complete_data.to_csv(f'Promedios_Cluster/Datos_{time}.csv',index=False)

#This is an attempt at relating the ancestors, under contruction
# calculate the Euclidean distances between all pairs of rows

for i in range(times-1):
    time=i*ticks
    
    if(time>0): 
       
        dataold = pd.read_csv(f'Ancestros/Datos_{time-ticks}.csv')
        datanew = pd.read_csv(f'Promedios_Cluster/Datos_{time+ticks}.csv')
        distances = np.sqrt(((datanew.iloc[:, -8:].values[:, np.newaxis, :] - dataold.iloc[:, -8:].values) ** 2).sum(axis=2))
        # find the row pairs with the smallest distances
        min_distances = np.min(distances, axis=1)
        argmin_distances = np.argmin(distances, axis=1)

        # print the most similar row pairs
        for i, j in enumerate(argmin_distances):      
            closest_cluster = dataold.iloc[j]  # extract the closest row from dold
            element_to_add = closest_cluster[0]
            datanew.at[i, "closest_ancestor"] = element_to_add
            datanew.to_csv(f'Ancestros/Datos_{time}.csv',index=False)
    else:
       dataold = pd.read_csv(f'Promedios_Cluster/Datos_{time}.csv')
       datanew = pd.read_csv(f'Promedios_Cluster/Datos_{time+ticks}.csv')
       distances = np.sqrt(((datanew.iloc[:, -8:].values[:, np.newaxis, :] - dataold.iloc[:, -8:].values) ** 2).sum(axis=2))
    # find the row pairs with the smallest distances
       min_distances = np.min(distances, axis=1)
       argmin_distances = np.argmin(distances, axis=1)

    # print the most similar row pairs
       for i, j in enumerate(argmin_distances):      
            closest_cluster = dataold.iloc[j]  # extract the closest row from dold
            element_to_add = closest_cluster[0]
            datanew.at[i, "closest_ancestor"] = element_to_add
            datanew.to_csv(f'Ancestros/Datos_{time}.csv',index=False)
    


for i in range(times):
    barWidth = 0.5
    time=i*ticks
    data = pd.read_csv(f'Promedios_Cluster/Datos_{time}.csv')
    num_clusters=data["Cluster"].max()
    data.set_index('Cluster', inplace=True) # Set the 'Cluster' column as the index for the DataFrame
    data = data.transpose() # Transpose the DataFrame so that each row represents a gene
    color_map = {0: 'red', 1: 'blue', 2: 'green', 3: 'purple', 4: 'orange', 5: 'pink', 6: 'brown', 7: 'gray'}

    # Create a single plot for all clusters
    if(num_clusters>3):
        fig, ax = plt.subplots(nrows=len(data.columns),sharex=True,figsize=(10, 2*num_clusters))
        fig.suptitle('Measurement by Gene and Cluster t='+str(time))
    else:
        fig, ax = plt.subplots(nrows=len(data.columns),sharex=True,figsize=(10, 7))
        fig.suptitle('Measurement by Gene and Cluster t='+str(time))
    # Loop over each cluster and add its bar plot to the shared axis
    for i, cluster in enumerate(data.columns):
        ax[i].bar(data.index, data[cluster], color=color_map[cluster], width=barWidth )
        ax[i].set_title(f'Cluster {cluster}')

    # Set shared axis labels
    fig.text(0.5, 0.04, 'Genes ', ha='center')
    fig.text(0.04, 0.5, 'Measurement', va='center', rotation='vertical')
    plt.subplots_adjust(hspace=1)
    plt.savefig(f'IMG/gene_avg_Clusters_{time}.png',bbox_inches="tight")
    


