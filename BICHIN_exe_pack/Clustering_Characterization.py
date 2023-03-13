import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
TMAX=int(sys.argv[1])
ticks=100
times=int(TMAX/ticks)
for i in range(times):
    time=i*ticks
    data = pd.read_csv(f'Datos_Completos/Datos_{time}.csv')

    genes_clusters = pd.DataFrame(columns=data.columns[[4,5,6,7,8,9,10,11]])
    n_clust_max=data["Cluster"].max()
    clusters=[]
    for i in range(n_clust_max):
        dfaux=data.loc[data['Cluster'] == i]
        avg=[]
        clusters.append(i)
        for i in range(8):
            avg.append(dfaux.iloc[:,4+i].sum())
        bichines_cluster = len(dfaux)
        avg=[x/bichines_cluster for x in avg]
        genes_clusters.loc[len(genes_clusters)] = avg
#dfaux.head()
        clustersdf=pd.DataFrame(clusters, columns=["Cluster"])
        complete_data=pd.concat([clustersdf,genes_clusters],axis=1)
        complete_data.to_csv(f'Promedios_Cluster/Datos_{time}.csv',index=False)


# calculate the Euclidean distances between all pairs of rows
for i in range(times-1):
    time=i*ticks
    if(times==0):  
        dataold = pd.read_csv(f'Promedios_Cluster/Datos_{time}.csv')
    else:
       dataold = pd.read_csv(f'Ancestros/Datos_{time}.csv')
    datanew = pd.read_csv(f'Promedios_Cluster/Datos_{time+ticks}.csv')
    distances = np.sqrt(((datanew.iloc[:, -8:].values[:, np.newaxis, :] - dataold.iloc[:, -8:].values) ** 2).sum(axis=2))
# find the row pairs with the smallest distances
    min_distances = np.min(distances, axis=1)
    argmin_distances = np.argmin(distances, axis=1)

# print the most similar row pairs
    for i, j in enumerate(argmin_distances):      
        closest_cluster = dataold.iloc[j]  # extract the closest row from df1
        element_to_add = closest_cluster[0]
        datanew.at[i, "closest_ancestor"] = element_to_add
        datanew.to_csv(f'Ancestros/Datos_{time}.csv',index=False)



barWidth = 0.25

colores=['darkviolet',
         'lightcoral','firebrick','darkred',
         'darkorange',
         'darkgreen','limegreen','lightgreen'
]
separacion=4

x=1

data = pd.read_csv(f'Promedios_Cluster/Datos_{800}.csv')

values = [[row[i] for i in range(1, 9)] for index, row in data.iterrows()]
values2=[]
for i in range(8):
    aux=[]
    for j in range(len(data)):
        aux.append(values[j][i])
    values2.append(aux)


# Create a bar plot with each row as a separate set of bars
fig, ax = plt.subplots(figsize=(10, 6))
#ax.bar(br1, values2[0], color=colores[i],width=barWidth, label=f"Gene {0}")
br1 = np.arange(len(values2[0]))  
print(values2[0]) 
#for i in range(len(values)):
x0=x+i*barWidth
br=[x0,x0+1*separacion]
for j in range(8):
    ax.bar(br, values2[j],color=colores[j], width=0.1, edgecolor ='grey', label=f"Gene {j}")
    

# Set axis labels and legend
ax.set_xlabel("Clusters")
ax.set_ylabel("Values")
ax.legend()

# Show the plot
plt.show()

   