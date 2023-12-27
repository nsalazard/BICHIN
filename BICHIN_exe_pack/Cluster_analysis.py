import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import igraph as ig


def find_distance(point1, point2):
    return np.linalg.norm(point1 - point2)

def load_centroids(times):
	centroid_list=[]
	for i in range(times):
		time=i*tick
		centroids=np.loadtxt(results_folder+f"centroids_{distribution}_{time}.csv",delimiter=',')
		centroid_list.append(centroids)

	return centroid_list


def find_ancestor(centroid,prev_centroids):
	distances=[]
	for ancestor in prev_centroids:
		distance=find_distance(centroid,ancestor)
		distances.append(distance)

	return np.argmin(distances),(1/(np.min(distances)+0.1))



	
	
Graph=nx.Graph()


config = pd.read_csv('config.csv')
TMAX=config['TMAX'][0]
tick=config['data_tick'][0]
distribution=str(config['food_dis'][0])
times=int(TMAX/tick)

results_folder='Clustering_results/'



centroid_list=load_centroids(times)
original_centers=centroid_list[0]
prev_centroids=original_centers



pos = {}

for time in range(times):
	centroids=centroid_list[time]
	ii=0
	for center in centroids:
		node=f"{time}_{ii}"
		Graph.add_node(node)
		pos[node] = (time, ii)
		if time>0:
			ancestor_index,strenght=find_ancestor(center,prev_centroids)
			Graph.add_edge(f"{time-1}_{ancestor_index}",node, weight=strenght)
		ii+=1

	prev_centroids=centroids
		

weights = [Graph[u][v]['weight'] for u, v in Graph.edges()]
pos=nx.spring_layout(Graph)
# # print(pos)

# nx.draw(Graph,pos,with_labels=True)
# plt.show()


subgraphs = nx.connected_components(Graph)

# Create a new igraph graph for each subgraph and create a layout for each
for i, nodes in enumerate(subgraphs):
	subgraph = Graph.subgraph(nodes)  # Create a subgraph from the nodes
	G_ig = ig.Graph.from_networkx(subgraph)
	G_ig.vs["label"] = list(subgraph.nodes())
	
	layout = G_ig.layout_reingold_tilford(mode="in", root=[0])
	plot = ig.plot(G_ig, layout=layout)
	plot.save(f'graph_{i}.png')