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

def get_subgraphs(G):
	subgraphs = []
	for c in nx.connected_components(G):
		subgraph = G.subgraph(c)
		subgraphs.append(subgraph)
	return subgraphs

def node_list(subgraph):
	nodes=[]
	for node in subgraph.nodes():
		nodes.append(node)
	return sorted(nodes)

def split_list_based_on_first_letter(lst):
    dict_letters = {}
    for word in lst:
        if word[0] not in dict_letters:
            dict_letters[word[0]] = [word]
        else:
            dict_letters[word[0]].append(word)
    return list(dict_letters.values())
	
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
		

subgraphs = get_subgraphs(Graph)

starting_height=0

for subgraph in subgraphs:
	
	nodes=node_list(subgraph)
	
	nodes_time_dependent=split_list_based_on_first_letter(nodes)
	

	for time_step in nodes_time_dependent:
		node_number=0
		for node in time_step:
			pos[node]=(pos[node][0],starting_height+node_number)
			node_number+=1



	lengths = [len(sublist) for sublist in nodes_time_dependent]


	
	subgraph_height=max(lengths)
	starting_height+=subgraph_height

	
nx.draw(Graph,pos=pos,with_labels=True)
plt.show()


	

	



