import numpy as np
import pandas as pd


config = pd.read_csv('config.csv')
TMAX=config['TMAX'][0]
tick=config['data_tick'][0]
distribution=str(config['food_dis'][0])
times=int(TMAX/tick)

results_folder='Clustering_results/'





	

	
