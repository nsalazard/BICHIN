# i dont like makefiles

time=20000
seed=2
clus= ./IMG/Clusters_0.png
config=2

all: Stv_genes.png promedios_genes.png $(clus) tasa_evolutiva.png dist_PCA.png


dist_PCA.png:
	python PCA_Data.py

tasa_evolutiva.png: haldanes.py
	python haldanes.py

$(clus):Clusters.py
	python .\Clusters.py $(time)

Clusters.py: PCA_datos\PCA_eden_0.csv PCA_Data.py
	python PCA_Data.py $(time)
PCA_datos\PCA_eden_0.csv:
	python PCA_Data.py $(time)
Stv_genes.png : Haldanes.txt haldanes.py
	python haldanes.py

Haldanes.txt: sim_bichin.exe
	./sim_bichin.exe $(config) $(seed) $(time)

sim_bichin.exe : Simulacion_Bichines.cpp | Genes_datos PCA_datos IMG
	g++ Simulacion_Bichines.cpp -O3 -o sim_bichin

Genes_datos:  
	mkdir Genes_datos
PCA_datos: 
	mkdir PCA_datos
per_var: 
	mkdir per_var
IMG:
	mkdir IMG

clean: 
	del /S *.png
	del /S *.txt
	del /S *.csv
	del /S *.gp
	rmdir IMG
	rmdir Genes_datos
	rmdir PCA_datos
