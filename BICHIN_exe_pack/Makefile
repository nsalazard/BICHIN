all: Basic_analysis.py PCA_Data.py
	
PCA_Data.py: Haldanes.txt | Genes_datos
	python PCA_Data.py

Basic_analysis.py: Haldanes.txt
	python Basic_analysis.py 

Haldanes.txt: sim_bichin.exe
	./sim_bichin.exe


sim_bichin.exe : Simulacion_Bichines.cpp | Genes_datos
	g++ Simulacion_Bichines.cpp -O3 -o sim_bichin

Genes_datos:  
	mkdir Genes_datos


clean: 
	del /S *.png
	del /S *.txt
	del /S *.csv
	del /S *.gp
