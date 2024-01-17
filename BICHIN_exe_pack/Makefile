all: Basic_analysis.py PCA_Data.py
	
PCA_Data.py: Haldanes.txt | Genes_datos
	python PCA_Data.py

Basic_analysis.py: Haldanes.txt
	python Basic_analysis.py 

Haldanes.txt: sim_bichin.exe
	./sim_bichin.exe -t 10000


sim_bichin.exe : Simulacion_Bichines.cpp | Genes_datos
	g++ Simulacion_Bichines.cpp -O3 -o sim_bichin

Genes_datos:  
	mkdir Genes_datos

clean:
	powershell -Command "Get-ChildItem -Path . -Include *.png -Recurse | Remove-Item -Force"
	powershell -Command "Get-ChildItem -Path . -Include *.txt -Recurse | Remove-Item -Force"
	powershell -Command "Get-ChildItem -Path . -Include *.csv -Recurse | Remove-Item -Force"
	powershell -Command "Get-ChildItem -Path . -Include *.gp -Recurse | Remove-Item -Force"