<div align="center">

# BICHIN
[![License: CC BY 4.0](https://licensebuttons.net/l/by/4.0/80x15.png)](https://creativecommons.org/licenses/by/4.0/)
Author1, Author2
<span style="color:red">B</span>iosimulación <span style="color:red">I</span>ntegrada de <span style="color:red">C</span>riaturas en <span style="color:red">H</span>ábitats <span style="color:red">IN</span>formáticos

</div>

This is a simulation that shows the effect of the environment on darwinian evolution. In particular it simulates macrophages in a sea of bacteria inspired by the work done in [Computer recreations](https://www.scientificamerican.com/article/computer-recreations-1989-05/) by A. K. Dewdney.

The simulation shows the behavior of evolution given a system where macrophages move according to a set of numbers that codify the probability of moving in a certain direction (Genes). A macrophage has the hability to eat (gaining energy), move (consuming energy) and reproduce (consuming a lot of energy). When a macrophage's internal energy reaches zero, it dies. When it reproduces, the child macrophage will inherit energy from the parent as well as their Genes, but with a random modificacion to the values (simulating genetic mutation). 

The simulation is finished and works as expected. It has been shown to replicate behaviors like r-k selection when changing gestation time, as wells a predator-prey population oscilation to name a few. 

## Use
To run the simulation, in the working directory you should have:
- Simulacion_Bichines.cpp &rightarrow; where the simulation runs
- Random64.h &rightarrow; The random number generator we used
- (optional) 2 empty folders named "Nodos" and "Edges" &rightarrow; Where nodes and edges data will be stored if the folders exist


The food distribution in the simulation depends on one parameter:
```cpp
int food_dis;
```
it uses 0 as default and is changed as the first command line argument.
| food_dis | Distribution  |
| -------- | --------------| 
| 0        | uniform       |
| 1        | gaussian      | 
| 2        | Garden of eden| 

if running several scenarios, you may need to change the name of the files everytime you want to try a new scenario given the files will be rewriten.

`TMAX` is the amount of time iterations the simulation will run for. 10 000 by default and is changed with the second command line argument.

`rand_seed` is the seed for the random generator and is changed with the third command line argument, it has to be a integer.



Compile with any c++ compiler and run as any other program. The executable outputs several files 


If you wanted to see the evolution process like this.

![Image of the simulation; Garden of eden](Resultados/Imagenes_readme/Jardin_eden.png "Map of the simulation")

You will need to open the .gp file created by the executable with [Gnuplot](http://www.gnuplot.info).


## Run Time and population
As with a real ecosystem there is a carrying capacity


## CANTIDAD DE BICHINES Y TIEMPO DE COMPUTO

El tiempo de simulación para una simulación grande es bastante costo. En Su configuración actual Nfood = 20000,Biome_energy=50000, el sistema tiene una capacidad de carga de unos 1700Bichines. Esto es bastante intensivo computacionalmente, unas 2h de procesador por t10000; Dada la secuencialidad de los ciclos for implementados, puede tomar mucho tiempo pero no mucho poder de procesamiento. Es decir, se puede correr varios escenarios en el mismo procesador simultáneamente sin problema. 

Si Biome_energy<50000 La población no sobrevive el colapso del jardín y se extinguen. Esto es especialmente importante para el caso de la distribución Eden (dis=2), las distribuciones con fronteras suaves como gaussianas dobles o distribuciones uniformes no sufren de un colapso poblacional de la misma manera que lo hace la distribución del Eden. Por consecuencia pueden usar un Biome_energy>30000. 


Para las distribuciones de Geometrías más complejas se espera que también sea necesaria una capacidad de carga alta del ecosistema. Para distribuciones como Eden2 se espera necesitar un Biome_energy aún más alto del probado 50000. 

