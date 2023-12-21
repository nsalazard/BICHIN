<div align="center">

# BICHIN
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



Compile with any c++ compiler and you should get an executable. 
If you wanted to see the evolution process like this.

![Example Image](Resultados/Imagenes_readme/Jardin_eden.png "This is an example image")

## USO


En tanto las Simulacion_Bichines.cpp, ya está terminada y funcional. Necesitamos correrla en diferentes escenarios y extraer los datos para análisis estadístico, el entorno de ejecución necesita Random64.h, 2 carpetas, Nodos, Edges. 

Para probar diferentes escenarios solo sería necesario cambiar la distribución de la comida. Esto tiene que suceder en 2 lugares, en Selection>food_distribution() y Food>ReStart(). Estas 2 funciones se encargan de la distribución inicial y redistribución de la comida en el mapa respectivamente y ambas dependen del parámetro. 

```php
int food_dis=I;
```

Por el momento existen 3 distribuciones implementadas para I=0; uniforme, I=1; gaussiana, I=2; Eden. 

Creen nuevas instrucciones de distribución para diferentes valores de I, para que se puedan integrar todos al final en un solo código. Corran sus escenarios y en la carpeta de resultados abran una carpeta con el nombre de su escenario y pongan ahí sus resultados.


## CANTIDAD DE BICHINES Y TIEMPO DE COMPUTO

El tiempo de simulación para una simulación grande es bastante costo. En Su configuración actual Nfood = 20000,Biome_energy=50000, el sistema tiene una capacidad de carga de unos 1700Bichines. Esto es bastante intensivo computacionalmente, unas 2h de procesador por t10000; Dada la secuencialidad de los ciclos for implementados, puede tomar mucho tiempo pero no mucho poder de procesamiento. Es decir, se puede correr varios escenarios en el mismo procesador simultáneamente sin problema. 

Si Biome_energy<50000 La población no sobrevive el colapso del jardín y se extinguen. Esto es especialmente importante para el caso de la distribución Eden (dis=2), las distribuciones con fronteras suaves como gaussianas dobles o distribuciones uniformes no sufren de un colapso poblacional de la misma manera que lo hace la distribución del Eden. Por consecuencia pueden usar un Biome_energy>30000. 


Para las distribuciones de Geometrías más complejas se espera que también sea necesaria una capacidad de carga alta del ecosistema. Para distribuciones como Eden2 se espera necesitar un Biome_energy aún más alto del probado 50000. 

