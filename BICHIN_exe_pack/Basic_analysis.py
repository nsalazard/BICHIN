import numpy as np
import math as m
import matplotlib.pylab as plt

output_folder='Basic_results/'

data=np.loadtxt('Haldanes.txt')


standart_d=[] #standard deviation
prom=[]
time=data[:,0]
pop=data[:,1]
for gene in range(8):
	prom.append(data[:,2*gene+2])
	standart_d.append(data[:,2*gene+3])

proms=np.array(prom)
Sd=np.array(standart_d)

def Haldanes(time,proms,sd,gene):
	Sd_pooled=[]
	for t in range(len(time)-1):
		q=((pop[t]-1)*Sd[gene,t]-(pop[t+1]-1)*Sd[gene,t+1])/(pop[t]+pop[t+1])
		Sd_pooled.append(q)
	H=[]
	for t in range(len(time)-1):
		q2=(m.log(proms[gene,t+1])-m.log(proms[gene,t]))/(Sd_pooled[t]*(time[t+1]-time[t]))
		H.append(q2)
	
	return H





for gene in range(8):
	plt.plot(time,proms[gene],label=str(gene)+' gene')
	
plt.xlabel('Time (clicks)')
plt.ylabel('Proportion of the gene in the genoma (%)')
plt.title('Evolution of the genomic distribution over time ')
plt.legend()
plt.savefig(output_folder+'promedios_genes.png',bbox_inches='tight')
plt.clf()
plt.close()

for gene in range(8):
	plt.plot(time,Sd[gene],label=str(gene)+' gene')
plt.xlabel('Time (clicks)')
plt.ylabel('Std deviation of the genes (u)')
plt.title('Evolution of the genomic distribution over time ')
plt.legend()
plt.savefig(output_folder+'Stv_genes.png',bbox_inches='tight')
plt.clf()
plt.close()

Hal_data=Haldanes(time,proms,Sd,1)

time_h=time[1:]
plt.plot(time_h,Haldanes(time,proms,Sd,0))
plt.xlabel('Time (clicks)')
plt.ylabel('Computed Haldanes (u)')
plt.savefig(output_folder+'tasa_evolutiva.png',bbox_inches='tight')
plt.close()


population=np.loadtxt('population.txt')[:,1]

plt.plot(np.arange(len(population)),population.T)
plt.xlabel('Time (clicks)')
plt.ylabel('Population')
plt.savefig(output_folder+'population2.png',bbox_inches='tight')
plt.clf()
plt.close()
