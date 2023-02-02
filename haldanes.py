import numpy as np
import math as m
import matplotlib.pylab as plt
data=np.loadtxt('Haldanes.txt')


sd=[]
prom=[]
time=data[:,0]
pop=data[:,1]
for gene in range(8):
    prom.append(data[:,2*gene+2])
    sd.append(data[:,2*gene+3])

proms=np.array(prom)
Sd=np.array(sd)

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

Hal_data=Haldanes(time,proms,Sd,1)
print(Hal_data)
time_h=time[1:]
plt.plot(time_h,Haldanes(time,proms,Sd,0))
plt.show()

'''
for gene in range(8):
    plt.plot(time,proms[gene],label=str(gene)+' gen')
plt.legend()
plt.show()
'''