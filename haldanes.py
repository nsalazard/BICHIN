import numpy as np
import matplotlib.pylab as plt
data=np.loadtxt('Haldanes.txt')

for i in range(1,8):
    plt.plot(data[:,0],data[:,i],label=str(i-1))
plt.legend()
plt.show()
print(data[0])
