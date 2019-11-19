import numpy as np
import matplotlib.pyplot as plt

E = []
Mabs = []

file = open("test.txt","r")

file.readline()

for line in file:
    words = line.split()
    E.append(float(words[0]))
    Mabs.append(float(words[-1]))

file.close()

file = open("test.txt","r")

E2 = []
Mabs2 = []

file.readline()
for line in file:
    words = line.split()
    E2.append(float(words[0]))
    Mabs2.append(float(words[-1]))


file.close()

t = np.linspace(0,len(E)-1,len(E))

plt.figure()
plt.plot(t,E,label="Ordered Initiation")
plt.plot(t,E2,label="Random Initiation")
plt.legend(fontsize="large")
plt.title("Energy per spin. T = 2.4",fontsize="x-large")
plt.xlabel("Number of MC Sweeps",fontsize="x-large")
plt.ylabel("Mean energy per spin [<E>/ J ]",fontsize="x-large")
#plt.xlim(-2000,600000)
plt.ylim(-1.24,-1.225)
plt.grid()

plt.figure()
plt.plot(t,Mabs,label="Ordered Initiation")
plt.plot(t,Mabs2,label="Random Initiation")
plt.legend(fontsize="large")
plt.title("Magnetization per spin. T = 2.4",fontsize="x-large")
plt.xlabel("Number of MC Sweeps",fontsize="x-large")
plt.ylabel("Mean Magnetization per spin [<|M|>/ J ]",fontsize="x-large")
#plt.xlim(-2000,600000)
plt.ylim(0.44,0.46)
plt.grid()

plt.figure()
plt.hist(E)
plt.show()

plt.show()
