import numpy as np
import matplotlib.pyplot as plt

DataExp = np.loadtxt("Explicit1Derror.txt")
DataImp = np.loadtxt("Implicit1Derror.txt")
DataCN = np.loadtxt("CrankNicolson1Derror.txt")

t1 = DataExp[:,0]
actualdata1 = DataExp[:,1]
t2 = DataImp[:,0]
actualdata2 = DataImp[:,1]
t3 = DataCN[:,0]
actualdata3 = DataCN[:,1]

plt.figure()
plt.plot(t1,actualdata1,label="Explicit method")
plt.plot(t2,actualdata2,label="Implicit method")
plt.plot(t3,actualdata3,label="Crank-Nicolson")
plt.semilogy()
plt.legend(fontsize="x-large")
plt.xlabel("Time (Rel Units)",fontsize="x-large")
plt.ylabel("$\log{\epsilon(t)}$",fontsize="x-large")
plt.grid()
plt.title("Finite Difference Error $\Delta x = 0.01$",fontsize="x-large")

DataExp = np.loadtxt("Explicit1Derror2.txt")
DataImp = np.loadtxt("Implicit1Derror2.txt")
DataCN = np.loadtxt("CrankNicolson1Derror2.txt")

t1 = DataExp[:,0]
actualdata1 = DataExp[:,1]
t2 = DataImp[:,0]
actualdata2 = DataImp[:,1]
t3 = DataCN[:,0]
actualdata3 = DataCN[:,1]

plt.figure()
plt.plot(t1,actualdata1,label="Explicit method")
plt.plot(t2,actualdata2,label="Implicit method")
plt.plot(t3,actualdata3,label="Crank-Nicolson")
plt.semilogy()
plt.legend(fontsize="x-large")
plt.xlabel("Time (Rel Units)",fontsize="x-large")
plt.ylabel("$\log{\epsilon(t)}$",fontsize="x-large")
plt.title("Finite Difference Error $\Delta x = 0.1$",fontsize="x-large")
plt.grid()
plt.show()
