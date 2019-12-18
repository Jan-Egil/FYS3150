import numpy as np
import matplotlib.pyplot as plt

DataExp = np.loadtxt("Explicit1DDelx10.txt")
DataExp2 = np.loadtxt("Explicit1DDelx100.txt")
DataImp = np.loadtxt("Implicit1DDelx10.txt")
DataImp2 = np.loadtxt("Implicit1DDelx100.txt")
DataCN = np.loadtxt("CrankNicolson1DDelx10.txt")
DataCN2 = np.loadtxt("CrankNicolson1DDelx100.txt")

t = DataExp[:,0]
x1 = np.linspace(0,1,11)
X,T = np.meshgrid(x1,t)
actualdata = DataExp[:,1:]

plt.figure()
plt.contourf(T,X,actualdata,100)
plt.colorbar(label="Value of function u")
plt.xlabel("Time t (Rel units)",fontsize="x-large")
plt.ylabel("Position x (Rel units)",fontsize="x-large")
plt.title("$u(x,t)$ visualized [$u(0,t) = 0$ and $u(1,t) = 1$] [$\Delta x = \\frac{1}{10}$]",fontsize="x-large")

t = DataExp2[:,0]
x1 = np.linspace(0,1,101)
X,T = np.meshgrid(x1,t)
actualdata = DataExp2[:,1:]

plt.figure()
plt.contourf(T,X,actualdata,100)
plt.colorbar(label="Value of function u")
plt.xlabel("Time t (Rel units)",fontsize="x-large")
plt.ylabel("Position x (Rel units)",fontsize="x-large")
plt.title("$u(x,t)$ visualized [$u(0,t) = 0$ and $u(1,t) = 1$] [$\Delta x = \\frac{1}{100}$]",fontsize="x-large")

t = DataImp[:,0]
x1 = np.linspace(0,1,11)
X,T = np.meshgrid(x1,t)
actualdata = DataImp[:,1:]

plt.figure()
plt.contourf(T,X,actualdata,100)
plt.colorbar(label="Value of function u")
plt.xlabel("Time t (Rel units)",fontsize="x-large")
plt.ylabel("Position x (Rel units)",fontsize="x-large")
plt.title("$u(x,t)$ visualized [$u(0,t) = 0$ and $u(1,t) = 1$] [$\Delta x = \\frac{1}{10}$]",fontsize="x-large")

t = DataImp2[:,0]
x1 = np.linspace(0,1,101)
X,T = np.meshgrid(x1,t)
actualdata = DataImp2[:,1:]

plt.figure()
plt.contourf(T,X,actualdata,100)
plt.colorbar(label="Value of function u")
plt.xlabel("Time t (Rel units)",fontsize="x-large")
plt.ylabel("Position x (Rel units)",fontsize="x-large")
plt.title("$u(x,t)$ visualized [$u(0,t) = 0$ and $u(1,t) = 1$] [$\Delta x = \\frac{1}{100}$]",fontsize="x-large")

t = DataCN[:,0]
x1 = np.linspace(0,1,11)
X,T = np.meshgrid(x1,t)
actualdata = DataCN[:,1:]

plt.figure()
plt.contourf(T,X,actualdata,100)
plt.colorbar(label="Value of function u")
plt.xlabel("Time t (Rel units)",fontsize="x-large")
plt.ylabel("Position x (Rel units)",fontsize="x-large")
plt.title("$u(x,t)$ visualized [$u(0,t) = 0$ and $u(1,t) = 1$] [$\Delta x = \\frac{1}{10}$]",fontsize="x-large")

t = DataCN2[:,0]
x1 = np.linspace(0,1,101)
X,T = np.meshgrid(x1,t)
actualdata = DataCN2[:,1:]

plt.figure()
plt.contourf(T,X,actualdata,100)
plt.colorbar(label="Value of function u")
plt.xlabel("Time t (Rel units)",fontsize="x-large")
plt.ylabel("Position x (Rel units)",fontsize="x-large")
plt.title("$u(x,t)$ visualized [$u(0,t) = 0$ and $u(1,t) = 1$] [$\Delta x = \\frac{1}{100}$]",fontsize="x-large")


plt.show()
