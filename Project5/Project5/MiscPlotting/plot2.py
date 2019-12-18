import numpy as np
import matplotlib.pyplot as plt

"""
data = np.zeros([40001,11,11])

file = open("Explicit2D.txt","r")

t = 0; x = 0

for line in file:
    elems = line.split()
    if elems == []:
        t += 0.5
        x = 0
    else:
        data[int(t),x] = elems
        x += 1

x = np.linspace(0,1,11); y = np.linspace(0,1,11)

X,Y = np.meshgrid(x,y)

plt.figure()
plt.contourf(X,Y,data[-1],100)
a = plt.colorbar()
a.set_label(label="Value of function $u(x,y,t\\rightarrow\infty$)",fontsize="x-large")
plt.xlabel("Position x (Rel units)",fontsize="x-large")
plt.ylabel("Position y (Rel units)",fontsize="x-large")
plt.title("$u(x,y,t\\rightarrow\infty)$ with boundaries: \n $u(0,y,t) = u(1,y,t) = 1$, \n $u(x,0,t) = u(x,1,t) = 0$.",fontsize="x-large")

plt.figure()
plt.contourf(X,Y,data[20],100)
b = plt.colorbar()
b.set_label(label="Value of function $u(x,y,t=20\cdot\Delta t)$",fontsize="x-large")
plt.xlabel("Position x (Rel units)",fontsize="x-large")
plt.ylabel("Position y (Rel units)",fontsize="x-large")
plt.title("$u(x,y,t=20\cdot\Delta t)$ with boundaries: \n $u(0,y,t) = u(1,y,t) = 1$, \n $u(x,0,t) = u(x,1,t) = 0$.",fontsize="x-large")


plt.figure()
plt.contourf(X,Y,data[3],100)
c = plt.colorbar()
c.set_label(label="Value of function $u(x,y,t=3\cdot\Delta t)$",fontsize="x-large")
plt.xlabel("Position x (Rel units)",fontsize="x-large")
plt.ylabel("Position y (Rel units)",fontsize="x-large")
plt.title("$u(x,y,t=3\cdot\Delta t)$ with boundaries: \n $u(0,y,t) = u(1,y,t) = 1$, \n $u(x,0,t) = u(x,1,t) = 0$.",fontsize="x-large")

plt.figure()
plt.contourf(X,Y,data[0],100)
a = plt.colorbar()
a.set_label(label="Value of function $u(x,y,t=0)$",fontsize="x-large")
plt.xlabel("Position x (Rel units)",fontsize="x-large")
plt.ylabel("Position y (Rel units)",fontsize="x-large")
plt.title("$u(x,y,t=0)$ with boundaries: \n $u(0,y,t) = u(1,y,t) = 1$, \n $u(x,0,t) = u(x,1,t) = 0$.",fontsize="x-large")

plt.show()
"""

data2 = np.zeros([2000,101,101])

file2 = open("Explicit2Dh100.txt","r")

t = 0; x = 0

for line in file2:
    if t >= 2000:
        break

    elems = line.split()
    if elems == []:
        t += 0.5
        x = 0
    else:
        data2[int(t),x] = elems
        x += 1

x = np.linspace(0,1,101); y = np.linspace(0,1,101)

X,Y = np.meshgrid(x,y)

plt.figure()
plt.contourf(X,Y,data2[-1],100)
a = plt.colorbar()
a.set_label(label="Value of function $u(x,y,t\\rightarrow\infty$)",fontsize="x-large")
plt.xlabel("Position x (Rel units)",fontsize="x-large")
plt.ylabel("Position y (Rel units)",fontsize="x-large")
plt.title("$u(x,y,t\\rightarrow\infty)$ with boundaries: \n $u(0,y,t) = u(1,y,t) = 0$, \n $u(x,0,t) = u(x,1,t) = 1$.",fontsize="x-large")

plt.figure()
plt.contourf(X,Y,data2[1500],100)
b = plt.colorbar()
b.set_label(label="Value of function $u(x,y,t=200\cdot\Delta t)$",fontsize="x-large")
plt.xlabel("Position x (Rel units)",fontsize="x-large")
plt.ylabel("Position y (Rel units)",fontsize="x-large")
plt.title("$u(x,y,t=200\cdot\Delta t)$ with boundaries: \n $u(0,y,t) = u(1,y,t) = 0$, \n $u(x,0,t) = u(x,1,t) = 1$.",fontsize="x-large")


plt.figure()
plt.contourf(X,Y,data2[30],100)
c = plt.colorbar()
c.set_label(label="Value of function $u(x,y,t=30\cdot\Delta t)$",fontsize="x-large")
plt.xlabel("Position x (Rel units)",fontsize="x-large")
plt.ylabel("Position y (Rel units)",fontsize="x-large")
plt.title("$u(x,y,t=30\cdot\Delta t)$ with boundaries: \n $u(0,y,t) = u(1,y,t) = 0$, \n $u(x,0,t) = u(x,1,t) = 1$.",fontsize="x-large")

plt.figure()
plt.contourf(X,Y,data2[0],100)
a = plt.colorbar()
a.set_label(label="Value of function $u(x,y,t=0)$",fontsize="x-large")
plt.xlabel("Position x (Rel units)",fontsize="x-large")
plt.ylabel("Position y (Rel units)",fontsize="x-large")
plt.title("$u(x,y,t=0)$ with boundaries: \n $u(0,y,t) = u(1,y,t) = 0$, \n $u(x,0,t) = u(x,1,t) = 1$.",fontsize="x-large")

plt.show()
