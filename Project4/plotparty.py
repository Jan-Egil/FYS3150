import numpy as np
import matplotlib.pyplot as plt
"""
Initiate data collection lists
"""
E1 = []; M1 = []; C1 = []; X1 = []
E2 = []; M2 = []; C2 = []; X2 = []
E3 = []; M3 = []; C3 = []; X3 = []
E4 = []; M4 = []; C4 = []; X4 = []
T = [2,2.05,2.1,2.15,2.2,2.25,2.3]

"""
Read files
"""

file = open("40X40_2.txt","r")

file.readline()

for line in file:
    words = line.split()
    E1.append(float(words[1]))
    M1.append(float(words[2]))
    C1.append(float(words[3]))
    X1.append(float(words[4]))

file.close()

file = open("60X60_2.txt","r")

file.readline()

for line in file:
    words = line.split()
    E2.append(float(words[1]))
    M2.append(float(words[2]))
    C2.append(float(words[3]))
    X2.append(float(words[4]))

file.close()

file = open("80X80_2.txt","r")

file.readline()

for line in file:
    words = line.split()
    E3.append(float(words[1]))
    M3.append(float(words[2]))
    C3.append(float(words[3]))
    X3.append(float(words[4]))

file.close()

file = open("100X100_2.txt","r")

file.readline()

for line in file:
    words = line.split()
    E4.append(float(words[1]))
    M4.append(float(words[2]))
    C4.append(float(words[3]))
    X4.append(float(words[4]))

file.close()

"""
Plot Energy
"""
plt.figure()
plt.title("Energy",fontsize="x-large")
plt.plot(T,E1,"o--",label="40X40-lattice")
plt.plot(T,E2,"o--",label="60X60-lattice")
plt.plot(T,E3,"o--",label="80X80-lattice")
plt.plot(T,E4,"o--",label="100X100-lattice")
plt.xlabel("Temperature",fontsize="x-large")
plt.ylabel("$\langle E \\rangle/N^2$",fontsize="x-large")
plt.legend(fontsize="large")
plt.grid()


"""
Plot Magnetizaton
"""
plt.figure()
plt.title("Magnetization",fontsize="x-large")
plt.plot(T,M1,"o--",label="40X40-lattice")
plt.plot(T,M2,"o--",label="60X60-lattice")
plt.plot(T,M3,"o--",label="80X80-lattice")
plt.plot(T,M4,"o--",label="100X100-lattice")
plt.xlabel("Temperature",fontsize="x-large")
plt.ylabel("$\langle |M| \\rangle/N^2$",fontsize="x-large")
plt.legend(fontsize="large")
plt.grid()


"""
Plot Heat Capacity
"""
plt.figure()
plt.title("Heat Capacity",fontsize="x-large")
plt.plot(T,C1,"o--",label="40X40-lattice")
plt.plot(T,C2,"o--",label="60X60-lattice")
plt.plot(T,C3,"o--",label="80X80-lattice")
plt.plot(T,C4,"o--",label="100X100-lattice")
plt.xlabel("Temperature",fontsize="x-large")
plt.ylabel("$C_v/N^2$",fontsize="x-large")
plt.legend(fontsize="large")
plt.grid()


"""
Plot Susceptibility
"""
plt.figure()
plt.title("Susceptibility",fontsize="x-large")
plt.plot(T,X1,"o--",label="40X40-lattice")
plt.plot(T,X2,"o--",label="60X60-lattice")
plt.plot(T,X3,"o--",label="80X80-lattice")
plt.plot(T,X4,"o--",label="100X100-lattice")
plt.xlabel("Temperature",fontsize="x-large")
plt.ylabel("$\chi/N^2$",fontsize="x-large")
plt.legend(fontsize="large")
plt.grid()

plt.show()
