import numpy as np
import matplotlib.pyplot as plt

file = open("histogramdata2.txt","r")

distribution = []

for line in file:
    words = line.split()
    distribution.append(float(words[0]))

file.close()

file = open("histogramdata1.txt","r")

distribution2 = []

for line in file:
    words = line.split()
    distribution2.append(float(words[0]))

file.close()

#247863

plt.figure()
plt.hist(distribution,range=(-800,800),bins=400)
plt.title("Histogram for temperature T = 2.4",fontsize="x-large")
plt.xlabel("Energy [1/J]",fontsize="x-large")
plt.ylabel("Counts",fontsize="x-large")

plt.figure()
plt.hist(distribution2,range=(-800,800),bins=400)
plt.title("Histogram for temperature T = 1.0",fontsize="x-large")
plt.xlabel("Energy [1/J]",fontsize="x-large")
plt.ylabel("Counts",fontsize="x-large")

plt.show()

print(distribution[-1])
print(distribution2[-1])
