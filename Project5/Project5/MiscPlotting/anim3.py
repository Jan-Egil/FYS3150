import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

data = np.zeros([351,121])

fig = plt.figure()
ax = plt.axes(xlim=(0,1),ylim=(0,1),xlabel="Position x",ylabel="Position y")


#line, = ax.contourf([],[],[],[])

t = 0

x = np.linspace(0,1,351); y = np.linspace(0,1,121)
X,Y = np.meshgrid(y,x)

iterate = 0

dt = 1/10000

file = open("PGPExplicit2D.txt","r")

def animate(i):
    iterate = 0
    for line in file:
        elems = line.split()
        if elems == []:
            iterate = 0
            break
        else:
            data[iterate] = elems
            iterate += 1

    contourplot = plt.contourf(X,Y,data,50)
    plt.title("T = %.4fgY"%(i*dt))
    return line

anim = animation.FuncAnimation(fig,animate,frames=9999,interval=10)
#anim.save("PGPHalfLife.mp4")
plt.show()
