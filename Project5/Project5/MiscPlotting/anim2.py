import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

data = np.zeros([11,11])

fig = plt.figure()
ax = plt.axes(xlim=(0,1),ylim=(0,1),xlabel="Position x",ylabel="Position y")

#line, = ax.contourf([],[],[],[])

t = 0

x = np.linspace(0,1,11); y = np.linspace(0,1,11)
X,Y = np.meshgrid(x,y)

iterate = 0

dt = 1/400

file = open("Explicit2D.txt","r")
"""
def init():
    line.set_data([],[],[[],[]])
    return line,
"""
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
    file.readline()

    plt.contourf(X,Y,data,60)
    plt.title("T = %.3f"%(i*dt))
    return line

anim = animation.FuncAnimation(fig,animate,frames=40001,interval=10)
anim.save("expdelxy10.mp4")
plt.show()


"""
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

plt.contourf(x,y,data[-1])
plt.show()
"""
