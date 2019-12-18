import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

DataExp = np.loadtxt("Explicit1DDelx10.txt")
DataExp2 = np.loadtxt("Explicit1DDelx100.txt")
DataImp = np.loadtxt("Implicit1DDelx10.txt")
DataImp2 = np.loadtxt("Implicit1DDelx100.txt")
DataCN = np.loadtxt("CrankNicolson1DDelx10.txt")
DataCN2 = np.loadtxt("CrankNicolson1DDelx100.txt")

"""
fig = plt.figure()
ax = plt.axes(xlim=(-0.2,1.2),ylim=(-0.2,1.1),xlabel="Position",ylabel="$u(x,t)$")
plt.grid()

line, = ax.plot([],[])

def init():
    line.set_data([],[])
    return line,
def animate(i):
    x = np.linspace(0,1,11)
    y = DataExp[i,1:]
    line.set_data(x,y)
    plt.title("T = %.3f"%DataExp[i,0])
    return line,

anim = animation.FuncAnimation(fig, animate, init_func=init,frames=199,interval=20)
anim.save('expdelx10.mp4')
plt.show()
"""

"""
fig2 = plt.figure()
ax2 = plt.axes(xlim=(-0.2,1.2),ylim=(-0.2,1.1),xlabel="Position",ylabel="$u(x,t)$")
plt.grid()

line2, = ax2.plot([],[])

def init2():
    line2.set_data([],[])
    return line2,
def animate2(i):
    x2 = np.linspace(0,1,101)
    y2 = DataExp2[i*100,1:]
    line2.set_data(x2,y2)
    plt.title("T = %.3f"%DataExp2[i*100,0])
    return line2,

anim2 = animation.FuncAnimation(fig2, animate2, init_func=init2,frames=199,interval=30)
anim2.save('expdelx100.mp4')
plt.show()
"""

"""
fig = plt.figure()
plt.grid()
ax = plt.axes(xlim=(-0.2,1.2),ylim=(-0.2,1.1),xlabel="Position")

line, = ax.plot([],[])

def init():
    line.set_data([],[])
    return line,
def animate(i):
    x = np.linspace(0,1,11)
    y = DataImp[i,1:]
    line.set_data(x,y)
    plt.title("T = %.3f"%DataImp[i,0])
    return line,

anim = animation.FuncAnimation(fig, animate, init_func=init,frames=199,interval=10)
plt.grid()
#anim.save('impdelx10.mp4')
plt.show()
"""

"""
fig = plt.figure()
plt.grid()
ax = plt.axes(xlim=(-0.2,1.2),ylim=(-0.2,1.1),xlabel="Position")

line, = ax.plot([],[])

def init():
    line.set_data([],[])
    return line,
def animate(i):
    x = np.linspace(0,1,101)
    y = DataImp2[i*100,1:]
    line.set_data(x,y)
    plt.title("T = %.3f"%DataImp2[i*100,0])
    return line,

anim = animation.FuncAnimation(fig, animate, init_func=init,frames=199,interval=10)
plt.grid()
#anim.save('impdelx100.mp4')
plt.show()
"""


fig = plt.figure()
plt.grid()
ax = plt.axes(xlim=(-0.2,1.2),ylim=(-0.2,1.1),xlabel="Position")

line, = ax.plot([],[])

def init():
    line.set_data([],[])
    return line,
def animate(i):
    x = np.linspace(0,1,11)
    y = DataCN[i,1:]
    line.set_data(x,y)
    plt.title("T = %.3f"%DataCN[i,0])
    return line,

anim = animation.FuncAnimation(fig, animate, init_func=init,frames=199,interval=10)
plt.grid()
anim.save('CNdelx10.mp4')
plt.show()


"""
fig = plt.figure()
plt.grid()
ax = plt.axes(xlim=(-0.2,1.2),ylim=(-0.2,1.1),xlabel="Position")

line, = ax.plot([],[])

def init():
    line.set_data([],[])
    return line,
def animate(i):
    x = np.linspace(0,1,101)
    y = DataCN2[i*100,1:]
    line.set_data(x,y)
    plt.title("T = %.3f"%DataCN2[i*100,0])
    return line,

anim = animation.FuncAnimation(fig, animate, init_func=init,frames=199,interval=10)
plt.grid()
#anim.save('CNdelx100.mp4')
plt.show()
"""
