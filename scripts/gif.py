import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import os

files = os.listdir("../saved/")
N = 1000

pos = np.zeros((len(files), N, 3))
ns = []
for filename in files:
    ns.append(int(filename.split(".")[0]))

ns.sort()

for i, n in enumerate(ns):
    pos[i] = np.fromfile(f"../saved/{n}.data").reshape((N, 3))

ims = []
fig = plt.figure()
fig.tight_layout()
ax = fig.add_subplot(projection="3d")

def animate(i):
    i = i
    ax.clear()
    # ax.set_xlim(-15, 15)
    # ax.set_ylim(-15, 15)
    # ax.set_zlim(-15, 15)
    ax.set_title(i)
    d = ax.scatter(*pos[i].T, marker=".", color="blue")
    return d

ani = animation.FuncAnimation(
    fig, animate, len(files)-1, interval=100
)
plt.show()
# ani.save("haha.gif")
