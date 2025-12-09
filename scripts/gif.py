import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import os

files = os.listdir("../saved/")
N = 1000
L = 34

pos = np.zeros((len(files), N, 3))
ns = []
for filename in files:
    ns.append(int(filename.split(".")[0]))

ns.sort()

for i, n in enumerate(ns):
    pos[i] = np.fromfile(f"../saved/{n}.data").reshape((N, 3))

pos_loc = np.fmod(pos, L)
pos_loc += L 
pos_loc = np.fmod(pos_loc, L)

ims = []
fig = plt.figure()
ax1 = fig.add_subplot(1, 2, 1, projection="3d")
ax2 = fig.add_subplot(1, 2, 2, projection="3d")
fig.set_size_inches(10, 6)
fig.tight_layout()


def animate(i):
    i = i
    ax1.clear()
    ax1.scatter(*pos[i].T, marker=".", color="blue")
    ax2.clear()
    ax2.set_title(i)
    d = ax2.scatter(*pos_loc[i].T, marker=".", color="blue")
    return d


ani = animation.FuncAnimation(
    fig, animate, len(files)-1, interval=100
)
# plt.show()
ani.save("haha.gif")
