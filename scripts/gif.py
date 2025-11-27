import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import os

files = os.listdir("../saved/")

pos = np.zeros((len(files), 1000, 3))

for filename in files:
    n = int(filename.split(".")[0])
    pos[n] = np.fromfile("../saved/" + filename).reshape((1000, 3))

ims = []
fig = plt.figure()
ax = fig.add_subplot(projection="3d")

def animate(i):
    d = ax.scatter(*pos[i].T, color="blue")
    return d

ani = animation.FuncAnimation(
    fig, animate, interval=20
)
plt.show()

