import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import os

files = os.listdir("../saved/")
N = 0
L = 0

ns = []
for filename in files:
    ns.append(int(filename.split(".")[0]))
ns.sort()

tmp = np.fromfile(f"../saved/{ns[0]}.data")
N = int(tmp[0])
L = int(tmp[1])

iters = np.zeros((len(files)))
U = np.zeros((len(files)))
T = np.zeros((len(files)))
E_k = np.zeros((len(files)))
center_of_mass = np.zeros((len(files), 3))
pos = np.zeros((len(files), N, 3))

for i, n in enumerate(ns):
    tmp = np.fromfile(f"../saved/{n}.data")
    iters[i] = tmp[2]
    U[i] = tmp[3]
    T[i] = tmp[4]
    E_k[i] = tmp[5]
    center_of_mass[i] = tmp[6:9]
    pos[i] =tmp[9:].reshape((N, 3))

pos_loc = np.fmod(pos, L)
pos_loc += L 
pos_loc = np.fmod(pos_loc, L)

ims = []
fig = plt.figure()
ax1 = fig.add_subplot(2, 2, 1, projection="3d")
ax2 = fig.add_subplot(2, 2, 2, projection="3d")
ax3 = fig.add_subplot(2, 2, 3)
ax4 = fig.add_subplot(2, 2, 4)
fig.set_size_inches(9, 9)
fig.tight_layout()


def animate(i):
    print(f"{round(i/len(files) * 100, 2)}%   ", end="\r")
    fig.suptitle(f"Iteration {int(iters[i])}")
    ax1.clear()
    ax1.scatter(*pos[i].T, marker=".", color="blue")
    ax1.scatter(*center_of_mass[i], color="tab:red")
    ax1.set_title("Global pos")
    ax1.set_xlim(np.min(pos[-1]), np.max(pos[-1]))
    ax1.set_ylim(np.min(pos[-1]), np.max(pos[-1]))
    ax1.set_zlim(np.min(pos[-1]), np.max(pos[-1]))
    ax2.clear()
    ax2.scatter(*pos_loc[i].T, marker=".", color="blue")
    ax2.scatter(*center_of_mass[i], color="tab:red")
    ax2.set_title("local pos")
    ax3.clear()
    ax3.plot(iters[1:i], U[1:i], color="tab:blue", label="Potential E.")
    ax3.plot(iters[:i], E_k[:i], color="tab:green", label="Kinetic E.")
    ax3.plot(iters[1:i], E_k[1:i] + U[1:i], color="tab:brown", label="Total E.")
    ax3.legend()
    ax4.clear()
    ax4.plot(iters[:i], T[:i], color="tab:red", label="Temp")
    ax4.legend()


ani = animation.FuncAnimation(
    fig, animate, len(files)-1, interval=100
)
# plt.show()
ani.save("haha.gif")
print("done          ")
