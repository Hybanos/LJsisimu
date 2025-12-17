import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import os

data = np.fromfile(f"../out.data")
N = int(data[0])
L = int(data[1])

# first 2 doubles are particule count and image count
iterations = data[2:]
# 7 doubles per iteration + position array
n_iterations = len(data[2:]) // (7 + N * 3)
iterations = iterations.reshape((n_iterations, 7 + N * 3))

iters = iterations.T[0]
U = iterations.T[1]
T = iterations.T[2]
E_k = iterations.T[3]
center_of_mass = iterations.T[4:7].T.reshape((n_iterations, 3))
pos = iterations.T[7:3008].T.reshape((n_iterations, N, 3))

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
    print(f"{round(i/n_iterations * 100, 2)}%   ", end="\r")
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
    ax3.plot(iters[:i], U[:i], color="tab:blue", label="Potential E.")
    ax3.plot(iters[:i], E_k[:i], color="tab:green", label="Kinetic E.")
    ax3.plot(iters[:i], E_k[:i] + U[:i], color="tab:brown", label="Total E.")
    ax3.legend()
    ax4.clear()
    ax4.plot(iters[:i], T[:i], color="tab:red", label="Temp")
    ax4.legend()

ani = animation.FuncAnimation(
    fig, animate, n_iterations-1, interval=100
)
# plt.show()
ani.save("haha.gif")
print("done          ")
