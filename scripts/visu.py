import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("../particule.xyz", usecols=(1, 2, 3), skiprows=1)

fig = plt.figure()
ax = fig.add_subplot(projection="3d")

ax.scatter(*data.T)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()