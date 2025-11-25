import numpy as np

data = np.loadtxt("../particule.xyz", usecols=(1, 2, 3), skiprows=1)
with open("../particule.data", "wb") as f:
    data.tofile(f, "")
