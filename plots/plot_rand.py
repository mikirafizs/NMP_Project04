import numpy as np
from matplotlib import pyplot as plt

plt.figure()
d = np.loadtxt(f"data/random_vec_sol.csv", delimiter=",")
plt.plot(d[:, 0], d[:, 2])
plt.xscale("log")
plt.yscale("log")
plt.xlabel("iteration")
plt.ylabel("||r||_2")
plt.legend()
plt.grid(True)
plt.savefig("plots/norm_r_random.png", dpi=200)

plt.figure()
d = np.loadtxt(f"data/random_vec_sol.csv", delimiter=",")
plt.plot(d[:, 0], d[:, 1])
plt.xscale("log")
plt.yscale("log")
plt.xlabel("iteration")
plt.ylabel("||v||_2")
plt.legend()
plt.grid(True)
plt.savefig("plots/norm_v_random.png", dpi=200)