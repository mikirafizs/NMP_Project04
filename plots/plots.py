import numpy as np
import matplotlib.pyplot as plt

omegas = [1, 1.5, 1.9, 1.99]

plt.figure()
for omega in omegas:
    d = np.loadtxt(f"data/sor_omega_{omega}.csv", delimiter=",")
    plt.plot(d[:, 0], d[:, 2], label=f"omega = {omega}")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("iteration")
plt.ylabel("||r||_2")
plt.legend()
plt.grid(True)
plt.savefig("plots/norm_r_all.png", dpi=200)

plt.figure()
for omega in omegas:
    d = np.loadtxt(f"data/sor_omega_{omega}.csv", delimiter=",")
    plt.plot(d[:, 0], d[:, 1], label=f"omega = {omega}")
plt.xscale("log")
plt.xlabel("iteration")
plt.ylabel("||v||_2")
plt.legend()
plt.grid(True)
plt.savefig("plots/norm_v_all.png", dpi=200)

vec = np.loadtxt("data/vector_omega_1.csv", delimiter=",")
plt.figure()
plt.plot(vec[:, 0], vec[:, 1])
plt.xlabel("x")
plt.ylabel("rho")
plt.title("Source term rho(x)")
plt.grid(True)
plt.savefig("plots/rho.png", dpi=200)

plt.figure()
for omega in omegas:
    vec = np.loadtxt(f"data/vector_omega_{omega}.csv", delimiter=",")
    plt.plot(vec[:, 0], vec[:, 2], label=f"omega = {omega}")
plt.xlabel("x")
plt.ylabel("v")
plt.title("Final solution v(x)")
plt.legend()
plt.grid(True)
plt.savefig("plots/v_final_all.png", dpi=200)