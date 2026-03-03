import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("energies.txt")
times = data[:, 0]
Ek    = data[:, 1]
Ep    = data[:, 2]
Etot  = data[:, 3]

plt.figure(figsize=(10, 6))
plt.plot(times, Ek, label="Ek")
plt.plot(times, Ep, label="Ep")
plt.plot(times, Etot, label="Etot", lw=2, color="black")
plt.xlabel("t")
plt.ylabel("E")
plt.title("E(t)")
plt.legend()
plt.tight_layout()
plt.savefig("energy.png", dpi=150)
print("Saved energy.png")
