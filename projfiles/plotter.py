import numpy as np
import matplotlib.pyplot as plt
import sys

omega = float(sys.argv[1]) if len(sys.argv) > 1 else "?"

data = np.loadtxt("data.csv", delimiter=",", ndmin=2)

alpha = data[:, 0]
E = data[:, 1]

# Energy plot
i = np.argmin(E)
plt.figure(figsize=(8, 6))
plt.plot(alpha, E, ".", markersize=4, label="Energy")
plt.plot(alpha[i], E[i], "ro",
         label=f"Min: Energy = {E[i]:.4f} at Alpha = {alpha[i]:.4f}")
plt.axvline(alpha[i], color="red", linestyle="--", linewidth=0.5)
plt.axhline(E[i], color="red", linestyle="--", linewidth=0.5)

plt.title(f"Energy vs. Alpha for ω = {omega}", fontsize=16)
plt.xlabel("Alpha", fontsize=14)
plt.ylabel("Energy (a.u.)", fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

# Separate variance plot, if the third column exists
if data.shape[1] >= 3:
    VarE = data[:, 2]
    var_opt = VarE[i]

    plt.figure(figsize=(8, 6))
    plt.plot(alpha, VarE, "x", markersize=3, color="orange", label="Var(E)")
    plt.plot(alpha[i], var_opt, "ro",
             label=f"Var(E) = {var_opt:.6f} at optimal Alpha = {alpha[i]:.4f}")
    plt.axvline(alpha[i], color="red", linestyle="--", linewidth=0.5)
    plt.axhline(var_opt, color="red", linestyle="--", linewidth=0.5)

    plt.title(f"Variance vs. Alpha for ω = {omega}", fontsize=16)
    plt.xlabel("Alpha", fontsize=14)
    plt.ylabel("Var(E) (a.u.)", fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

plt.show()