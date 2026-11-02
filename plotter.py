import numpy as np
import matplotlib.pyplot as plt
import sys

# usage: uv run plotter.py <omega>
# data.csv columns: 3 columns -> alpha, E, omega
#                   4 columns -> alpha, E, Var(E), omega

if len(sys.argv) < 2:
    sys.exit("usage: uv run plotter.py <omega>   (e.g. 0.01, 0.5 or 1.0)")
omega = float(sys.argv[1])

data = np.loadtxt("data.csv", delimiter=",", ndmin=2)
ncols = data.shape[1]
if ncols not in (3, 4):
    sys.exit(f"data.csv has {ncols} columns; expected 3 (alpha,E,omega) or 4 (alpha,E,Var,omega)")

# keep only the rows belonging to the requested omega (last column)
data = data[np.isclose(data[:, -1], omega)]
if data.shape[0] == 0:
    sys.exit(f"no rows found for omega = {omega}")

alpha = data[:, 0]
E = data[:, 1]

# Energy plot
i = np.argmin(E)
plt.figure(figsize=(8, 6))
plt.plot(alpha, E, ".", markersize=4, label="Energy")
plt.plot(alpha[i], E[i], "ro",
         label=f"Min: Energy = {E[i]:.6g} at Alpha = {alpha[i]:.4f}")
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

# Separate variance plot (4-column files only)
if ncols == 4:
    VarE = data[:, 2]
    var_opt = VarE[i]  # variance at the energy minimum

    plt.figure(figsize=(8, 6))
    plt.plot(alpha, VarE, "x", markersize=3, color="orange", label="Var(E)")
    plt.plot(alpha[i], var_opt, "ro",
             label=f"Var(E) = {var_opt:.6g} at optimal Alpha = {alpha[i]:.4f}")
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