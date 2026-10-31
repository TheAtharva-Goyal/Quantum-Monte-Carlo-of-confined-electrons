import numpy as np
import matplotlib.pyplot as plt
import sys

omega = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0

alpha_data = np.atleast_2d(
    np.loadtxt("energy_vs_alpha.csv", delimiter=",", skiprows=1)
)
beta_data = np.atleast_2d(
    np.loadtxt("energy_vs_beta.csv", delimiter=",", skiprows=1)
)

# Keep only the selected frequency.
alpha_data = alpha_data[np.isclose(alpha_data[:, 0], omega)]
beta_data = beta_data[np.isclose(beta_data[:, 0], omega)]

if len(alpha_data) == 0 or len(beta_data) == 0:
    sys.exit(f"No data found for omega = {omega}")


def make_plot(data, variable, fixed):
    # Both files use the same column positions:
    # 0: omega, 1: varying parameter, 2: fixed parameter,
    # 3: final optimized varying parameter, 4: energy
    x = data[:, 1]
    y = data[:, 4]
    fixed_value = data[0, 2]
    final_value = data[0, 3]

    # Sort by the varying parameter.
    order = np.argsort(x)
    x = x[order]
    y = y[order]

    i = np.argmin(y)
    j = np.argmin(np.abs(x - final_value))

    plt.figure(figsize=(8, 6), dpi=100)

    plt.plot(x, y, ".", markersize=2, label="Data points")

    #====FITTING====
    # Choose the fitting strategy.
    if np.isclose(omega, 0.01) and variable == "Beta":
        idx = np.argsort(np.abs(x - final_value))[:151]
        fit_label = "Local degree-4 polynomial fit"
    else:
        idx = np.arange(len(x))
        fit_label = "Degree-4 polynomial fit"

    degree = 4

    fit = np.polynomial.Polynomial.fit(
        x[idx], y[idx], deg=degree
    )

    # Plot over the range used for fitting.
    x_lo = x[idx].min()
    x_hi = x[idx].max()

    x_fit = np.linspace(x_lo, x_hi, 500)
    y_fit = fit(x_fit)

    plt.plot(
        x_fit, y_fit,
        color="black",
        linestyle="--",
        linewidth=1,
        alpha=0.65,
        label=fit_label,
        zorder=4
    )

    # Find and mark fitted local minima inside the fitting interval.
    d1 = fit.deriv()
    d2 = fit.deriv(2)

    roots = d1.roots()
    candidates = [
        root.real
        for root in roots
        if abs(root.imag) < 1e-8
        and x_lo <= root.real <= x_hi
        and d2(root.real) > 0
    ]

    if candidates:
        x_min = min(candidates, key=lambda xm: fit(xm))
        E_min = fit(x_min)

        plt.scatter(
            x_min, E_min,
            color="magenta",
            marker="s",
            s=65,
            label=f"Fitted minimum: {variable}={x_min:.6f}",
            zorder=10
        )
    #---------


    plt.plot(
        x[i], y[i], "ro",
        label=f"Lowest sampled E = {y[i]:.6f} at {variable} = {x[i]:.6f}"
    )
    plt.axvline(x[i], color="red", linestyle="--", linewidth=0.5)
    plt.axhline(y[i], color="red", linestyle="--", linewidth=0.5)

    plt.plot(
        x[j], y[j], "go",
        label=f"Final optimized {variable} = {final_value:.6f}"
    )
    plt.axvline(final_value, color="green",
                linestyle="--", linewidth=0.5)

    plt.title(
        f"Energy vs. {variable} for ω = {omega}\n"
        f"{fixed} fixed at {fixed_value:.6f}",
        fontsize=15
    )
    plt.xlabel(variable, fontsize=14)
    plt.ylabel("Energy (a.u.)", fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.legend()


make_plot(alpha_data, "Alpha", "Beta")
make_plot(beta_data, "Beta", "Alpha")

plt.show()