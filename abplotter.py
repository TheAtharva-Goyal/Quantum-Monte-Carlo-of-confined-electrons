import numpy as np
import matplotlib.pyplot as plt
import sys

DEGREE = 4  # degree of the polynomial fit

omega = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0

alpha_data = np.loadtxt("energy_vs_alpha.csv", delimiter=",", ndmin=2)
beta_data = np.loadtxt("energy_vs_beta.csv", delimiter=",", ndmin=2)

# Keep only the selected frequency.
alpha_data = alpha_data[np.isclose(alpha_data[:, 0], omega)]
beta_data = beta_data[np.isclose(beta_data[:, 0], omega)]

if len(alpha_data) == 0 or len(beta_data) == 0:
    sys.exit(f"No data found for omega = {omega}")


def mark(x, y, colour, marker, label, size=36, zorder=5):
    """Draw a marked point together with dashed guide lines to both axes."""
    plt.scatter(x, y, color=colour, marker=marker, s=size, label=label, zorder=zorder)
    plt.axvline(x, color=colour, linestyle="--", linewidth=0.5)
    plt.axhline(y, color=colour, linestyle="--", linewidth=0.5)


def fit_minimum(x, y):
    """Fit a polynomial to (x, y); return the curve and its lowest local minimum.

    Returns x_fit, y_fit, x_min, E_min (x_min and E_min are None if the fit
    has no local minimum inside the fitted interval).
    """
    fit = np.polynomial.Polynomial.fit(x, y, deg=DEGREE)
    x_lo, x_hi = x.min(), x.max()

    x_fit = np.linspace(x_lo, x_hi, 500)
    y_fit = fit(x_fit)

    d2 = fit.deriv(2)
    candidates = [
        root.real
        for root in fit.deriv().roots()
        if abs(root.imag) < 1e-8
        and x_lo <= root.real <= x_hi
        and d2(root.real) > 0
    ]
    if not candidates:
        return x_fit, y_fit, None, None

    x_min = min(candidates, key=fit)
    return x_fit, y_fit, x_min, fit(x_min)


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

    i = np.argmin(y)                      # lowest sampled point
    j = np.argmin(np.abs(x - final_value))  # sample closest to the final value

    plt.figure(figsize=(8, 6))
    plt.plot(x, y, ".", markersize=2, label="Data points")

    # Choose the fitting strategy.
    if np.isclose(omega, 0.01) and variable == "Beta":
        idx = np.argsort(np.abs(x - final_value))[:151]
        fit_label = "Local degree-4 polynomial fit"
    else:
        idx = np.arange(len(x))
        fit_label = "Degree-4 polynomial fit"

    x_fit, y_fit, x_min, E_min = fit_minimum(x[idx], y[idx])
    plt.plot(x_fit, y_fit, color="black", linestyle="--", linewidth=1,
             alpha=0.65, label=fit_label, zorder=4)

    if x_min is not None:
        mark(x_min, E_min, "magenta", "s",
             f"Fitted minimum: {variable}={x_min:.6f}", size=65, zorder=10)

    mark(x[i], y[i], "red", "o",
         f"Lowest sampled E = {y[i]:.6f} at {variable} = {x[i]:.6f}")
    mark(x[j], y[j], "green", "o",
         f"Final optimized {variable} = {final_value:.6f}")

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