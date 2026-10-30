# Quantum Monte Carlo of confined electrons

![License](https://img.shields.io/badge/license-MIT-green)

The aim of this project is to use the Variational Monte Carlo method to evaluate the ground state energy and the relative distance between two electrons for a quantum dot with two electrons in three dimensions.
Our aim here is to use Variational Monte Carlo method to study such systems, with an emphasis on understanding correlations due to the repulsive interaction between electrons.

We consider a system of electrons confined in a pure 3D harmonic oscillator potential of frequency $\omega$, with an idealised total Hamiltonian given by:

$$
\hat{H} = \sum_{i=1}^{2} \left(-\frac{1}{2}\nabla_i^2 + \frac{1}{2}\omega^2 r_i^2 \right) + \frac{1}{r_{12}}
$$

where all natural units ($\hbar$ = c = e = $m_e$ = 1) are used. (energies are in atomic units, distances are dimensionless)

### Trial functions and corresponding energies
We use trial functions:

$$
\Psi_{T1}(\mathbf{r}_1, \mathbf{r}_2) = C \exp\left(-\frac{\alpha\omega}{2} \left(r_1^2 + r_2^2\right)\right),
$$

$$
\Psi_{T2}(\mathbf{r}_1, \mathbf{r}_2) = C \exp\left(-\frac{\alpha\omega}{2} \left(r_1^2 + r_2^2\right)\right) \exp\left(\frac{r_{12}}{2(1 + \beta r_{12})}\right)
$$

With local energy (including Coulomb interaction) corresponding to $\Psi_{T1}$ given by:

$$
E_{L1} = \frac{1}{2} \omega^2 \left(r_1^2 + r_2^2\right) \left(1-\alpha^2\right) + 3\alpha\omega + \frac{1}{r_{12}}
$$

And local energy corresponding to $\Psi_{T2}$ given by:

$$
E_{L2} = E_{L1} + \frac{1}{2(1+\beta r_{12})^2} \left[\alpha\omega r_{12} - \frac{1}{2(1+\beta r_{12})^2} - \frac{2}{r_{12}} + \frac{2\beta}{1+\beta r_{12}}\right]
$$

## Method
Given the Hamiltonian $\hat{H}$, and the wavefunction depending on some parameters $\mathbf{\theta}$, the expectation of the energy is given by:

$$
E(\theta)=\frac{\langle\Psi(\theta)|\hat{H}|\Psi(\theta)\rangle}{\langle\Psi(\theta)|\Psi(\theta)\rangle} = \frac{\int|\Psi(R,\theta)|^2 \frac{\hat{H}\Psi(R,\theta)}{\Psi(R,\theta)} dR}{\int|\Psi(R,\theta)|^2 dR}
$$

Upon solving this, the energy we get is $E(\theta)\geq E_0$ where the (unknown) ground state energy $E_0$ is defined by $\hat{H}|i_0\rangle=E_0|i_0\rangle$

We try to minimise the value of $E(\theta)$ to make it as close as possible to the true energy $E_0$.

For that, we perform a Variational Monte Carlo calculation of the ground state of two electrons in a quantum dot well with different oscillator energies, taking total spin $s=0$ (since both electrons occupy the same lowest energy orbital, they must have opposite spins to obey Pauli's exclusion principle). Compute:

$$
\langle H\rangle = \int P(\mathbf{R}) E_{L}(\mathbf{R}) d\mathbf{R}
$$

where

$$
P(\mathbf{R}) = \frac{\left|\Psi(R)\right|^2}{\int\left|\Psi(R)\right|^2 d\mathbf{R}} \quad \text{and} \quad E_L=\frac{\hat{H}\Psi(R,\theta)}{\Psi(R,\theta)}
$$

where $E_L$ is the local energy. The programme calculates the energy by averaging the local energy over a number of Monte Carlo samples (and hence may also fluctuate slightly between runs).

Our Monte Carlo moves are determined by $\mathbf{R'} = \mathbf{R} + \delta r$, where $\delta$ is a fixed step size, $r$ is a random number from the uniform distribution, and $\mathbf{R}$=($x_1, y_1, z_1, x_2, y_2, z_2$).

We first use this to find the energy minima and the corresponding alpha for the trial function $\Psi_{T1}$. We also find the expectation value of the mean distance between the electrons at the energy minimum for the optimal set of variational parameters. This is done for $\omega$ = 0.01, 0.5 and 1.0

We thereafter use the optimal value of $\alpha$ as the starting point for computing the ground state energy of the system for the trial function $\Psi_{T2}(\mathbf{r_1,r_2,r_{12}})$, and its corresponding local energy $E_{L2}$. The strategy is to use the optimal $\alpha$, vary $\beta$ to find the lowest energy, and then change $\alpha$ to see if we can find an even lower energy, and so on. This is again done for $\omega$ = 0.01, 0.5 and 1.0. Every time $\alpha$ or $\beta$ is changed, the value is accepted if it gives a smaller value of energy.

## Outputs

### Part 1: Minimum Energy for Non-interacting System
The ground state energy for two-electron, non-interacting system is $3\omega$ (since the ground state energy for a single electron is $\frac{3\omega}{2}$) at $\alpha=1.0$. Hence, the exact result is $E=3$ atomic units at $\omega=1.0$. This acts as a benchmark for the code.

Running `projpart1.c` gives $E_{min}=2.999771$ corresponding to $\alpha=1.011371$ at $\omega=1.0$, with acceptance rate of 48.47%.


### Part 2: Plotting $E$ vs. $\alpha$ for $\Psi_{T1}$ and $E_{L1}$
The code `projpart2.c` calculates values of $E$ at uniform intervals of $\alpha$ and plots the values using `plotter.py`

#### Graph
<!--link to graph image-->


### Part 3: Calculating $E_{min}$, Var($E$) and $\langle\mathbf{r_{12}}\rangle$ for Uncorrelated System
The code `projpart3.c` calculates the values of $E_{min}$, Var($E$) and $\langle\mathbf{r_{12}}\rangle$ for $\Psi_{T1}$ and $E_{L1}$. Then we plot $E$ and Var($E$) vs. $\alpha$ using `plotter.py`

Note: The programme runs with acceptance of ~50% and gives standard deviation of Energy $< 6.37\times 10^{-4}$ and standard deviation of Separation $< 2.06\times 10^{-2}$ for all three values of $\omega$ in 100 runs of 5 million Monte Carlo steps each.

The deviation of calculated $E_{min}$ from analytical expression for ground state energy for trial function $\Psi_{T1}$

$$
E_{ana}(\alpha)=\frac{3\omega}{2}\left(\alpha +\frac{1}{\alpha}\right)+\sqrt{\frac{2\alpha\omega}{\pi}}
$$

in 100 runs of 5 million Monte Carlo steps each is $< 6.49\times 10^{-4}$

Results:

| $\omega$ | Optimal $\alpha$ | Estimated $E_{min}$ | $\langle\mathbf{r_{12}}\rangle$ | $Var(E_{min})$ |
|----------|------------------|---------------------|---------------------------------|----------------|
| 0.01     | 0.445226         | 0.093426            | 23.895274                       | 0.001315       |
| 0.5      | 0.825126         | 2.038693            | 2.496276                        | 0.112110       |
| 1.0      | 0.891457         | 3.769500            | 1.690311                        | 0.247102       | 

#### Graphs
| $\omega$ | 0.01 | 0.5 | 1.0 |
|----------|------|-----|-----|
| $E$ vs. $\alpha$ | <!--link to graph image--> | <!--link to graph image--> | <!--link to graph image--> |
| Var($E$) vs. $\alpha$ | <!--link to graph image--> | <!--link to graph image--> | <!--link to graph image--> |


### Part 4: Calculating $E_{min}$, Var($E$) and $\langle\mathbf{r_{12}}\rangle$ for Correlated System
The code `projpart4.c` calculates the most optimal values of $\alpha$ and $\beta$ for each of the three $\omega$'s. At those values, it calculates ground state energy $E$, Var($E$), and mean separation between electrons $\langle\mathbf{r_{12}}\rangle$. The initial values of $\alpha$ are taken from the previous part.

Note: The programme runs with acceptance of ~50% and gives standard deviation of Energy $< 2.65\times 10^{-5}$ and standard deviation of Separation $< 9.43\times 10^{-3}$ for all three values of $\omega$ in 100 runs of 5 million Monte Carlo steps each.

Results for the trial function $\Psi_{T2}$:

| $\omega$ | Optimal $\beta$ | Optimal $\alpha$ | Estimated $E_{min}$ | $\langle\mathbf{r_{12}}\rangle$ | $Var(E_{min})$ |
|----------|-----------------|------------------|---------------------|---------------------------------|----------------|
| 0.01     | 0.048438        | 0.878878         | 0.079359            | 29.206342                       | 0.000004       |
| 0.5      | 0.212500        | 0.993868         | 2.000076            | 2.671733                        | 0.000144       |
| 1.0      | 0.330078        | 0.991940         | 3.730662            | 1.795705                        | 0.001040       |


### Part 5: (Extra) Plotting graphs of $E$ vs $\alpha$ at optimal $\beta$, $E$ vs $\beta$ at optimal $\alpha$
We plot the graphs of $E$ vs $\alpha$ at optimal $\beta$ and $E$ vs $\beta$ at optimal $\alpha$ for each $\omega$. Then we fit a 4-degree polynomial to the resultant data points. Note: For $\omega$=0.01, since graph of $E$ vs $\beta$ is not well-behaving, we only fit the polynomial to the 150 points around the energy minimum

#### Graphs
| $\omega$ | 0.01 | 0.5 | 1.0 |
|----------|------|-----|-----|
| $E$ vs. $\alpha$, $\beta=\beta_{optimal}$ | <!--link to graph image--> | <!--link to graph image--> | <!--link to graph image--> |
| $E$ vs. $\beta$, $\alpha=\alpha_{optimal}$ | <!--link to graph image--> | <!--link to graph image--> | <!--link to graph image--> |

## Results
It can be seen that upon introducing the Jastrow factor

$$
\exp\left(\frac{r_{12}}{2(1 + \beta r_{12})}\right)
$$

to the uncorrelated wavefunction $\Psi_{T1}$, we get the following changes:

| $\omega$ | Estimated $E_{min}$     | $Var(E_{min})$          | $\langle\mathbf{r_{12}}\rangle$ |
|----------|-------------------------|-------------------------|---------------------------------|
| 0.01     | 0.093426 $\to$ 0.079359 | 0.001315 $\to$ 0.000004 | 23.895274 $\to$ 29.206342       |
| 0.5      | 2.038693 $\to$ 2.000076 | 0.112110 $\to$ 0.000144 | 2.496276 $\to$ 2.671733         |
| 1.0      | 3.769500 $\to$ 3.730662 | 0.247102 $\to$ 0.001040 | 1.690311 $\to$ 1.795705         |

#### Why the Jastrow factor decreases $E_{min}$ and increases $\langle\mathbf{r_{12}}\rangle$
Introducing the Jastrow factor to $\Psi_{T1}$ accounts for the inter-electron repulsion, which is obviously a more accurate description of the system. Hence we are able to model the ground state energy better, which gives a lower (and closer to actual) estimate for the energy. This can be seen by the fact that the variance of minimum energy drops substantially, which indicates that the calculated energy is more stable.

We get an increase in the mean inter-electron distance upon accounting for their interaction since the electrons will repel each other, and hence will tend to stay further apart.

Thus it can be seen that the Jastrow factor improves the uncorrelated trial wavefunction.

## Conclusion
We used the Variational Monte Carlo method on the system of two electrons in a 3-D quantum dot to estimate their ground state energy and mean separation at three values of oscillator frequency $\omega$. We used two trial wavefunctions, $\Psi_{T1}(\mathbf{R},\alpha)$ and $\Psi_{T2}(\mathbf{R},\alpha,\beta)$, where the second (correlated) function incorporates a Jastrow factor to account for the inter-electron repulsion.

For both trial functions, $E_{min}$ is lowest for $\omega=0.01$ since due to the low strength of confinement, the electrons will tend to wander farther from the origin. For $\omega=1.0$, since the strength of confinement is higher, the electrons will tend to gather closer to the origin, hence the higher energy.

Similarly, $\langle\mathbf{r_{12}}\rangle$ is highest for  $\omega=0.01$ since due to the low strength of confinement, the electrons will tend to wander outwards more and hence will tend to be spaced further apart. For $\omega=1.0$, the higher strength of confinement will bring the electrons closer to the origin, and subsequently, each other.

The higher and lower values for $E_{min}$ and $\langle\mathbf{r_{12}}\rangle$ respectively for the second trial function compared to the first is explained by the role of the Jastrow factor. Since it accounts for the inter-electron repulsion, it better models the system, and hence the VMC estimate of $E_{min}$ is lower. This inter-electron repulsion also leads to the electrons tending to stay comparatively further from each other, leading to an increased value of $\langle\mathbf{r_{12}}\rangle$ (especially when $\omega=0.01$, i.e. when the confinement is weak).

The reliability of the results for $\Psi_{T1}$ were calculated by running 100 loops of 5 million steps each. For the three values of $\omega$, the deviation from the analytical value of $E_{min}$ was less than $6.49\times 10^{-4}$, and the inter-run standard deviation for Energy and Separation was less than $6.37\times 10^{-4}$ and $2.06\times 10^{-2}$ respectively.

Similarly, 100 loops of 5 million steps each were ran for $\Psi_{T2}$. Since no suitable analytical expression exists for it, we calculated only the inter-run standard deviations for Energy and Separation. These were respectively less than $2.65\times 10^{-5}$ and $9.43\times 10^{-3}$.

---
