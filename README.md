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
and
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

For that, we perform a Variational Monte Carlo calculation of the ground state of two electrons in a quantum dot well with different oscillator energies, assuming total spin $s=0$. Compute:

$$
\langle H\rangle = \int P(\mathbf{R}) E_{L}(\mathbf{R}) d\mathbf{R}
$$

where

$$
P(\mathbf{R}) = \frac{\left|\Psi(R)\right|^2}{\int\left|\Psi(R)\right|^2 d\mathbf{R}} \quad \text{and} \quad E_L=\frac{\hat{H}\Psi(R,\theta)}{\Psi(R,\theta)}
$$

where $E_L$ is the local energy. The programme calculates the energy by averaging the local energy over a number of Monte Carlo samples (and hence may also fluctuate slightly between runs).

Our Monte Carlo moves are determined by $\mathbf{R'} = \mathbf{R} + \delta r$, where $\delta$ is a fixed step size for a particular $\alpha$, $r$ is a random number from the uniform distribution, and $\mathbf{R}$=($x_1, y_1, z_1, x_2, y_2, x_2$).

We first use this to find the energy minima and the corresponding alpha for the trial function $\Psi_{T1}$. We also find the expectation value of the mean distance between the electrons at the energy minimum for the optimal set of variational parameters. This is done for $\omega$ = 0.01, 0.5 and 1.0

We thereafter use the optimal value of $\alpha$ as the starting point for computing the ground state energy of the system for the trial function $\Psi_{T2}(\mathbf{r_1,r_2,r_{12}})$, and its corresponding local energy $E_{L2}$. The strategy is to use the optimal $\alpha$, vary $\beta$ to find the lowest energy, and then change $\alpha$ to see if we can find an even lower energy, and so on. This is again done for $\omega$ = 0.01, 0.5 and 1.0. Every time $\alpha$ or $\beta$ is changed, the value is accepted if it gives a smaller value of energy.

## Outputs

### Part 1: Minimum Energy for Non-interacting System
The ground state energy for two-electron, non-interacting system is $3\omega$ . This acts as a benchmark for the code.

Running `projpart1.c` gives $E_{min}$=2.999771 corresponding to $\alpha$=1.011371 at $\omega$=1.0, with acceptance rate of 48.47%.


### Part 2: Plotting $E$ vs. $\alpha$ for $\Psi_{T1}$ and $E_{L1}$
The code `projpart2.c` calculates values of $E$ at fixed intervals of $\alpha$ and plots the graph using `plotter.py`

#### Graph
<!--link to graph image-->

### Part 3: Calculating $E_{min}$, Var($E$) and $\langle\mathbf{r_{12}}\rangle$ for Uncorrelated System
The code `projpart3.c` calculates the values of $E_{min}$, Var($E$) and $\langle\mathbf{r_{12}}\rangle$ for $\Psi_{T1}$ and $E_{L1}$. Then we plot $E$ and Var($E$) vs. $\alpha$ using `plotter.py`

Results (with small margin of error):

| $\omega$ | Optimal $\alpha$ | Estimated $E_{min}$ | $\langle \mathbf{r_{12}} \rangle$ | $Var(E_{min})$ |
|----------|------------------|---------------------|-----------------------------------|----------------|
| 0.01     | 0.445226         | 0.093426            | 23.895274                         | 0.001315       |
| 0.5      | 0.825126         | 2.038693            | 2.496276                          | 0.112110       |
| 1.0      | 0.891457         | 3.769500            | 1.690311                          | 0.247102       | 

#### Graphs of $E$ and Var($E$) vs. $\alpha$
| $\omega$ | 0.01 | 0.5 | 1.0 |
|----------|------|-----|-----|
| $E$ vs. $\alpha$ | <!--link to graph image--> | <!--link to graph image--> | <!--link to graph image--> |
| Var($E$) vs. $\alpha$ | <!--link to graph image--> | <!--link to graph image--> | <!--link to graph image--> |

### Part 4: Calculating $E_{min}$, Var($E$) and $\langle\mathbf{r_{12}}\rangle$ for Uncorrelated System
The code `projpart4.c` calculates the most optimal values of $\alpha$ and $\beta$ for each of the three $\omega$'s. At those values, it calculates ground state energy $E$, Var($E$), and mean separation between electrons $\langle\mathbf{r_{12}}\rangle$. The initial values of $\alpha$ are taken from the previous part.

Results for the trial function $\Psi_{T2}$ (with small margin of error):

| $\omega$ | Optimal $\beta$ | Optimal $\alpha$ | Estimated $E_{min}$ | $\langle \mathbf{r_{12}} \rangle$ | $Var(E_{min})$ |
|----------|-----------------|------------------|---------------------|-----------------------------------|----------------|
| 0.01     | 0.048438        | 0.878878         | 0.079359            | 29.206342                         | 0.000004       |
| 0.5      | 0.212500        | 0.993868         | 2.000076            | 2.671733                          | 0.000144       |
| 1.0      | 0.330078        | 0.991940         | 3.730662            | 1.795705                          | 0.001040       |

### Part 5: (Extra) Plotting graphs of $E$ vs $\alpha$ at optimal $\beta$, $E$ vs $\beta$ at optimal $\alpha$
We plot the graphs of $E$ vs $\alpha$ at optimal $\beta$ and $E$ vs $\beta$ at optimal $\alpha$ for each $\omega$. Then we fit a 4-degree polynomial to the resultant data points. Note: For $\omega$=0.01, since graph of $E$ vs $\beta$ is not well-behaving, we only fit the polynomial to the 150 points around the energy minimum

#### Graphs
| $\omega$ | 0.01 | 0.5 | 1.0 |
|----------|------|-----|-----|
| $E$ vs. $\alpha$, $\beta_{optimal}$ | <!--link to graph image--> | <!--link to graph image--> | <!--link to graph image--> |
| $E$ vs. $\beta$, $\alpha_{optimal}$ | <!--link to graph image--> | <!--link to graph image--> | <!--link to graph image--> |

## Results
It can be seen that upon introducing the Jastrow factor-

$$
\exp\left(\frac{r_{12}}{2(1 + \beta r_{12})}\right)
$$

to the uncorrelated wavefunction $\Psi_{T1}$, we get the following changes:

| $\omega$ | Estimated $E_{min}$ | | $\langle \mathbf{r_{12}} \rangle$ |
|----------|---------------------|-|-----------------------------------|
| 0.01 | 0.093426 $\to$ 0.079359 | | 23.895274 $\to$ 29.206342 |
| 0.5  | 2.038693 $\to$ 2.000076 | | 2.496276 $\to$ 2.671733 |
| 1.0  | 3.769500 $\to$ 3.730662 | | 1.690311 $\to$ 1.795705 |




