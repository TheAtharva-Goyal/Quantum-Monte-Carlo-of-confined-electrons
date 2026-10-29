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
We want to perform a Variational Monte Carlo calculation of the ground state of two electrons in a quantum dot well with different oscillator energies, assuming total spin S=0. Compute:

$$
\langle H\rangle = \int P(\mathbf{R}) E_{L}(\mathbf{R}) d\mathbf{R}
$$

where

$$
P(\mathbf{R}) = \frac{\left|\Psi(R)\right|^2}{\int\left|\Psi(R)\right|^2 d\mathbf{R}}
$$

and $E_L$ is the local energy. The programme calculates the energy by averaging the local energy over a number of Monte Carlo samples (and hence may also fluctuate slightly between runs).

The Monte Carlo moves are determined by $\mathbf{R'} = \mathbf{R} + \delta r$, where $\delta$ is a fixed step size for a particular $\alpha$, and $r$ is a random number from the uniform distribution.

We first use this to find the energy minima and the corresponding alpha for the trial function $\Psi_{T1}$. We also find the expectation value of the mean distance between the electrons at the energy minimum for the optimal set of variational parameters. This is done for $\omega$ = 0.01, 0.5 and 1.0

We thereafter use the optimal value of $\alpha$ as the starting point for computing the ground state energy of the system for the trial function $\Psi_{T2}(\mathbf{r_1,r_2,r_{12}})$, and its corresponding local energy $E_{L2}$. The strategy is to use the optimal $\alpha$, vary $\beta$ to find the lowest energy, and then change $\alpha$ to see if we can find an even lower energy, and so on. This is again done for $\omega$ = 0.01, 0.5 and 1.0. Every time $\alpha$ or $\beta$ is changed, the value is accepted if it gives a smaller value of energy.

## Results
The programme runs with $(50{\pm}0.2)$% acceptance.

It finds the most optimal values of $\alpha$ and $\beta$ for each $\omega$. At those values, it calculates ground state energy $E$, Var($E$), and mean separation between electrons $\langle \mathbf{r_{12}} \rangle$.

Results for the trial function $\Psi_{T1}$ (With some margin of error):

| $\omega$ | Optimal $\alpha$ | Estimated energy $E$ | Mean separation $\langle \mathbf{r_{12}} \rangle$ |
|----------|------------------|----------------------|---------------------------------------------------|
| 0.01     | 0.454271         | 0.093298             | 23.670367                                         |
| 0.50     | 0.843216         | 2.037538             | 2.459868                                          |
| 1.00     | 0.885427         | 3.770180             | 1.692222                                          |


It also plots the graphs of $E$ vs $\alpha$ at optimal $\beta$ and $E$ vs $\beta$ at optimal $\alpha$ for each $\omega$, and fits a polynomial to the resultant data points.
