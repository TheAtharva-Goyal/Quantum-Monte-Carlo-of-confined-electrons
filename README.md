# Quantum Monte Carlo of confined electrons
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
