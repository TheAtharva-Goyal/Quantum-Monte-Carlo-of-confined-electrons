# Quantum-Monte-Carlo-of-confined-electrons
The aim of this project is to use the Variational Monte Carlo method to evaluate the ground state energy and the relative distance between two electrons for a quantum dot with two electrons in three dimensions.
Our aim here is to use Variational Monte Carlo method to study such systems, with an emphasis on understanding correlations due to the repulsive interaction between electrons.

We consider a system of electrons confined in a pure 3D harmonic oscillator potential, with an idealised total Hamiltonian given by:

$$
\hat{H} = \sum_{i=1}^{2} \left(-\frac{1}{2}\nabla_i^2 + \frac{1}{2}\omega^2 r_i^2 \right) + \frac{1}{r_{12}}
$$

where all natural units ($\hbar$ = c = e = $m_e$ = 1) are used. (energies are in atomic units, distances are dimensionless)
