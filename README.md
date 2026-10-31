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

### Why numerical methods are necessary
While the analytical form of energy for $E_{L1}$ exists, no such universal analytical form exists for $E_{L2}$, which is why we use the Monte Carlo method to calculate ground state energy.

$$
E_{analytical}(\alpha)=\frac{3\omega}{2}\left(\alpha +\frac{1}{\alpha}\right)+\sqrt{\frac{2\alpha\omega}{\pi}}
$$

The primary shortcoming of analytical approach on expectation energy corresponding to $\Psi_{T2}$ is that the subsequent integral is generally non-integrable due to the Coulomb interaction coupling the radial and axial coordinates (although the referenced article by M. Taut shows that for a particular denumerably infinite set of oscillator frequencies, analytical solutions exist). Hence, we are motivated to use the numerical approach.

## Method
Given the Hamiltonian $\hat{H}$, and the wavefunction depending on some parameters $\mathbf{\theta}=\left(\alpha,\beta,\gamma,...\right)$, the expectation of the energy is given by:

$$
E(\theta)=\frac{\langle\Psi(\theta)|\hat{H}|\Psi(\theta)\rangle}{\langle\Psi(\theta)|\Psi(\theta)\rangle} = \frac{\int|\Psi(R,\theta)|^2 \frac{\hat{H}\Psi(R,\theta)}{\Psi(R,\theta)} dR}{\int|\Psi(R,\theta)|^2 dR}
$$

Upon solving this, the energy we get is $E(\theta)\geq E_0$ where the (unknown) ground state energy $E_0$ is defined by $\hat{H}|i_0\rangle=E_0|i_0\rangle$

We try to minimise the value of $E(\theta)$ to get it as close as possible to the true energy $E_0$.

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

#### Graphs
![Part2_w=0.01](https://github.com/TheAtharva-Goyal/Quantum-Monte-Carlo-of-confined-electrons/blob/main/graphs/projpart2/EvsA_w001.jpeg)

![Part2_w=0.5](https://github.com/TheAtharva-Goyal/Quantum-Monte-Carlo-of-confined-electrons/blob/main/graphs/projpart2/EvsA_w050.jpeg)

![Part2_w=1.0](https://github.com/TheAtharva-Goyal/Quantum-Monte-Carlo-of-confined-electrons/blob/main/graphs/projpart2/EvsA_w100.jpeg)


### Part 3: Calculating $E_{min}$, Var($E$) and $\langle\mathbf{r_{12}}\rangle$ for Uncorrelated System
The code `projpart3.c` calculates the values of $E_{min}$, Var($E$) and $\langle\mathbf{r_{12}}\rangle$ for $\Psi_{T1}$ and $E_{L1}$. Then we plot $E$ and Var($E$) vs. $\alpha$ using `plotter.py`

Note: The programme runs with acceptance of ~50% and gives standard deviation of Energy $< 6.37\times 10^{-4}$ (not related to $Var(E_{min})$ ) and standard deviation of Separation $< 2.06\times 10^{-2}$ for all three values of $\omega$ in 100 runs of 5 million Monte Carlo steps each at the minima.

The deviation of calculated $E_{min}$ from analytical expression of energy expectation for trial function $\Psi_{T1}$ in 100 runs of 5 million Monte Carlo steps each at the minima is $< 6.50\times 10^{-4}$

Results:

| $\omega$ | Optimal $\alpha$ | Estimated $E_{min}$ | $\langle\mathbf{r_{12}}\rangle$ | $Var(E_{min})$ |
|----------|------------------|---------------------|---------------------------------|----------------|
| 0.01 | 0.437310  | 0.093629        | 24.130319 | 0.001365   |
| 0.5  | 0.817296  | 2.040659        | 2.496515  | 0.126559   |
| 1.0  | 0.878946  | 3.773088        | 1.702181  | 0.269931   |

#### Graphs
| $\omega$ | 0.01 | 0.5 | 1.0 |
|----------|------|-----|-----|
|  | ![Part3_E_w=0.01](https://github.com/TheAtharva-Goyal/Quantum-Monte-Carlo-of-confined-electrons/blob/main/graphs/projpart3/EvsA_w001.jpeg) | ![Part3_E_w=0.5](https://github.com/TheAtharva-Goyal/Quantum-Monte-Carlo-of-confined-electrons/blob/main/graphs/projpart3/EvsA_w050.jpeg) | ![Part3_E_w=1.0](https://github.com/TheAtharva-Goyal/Quantum-Monte-Carlo-of-confined-electrons/blob/main/graphs/projpart3/EvsA_w100.jpeg) |
| Var($E$) vs. $\alpha$ | ![Part3_V_w=0.01](https://github.com/TheAtharva-Goyal/Quantum-Monte-Carlo-of-confined-electrons/blob/main/graphs/projpart3/VvsA_w001.jpeg) | ![Part3_V_w=0.5](https://github.com/TheAtharva-Goyal/Quantum-Monte-Carlo-of-confined-electrons/blob/main/graphs/projpart3/VvsA_w050.jpeg) | ![Part3_V_w=1.0](https://github.com/TheAtharva-Goyal/Quantum-Monte-Carlo-of-confined-electrons/blob/main/graphs/projpart3/VvsA_w100.jpeg) |


### Part 4: Calculating $E_{min}$, Var($E$) and $\langle\mathbf{r_{12}}\rangle$ for Correlated System
The code `projpart4.c` calculates the most optimal values of $\alpha$ and $\beta$ for each of the three $\omega$'s. At those values, it calculates ground state energy $E$, Var($E$), and mean separation between electrons $\langle\mathbf{r_{12}}\rangle$. The initial values of $\alpha$ are taken from the previous part.

Note: The programme runs with acceptance of ~50% and gives standard deviation of Energy $< 2.65\times 10^{-5}$ (not related to $Var(E_{min})$ ) and standard deviation of Separation $< 9.43\times 10^{-3}$ for all three values of $\omega$ in 100 runs of 5 million Monte Carlo steps each at the minima.

Results for the trial function $\Psi_{T2}$:

| $\omega$ | Optimal $\beta$ | Optimal $\alpha$ | Estimated $E_{min}$ | $\langle\mathbf{r_{12}}\rangle$ | $Var(E_{min})$ |
|----------|-----------------|------------------|---------------------|---------------------------------|----------------|
| 0.01 | 0.044802  | 0.933212  | 0.079295        | 29.178112 | 0.000003   |
| 0.5  | 0.199826  | 0.996432  | 2.000065        | 2.684703  | 0.000164   |
| 1.0  | 0.275797  | 0.999145  | 3.730192        | 1.811207  | 0.000263   |


### Part 5: (Extra) Plotting graphs of $E$ vs $\alpha$ at optimal $\beta$, $E$ vs $\beta$ at optimal $\alpha$
We plot the graphs of $E$ vs $\alpha$ at optimal $\beta$ and $E$ vs $\beta$ at optimal $\alpha$ for each $\omega$. Then we fit a 4-degree polynomial to the resultant data points. Note: For $\omega$=0.01, since graph of $E$ vs $\beta$ is not well-behaving, we only fit the polynomial to the 150 points around the energy minimum

#### Graphs
| $\omega$ | 0.01 | 0.5 | 1.0 |
|----------|------|-----|-----|
| $E$ vs. $\alpha$, $\beta=\beta_{optimal}$ | ![Part5_EA_w=0.01](https://github.com/TheAtharva-Goyal/Quantum-Monte-Carlo-of-confined-electrons/blob/main/graphs/projpart5/EvsA_B_w001.jpeg) | ![Part5_EA_w=0.5](https://github.com/TheAtharva-Goyal/Quantum-Monte-Carlo-of-confined-electrons/blob/main/graphs/projpart5/EvsA_B_w050.jpeg) | ![Part5_EA_w=1.0](https://github.com/TheAtharva-Goyal/Quantum-Monte-Carlo-of-confined-electrons/blob/main/graphs/projpart5/EvsA_B_w100.jpeg) |
| $E$ vs. $\beta$, $\alpha=\alpha_{optimal}$ | ![Part5_EB_w=0.01](https://github.com/TheAtharva-Goyal/Quantum-Monte-Carlo-of-confined-electrons/blob/main/graphs/projpart5/EvsB_A_w001.jpeg) | ![Part5_EB_w=0.5](https://github.com/TheAtharva-Goyal/Quantum-Monte-Carlo-of-confined-electrons/blob/main/graphs/projpart5/EvsB_A_w050.jpeg) | ![Part5_EB_w=1.0](https://github.com/TheAtharva-Goyal/Quantum-Monte-Carlo-of-confined-electrons/blob/main/graphs/projpart5/EvsB_A_w100.jpeg) |

#### Limitation of Variational Monte Carlo
It can be seen that often times the minimum energy calculated by this approach lies below the exact energy given by the analytical approach, which should not be possible since it is established that the calculated energy should always lie above or equal to the analytical energy. It is not entirely right to say that such values are wrong, although they are not right either. This is a shortcoming of this numerical method, that each data point has some amount of noise and error, and taking the minimum over them often leads to picking a data point which lies below the exact energy. Hence, it is necessary to proceed with caution and properly learn the error in the results.

#### Advantage of Variational Monte Carlo method
The main advantage of this method is the freedom to choose complex trial functions which properly model the physical system without having to worry about solving their corresponding integrals and such. Which is because this method mainly requires the value of these functions at different points, which makes it relatively simple compared to something like matrix diagonalisation, and yields better results with significantly less CPU time. It is also easily parallelisable, which makes computation very efficient and quick. Another advantage in numerical integration is that VMC is able to handle the exponentially growing dimensionality of the multi-body Hilbert space (for N particles, it is $3^N$), which is where traditional numerical methods like Gauss-Legendre fail.

## Results
It can be seen that upon introducing the Jastrow factor

$$
\exp\left(\frac{r_{12}}{2(1 + \beta r_{12})}\right)
$$

to the uncorrelated wavefunction $\Psi_{T1}$, we get the following changes:

| $\omega$ | Estimated $E_{min}$     | $Var(E_{min})$          | $\langle\mathbf{r_{12}}\rangle$ |
|----------|-------------------------|-------------------------|---------------------------------|
| 0.01     | 0.093629 $\to$ 0.079295 | 0.001365 $\to$ 0.000003 | 24.130319 $\to$ 29.178112       |
| 0.5      | 2.040659 $\to$ 2.000065 | 0.126559 $\to$ 0.000164 | 2.496515 $\to$ 2.684703         |
| 1.0      | 3.773088 $\to$ 3.730192 | 0.269931 $\to$ 0.000263 | 1.702181 $\to$ 1.811207         |

#### Why the Jastrow factor decreases $E_{min}$ and increases $\langle\mathbf{r_{12}}\rangle$
Introducing the Jastrow factor to $\Psi_{T1}$ accounts for the inter-electron repulsion, which is obviously a more accurate description of the system. Hence we are able to model the ground state energy better, which gives a lower (and closer to actual) estimate for the energy. This can be seen by the fact that the variance of minimum energy drops substantially, which indicates that the calculated energy is more stable.

We get an increase in the mean inter-electron distance upon accounting for their interaction since the electrons will repel each other, and hence will tend to stay further apart.

Thus the Jastrow factor improves the uncorrelated trial wavefunction.

## Conclusion
We used the Variational Monte Carlo method on the system of two electrons in a 3-D quantum dot to estimate their ground state energy and mean separation at three values of oscillator frequency $\omega$. We used two trial wavefunctions, $\Psi_{T1}(\mathbf{R},\alpha)$ and $\Psi_{T2}(\mathbf{R},\alpha,\beta)$, where the second (correlated) function incorporates a Jastrow factor to account for the inter-electron repulsion.

For both trial functions, $E_{min}$ is lowest for $\omega=0.01$ since due to the low strength of confinement, the electrons will tend to wander farther from the origin. For $\omega=1.0$, since the strength of confinement is higher, the electrons will tend to gather closer to the origin, hence the higher energy.

Similarly, $\langle\mathbf{r_{12}}\rangle$ is highest for  $\omega=0.01$ since due to the low strength of confinement, the electrons will tend to wander outwards more and hence will tend to be spaced further apart. For $\omega=1.0$, the higher strength of confinement will bring the electrons closer to the origin, and subsequently, each other.

The lower and higher values of $E_{min}$ and $\langle\mathbf{r_{12}}\rangle$ respectively for the second trial function compared to the first is explained by the role of the Jastrow factor. Since it accounts for the inter-electron repulsion, it better models the system, and hence the VMC estimate of $E_{min}$ is lower. This inter-electron repulsion also leads to the electrons tending to stay comparatively further from each other, leading to an increased value of $\langle\mathbf{r_{12}}\rangle$ (especially when $\omega=0.01$, i.e. when the confinement is weak).

The reliability of the results for $\Psi_{T1}$ were calculated by running 100 loops of 5 million steps each. For the three values of $\omega$, the deviation from the analytical value of $E_{min}$ was less than $6.49\times 10^{-4}$, and the inter-run standard deviation for Energy and Separation was less than $6.37\times 10^{-4}$ and $2.06\times 10^{-2}$ respectively.

Similarly, 100 loops of 5 million steps each were ran for $\Psi_{T2}$. Since no suitable analytical expression exists for it, we calculated only the inter-run standard deviations for Energy and Separation. These were respectively less than $2.65\times 10^{-5}$ and $9.43\times 10^{-3}$.

#### An extension of this problem
This method can be used to find the shape and size of the orbitals of He atom, where in the Hamiltonian, the oscillator term $\frac{1}{2}\omega^2r_i^2$ will be replaced by the Coulomb potential between the electron and nucleus, i.e. $\frac{2}{r_i}$, and subsequently the trial functions and their corresponding local energies will change.

---
## Running the programme
The following is for Windows only. MacOS/Linux support is not guaranteed

#### Prerequisites: (preferrably, the latest versions of)
- Python
- GCC or any compatible C compiler
- Git (only for cloning the repository)

### Using uv (recommended):
If you have not done so, please install uv through the official windows instructions:

https://docs.astral.sh/uv/getting-started/installation/

Verify that the installation was successful before proceeding.

1. Clone the GitHub repository

```powershell
    git clone https://github.com/TheAtharva-Goyal/Quantum-Monte-Carlo-of-confined-electrons

    cd Quantum-Monte-Carlo-of-confined-electrons
```

2. Set up a Python virtual environment and activate it

```powershell
    uv sync
```

This creates/updates the project's .venv environment and installs the dependencies. Activate the virtual environment using

```powershell
    # for windows powershell
    .\.venv\Scripts\activate.ps1

    #for windows command prompt
    .\.venv\Scripts\activate.bat
```

3. Compile and run the C files

For example, using GCC: (ensure that the terminal is open in the folder containing the C files and that the venv environment is activated)

```powershell
    gcc cfile.c -O3 -march=native -o cfile
    if($?){ .\cfile }
```
4. To run the Python code independently

```powershell
    uv run pythonfile.py omega
```

Replace omega with your desired value of omega- 0.01, 0.5 or 1.0

(Note: the .csv files must exist in the same folder as the Python and C files)

5. Deactivate once you are done

```powershell
    deactivate
```

### Without using uv:
This method uses Python's in-built venv module and pip.

1. Clone the GitHub repository

```powershell
    git clone https://github.com/TheAtharva-Goyal/Quantum-Monte-Carlo-of-confined-electrons

    cd Quantum-Monte-Carlo-of-confined-electrons
```

2. Create a Python virtual environment

```powershell
    py -m venv .venv
```

Activate the virtual environment

```powershell
    # for windows powershell
    .\.venv\Scripts\activate.ps1

    #for windows command prompt
    .\.venv\Scripts\activate.bat
```

3. Install the dependencies

```powershell
    .\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

4. Compile and run the C files

For example, using GCC: (ensure that the terminal is open in the folder containing the C files and that the venv environment is activated)

```powershell
    gcc cfile.c -O3 -march=native -o cfile
    if($?){ .\cfile }
```

5. To run the Python code independently

```powershell
    python pythonfile.py omega
```

Replace omega with your desired value of omega- 0.01, 0.5 or 1.0

(Note: the .csv files must exist in the same folder as the Python and C files)

6. Deactivate once you are done

```powershell
    deactivate
```


---
### References
- Metropolis et al. (1953), Equation of State Calculations by Fast Computing Machines
- Hastings (1970), Monte Carlo Sampling Methods Using Markov Chains
- [M. Taut, Phys. Rev. A 48, 3561 (1993)](https://journals.aps.org/pra/pdf/10.1103/PhysRevA.48.3561)
- [Variational Principle short notes, Section 8.1](https://www2.ph.ed.ac.uk/~gja/qp/qp8.pdf)
- [Variational Monte Carlo wiki](https://en.wikipedia.org/wiki/Variational_Monte_Carlo)
- [Computational Physics GitHub repository](https://github.com/CompPhysics/ComputationalPhysics)

---
### Authors

[Amey Yadav](https://github.com/Amey-yadav26) -- [Atharva Goyal](https://github.com/TheAtharva-Goyal)

This is our submission for the final project for **Dr. Kirit Makwana**'s Computational Physics (EP4120) course in Jul-Nov '26, at IIT Hyderabad