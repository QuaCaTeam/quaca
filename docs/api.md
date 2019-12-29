# API

QuaCa implements several physical quantities such as the polarizabilities and Green's tensors.


## Polarizability

### PolarizabilityBath
Implements a polarizability with an internal bath.
Obeys the equation
$ \underline{\alpha}(\omega) = \alpha_0 \left( \omega_a^2 - \omega - i \omega \mu(\omega) + \int \frac{d^2 \mathbf{k}}{(2 \pi)^2} \underline{G}(\mathbf{k}, \omega) \right)^{-1}$, where $\alpha_0$ is the vacuum polarizability, $\omega_a$ is the dipole resonance frequency, $\mu(\omega)$ is the memory kernel in Fourier space and $\underline{G}$ is the Green's tensor.


### PolarizabilityNoBath
Implements a polarizability with no internal bath.
Obeys the equation
$ \underline{\alpha}(\omega) = \alpha_0 \left( \omega_a^2 - \omega + \int \frac{d^2 \mathbf{k}}{(2 \pi)^2} \underline{G}(\mathbf{k}, \omega) \right)^{-1}$, where $\alpha_0$ is the vacuum polarizability, $\omega_a$ is the dipole resonance frequency and $\underline{G}$ is the Green's tensor.


### Examples
<!-- tabs:start -->

#### ** Example 1 **

We want to calculate the polarizability of a particle with an internal bath, above a plate with distance $z_a = 10$ and at frequency $3$.
First let's define all quantities that we need in order to construct the polarizability.
We assume that all relevant parameters are defined in a file called `parameters.ini`.
```cpp
OhmicMemoryKernel mu(parameters.ini);

PolarizabilityBath();
```


#### ** Example 2 **

Let's now calculate the polarizability of a particle with no internal bath, in the vacuum.

<!-- tabs:end -->





## Green's tensor

### GreensTensorVacuum
Implements the vacuum Green's tensor.
