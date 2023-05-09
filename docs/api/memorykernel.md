# MemoryKernel {docsify-ignore-all}
This is an abstract class that defines the Fourier transform of the memory kernel, i.e. $\mu(\omega)$.
A specific model for this will be a child of this class.


## Member functions
### `std::complex<double> calculate(double omega)`
Returns the value $\mu(\omega)$ which in general is a complex number.
* Input parameters:
    * `double omega`: frequency, at which the memory kernel is evaluated.
* Return value:
    * `std::complex<double>`: complex number, representing the real and imaginary part of the memory kernel at the given frequency.

# OhmicMemoryKernel
Implements an ohmic bath memory kernel of the form
$$
\mu(\omega) = \gamma,
$$
where $\gamma$ is the damping coefficient of the bath.

## Member functions
### `OhmicMemoryKernel(double gamma)`
Direct constructor for the class.
* Input parameters:
    * `double gamma`: damping coefficent $\gamma$.
* Return value:
    * `OhmicMemoryKernel`: class instance.

### `OhmicMemoryKernel(std::string input_file)`
Input file constructor for the class.
* Input parameters:
    * `std::string input_file`: json-formatted file with all relevant quantities. See the final section of this page for an example.
* Return value:
    * `OhmicMemoryKernel`: class instance.

### `std::complex<double> calculate(double omega)`
See [MemoryKernel](#MemoryKernel).

### `double get_gamma()`
Getter function, which returns the constant $\gamma$

# SinglePhononMemoryKernel
Represents a memory kernel of an heat bath with one distinct (phonon) mode, as for example introduced in [this publication](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.98.155405). The memory kernel reads
$$
\mu(\omega) = \gamma +  \frac{g^2\omega_\mathrm{phon}^4}{\mathrm{i}\omega(\omega^2_\mathrm{phon} - \omega^2 - \mathrm{i}\gamma_\mathrm{phon} \omega)} ,
$$
where $\gamma$ is the ohmic damping, $g$ is the coupling constant of the distinct (phonon) mode to the dipole oscillation, $\omega_\mathrm{phon}$ is the frequency and $\gamma_\mathrm{phon}$ the damping of the distinct (phonon) mode. The respective header of this child class reads
## Member functions
### `SinglePhononMemoryKernel(double gamma, double gamma_phon, double omega_phon, double coupling);`
Direct constructor for the class.
* Input parameters:
    * `double gamma`: ohmic damping coefficent $\gamma$.
    * `double gamma_phon`: damping coefficent of the distinct mode $\gamma_{phon}$.
    * `double omega_phon`: frequency of the distinct mode $\omega_{phon}$.
    * `double coupling`: coupling constant $g$, which couples the distinct mode to the dipole oscialltion.
* Return value:
    * `SinglePhononMemoryKernel`: class instance.

### `SinglePhononKernel(std::string input_file)`
Input file constructor for the class.
* Input parameters:
    * `std::string input_file`: json-formatted file with all relevant quantities. See the final section of this page for an example.
* Return value:
    * `SinglePhononMemoryKernel`: class instance.

### `std::complex<double> mu(double omega)`
See [MemoryKernel](#MemoryKernel).

### `double get_...`
Getter functions, which return the respective quantity (`gamma`, `gamma_phon`, `omega_phon` or `coupling`)

## Input file
The input file sections for the memory kernels look like this

<!-- tabs:start -->
#### **OhmicMemoryKernel**
```json
{
  "MemoryKernel" : {
    "type" : "ohmic",
    "gamma" : 
  }
}
```
#### **SinglePhononMemoryKernel**
```json
{
    "MemoryKernel": {
        "type": "single_phonon",
        "gamma": ,
        "gamma_phon": ,
        "omega_phon": ,
        "coupling": 
    }
}
```
<!-- tabs:end -->

## Examples
The input file sections for the memory kernels look like this

<!-- tabs:start -->
#### **Example: Ohmice Memory Kernel**
If we want to construct an Ohmic memory kernel with $\gamma=0.1\,\mathrm{eV}$, we can employ the direct constructor in following form
```cpp
    double gamma = 0.1;
    OhmicMemoryKernel memorykernel(gamma);
    std::cout << memorykernel.get_gamma() << std::endl;
```
where we use the ```get_gamma()``` function to access the attribut of the class.

#### **Example (.json): Ohmice Memory Kernel**
Again, we can easily define a memory kernel with e.g. $\gamma=0.1\,\mathrm{eV}$ by employing a parameter file ```parameters.json```
```json
{
  "MemoryKernel" : {
    "type" : "ohmic",
    "gamma" : 0.1
  }
}
```
If we want access the memory kernel, we can use the ```get_gamma()``` function as follows
```cpp
    OhmicMemoryKernel memorykernel("parameters.json");
    std::cout << memorykernel.get_gamma() << std::endl;
```
#### **Example: Single Phonon Memory Kernel**
If we want to construct an single phonon memory kernel with $\gamma=0.1\,\mathrm{eV}$, $\gamma_\mathrm{phon}=10^{-5}\,\mathrm{eV}$, $\omega_\mathrm{phon}=4.34\,\mathrm{eV}$, and $g=10^{-5}$, we can employ the direct constructor in following form
```cpp
    double gamma = 0.1;
    double gamma_phon = 1e-5;
    double omega_phon = 4.34;
    double coupling = 1e-5;
    SinglePhononMemoryKernel memorykernel(gamma, gamma_phon, omega_phon, coupling);
    std::cout << memorykernel.get_gamma() << std::endl;
```
where we use the ```get_gamma()``` function to access the attribut of the class.

#### **Example (.json): Single Phonon Memory Kernel**
Again, we can easy define the singel phonon memory kernel given in the previous tab by employing a parameter file ```parameters.json```
```json
{
    "MemoryKernel": {
        "type": "single_phonon",
        "gamma": 1e-1,
        "gamma_phon": 1e-5,
        "omega_phon": 4.34,
        "coupling": 1E-5
    }
}
```
The construction of the single phonon memory kernel then simply reads
```cpp
    SinglePhononMemoryKernel memorykernel("parameters.json");
```
<!-- tabs:end -->
