# ReflectionCoefficients {docsify-ignore-all}
This is an abstract class that defines the reflection coefficients $r^{s/p}(\omega, \kappa)$ as a function of the frequency $\omega$ and $\kappa=\sqrt{k^2-\omega^2/c^2}$ (the following definition of the square root is used, $\mathrm{Re}\{\kappa\}\geq 0$ and $\mathrm{Im}\{\kappa\}<0$), where $k$ is the modulus of the two-dimensional wavevector $\mathbf{k}$. The specific reflection coefficients with respect to different systems are defined in the corresponding child of this class.

## Member functions
### `virtual void calculate  (double omega, std::complex<double> kappa, std::complex<double> &r_p, std::complex<double> &r_s) = 0;`
Writes the complex reflection coefficients $r^{s/p}$ in the give variable.
* Input parameters:
    * `double omega`: Frequency, at which the reflection coefficient is evaluated
    * `std::complex<double> kappa`: value of $\kappa$, at which the reflection coefficient is evaluated
    * `std::complex<double> &r_p$`: two dimensional vector of complex numbers, which stores the real and imaginary part of $r^p$
    * `std::complex<double> &r_s$`: two dimensional vector of complex numbers, which stores the real and imaginary part of $r^s$
* Return value: `void`

# ReflectionCoefficientsLocBulk
Implements the reflection coefficient of a local bulk material according to
$$
r^s(\omega,\kappa) = \frac{\kappa -\kappa_\epsilon}{\kappa -\kappa_\epsilon}
\quad \text{and} \quad
r^p(\omega,\kappa) = \frac{\epsilon(\omega)\kappa -\kappa_\epsilon}{\epsilon(\omega)\kappa -\kappa_\epsilon}
$$
where $\epsilon(\omega)$ is the local permittivity of the bulk material and $\kappa_\epsilon=\sqrt{k^2-\epsilon(\omega)\omega^2/c^2}$ is defined analogously to $\kappa$. Here, a [Permittivity](api/permittivity) class has to be defined.


## Member functions
### `ReflectionCoefficientsLocBulk(Permittivity *permittivity)`
Direct constructor for the class.
* Input parameters:
    * `Permittivity *permittivity`: Pointer to a class instant of `Permittivity`.
* Return value:
    * `ReflectionCoefficientLocBulk`: class instance.

### `ReflectionCoefficientsLocBulk(std::string input_file)`
Input file constructor for the class.
* Input parameters:
    * `std::string input_file`: json-formatted file with all relevant quantities. See the final section of this page for an example.
* Return value:
    * `ReflectionCoefficientLocBulk`: class instance.

###  `void calculate(double omega, std::complex<double> kappa,std::complex<double> &r_p, std::complex<double> &r_s) const;`
See [ReflectionCoefficients](#ReflectionCoefficients).

# ReflectionCoefficientsLocSlab
Implements the reflection coefficient of a plate of a local material and of finite thickness according to
$$
r^{s/p}(\omega,\kappa) = r^{s/p}_\mathrm{bulk}(\omega,\kappa) \frac{1-\exp(-2\kappa_\epsilon d)}{1-\left(r^{s/p}(\omega,\kappa) \exp(-\kappa_\epsilon d)\right)^2}
$$
where $r^{s/p}_\mathrm{bulk}$ refers to the corresponding bulk refelction coefficient, as defined in [ReflectionCoefficientsLocBulk](#ReflectionCoefficientsLocBulk), and $d$ is the thickness of the plate. 

## Member functions
### `ReflectionCoefficientsLocSlab(Permittivity *permittivity, double thickness)`
Direct constructor for the class.
* Input parameters:
    * `Permittivity *permittivity`: Pointer to a class instant of `Permittivity`.
* Return value:
    * `ReflectionCoefficientLocBulk`: class instance.

### `ReflectionCoefficientsLocSlab(std::string input_file)`
Input file constructor for the class.
* Input parameters:
    * `std::string input_file`: json-formatted file with all relevant quantities. See the final section of this page for an example.
* Return value:
    * `ReflectionCoefficientLocBulk`: class instance.

###  `void calculate(double omega, std::complex<double> kappa,std::complex<double> &r_p, std::complex<double> &r_s) const;`
See [ReflectionCoefficients](#ReflectionCoefficients).
## Input file
The input file sections for the permittivities look like this

<!-- tabs:start -->
### **ReflectionCoefficientsLocBulk**
```json
{
  "ReflectionCoefficients" : {
    "type" : "local bulk"
  }
}
```

### **ReflectionCoefficientsLocSlab**
```json
{
  "ReflectionCoefficients" : {
    "type" : "local slab",
    "thickness" : 0.05
  }
}
```
<!-- tabs:end -->

## Examples

<!-- tabs:start -->
### **Example: Local Bulk**
In order to construct and calculate $s$ and $p$ polarized reflection coefficient above a local Drude material at $\omega=1\,\mathrm{eV}$ and $\kappa = 10^{-3}\,\mathrm{eV}$, we first construct a permittivity with the chose parameters $\omega_p=9\,\mathrm{eV}$ and $\gamma=0.1\,\mathrm{eV}$.
```cpp
double gamma = 0.1;
double omega_p = 9;

// Define permittivity class
auto perm std::shared_ptr<PermittivityDrude>(omega_p,gamma);
```
Now we can construct and calculate the reflection coefficients at the desired $\omega$ and $\kappa$
```cpp
// Define reflection coefficient class
ReflectionCoefficientsLocBulk refl(perm); 

double omega =1e0;
std::complex<double> kappa = 1e-3;
std::complex<double> rp, rs;


refl.calculate(omega, kappa, rp, rs);
```

### **Example (.json): Local Bulk**
In order to construct and calculate $s$ and $p$ polarized reflection coefficient above a local Drude material at $\omega=1\,\mathrm{eV}$ and $\kappa = 10^{-3}\,\mathrm{eV}$, with a permittivity of the chose parameters $\omega_p=9\,\mathrm{eV}$ and $\gamma=0.1\,\mathrm{eV}$, we employ the following `parameter.json` file.
```json
{
    "ReflectionCoefficients": {
        "type": "local bulk",
    },
    "Permittivity": {
        "type": "drude",
        "gamma": 3.5e-2,
        "omega_p": 9
    }
}
```
Now we can construct and calculate the reflection coefficients at the desired $\omega$ and $\kappa$
```cpp
// Define reflection coefficient class
ReflectionCoefficientsLocBulk refl("parameter.json"); 

double omega =1e0;
std::complex<double> kappa = 1e-3;
std::complex<double> rp, rs;


refl.calculate(omega, kappa, rp, rs);
```
### **Example: Local Slab**
In order to construct and calculate $s$ and $p$ polarized reflection coefficient above a finite local Drude slab with a thickness of $d=10\,\mathrm{nm}=0.05\,\mathrm{eV}$ at $\omega=1\,\mathrm{eV}$ and $\kappa = 10^{-3}\,\mathrm{eV}$, we first construct a permittivity with the chosen parameters $\omega_p=9\,\mathrm{eV}$ and $\gamma=0.1\,\mathrm{eV}$.
```cpp
double gamma = 0.1;
double omega_p = 9;

// Define permittivity class
auto perm std::shared_ptr<PermittivityDrude>(omega_p,gamma);
```
Now we can construct and calculate the reflection coefficients at the desired $\omega$ and $\kappa$
```cpp
// Define reflection coefficient class
double thickness = 0.05;

ReflectionCoefficientsLocSlab refl(perm,thickness); 

double omega =1e0;
std::complex<double> kappa = 1e-3;
std::complex<double> rp, rs;


refl.calculate(omega, kappa, rp, rs);
```

### **Example (.json): Local Slab**
In order to construct and calculate $s$ and $p$ polarized reflection coefficient above a finite local Drude slab with a thickness of $d=10\,\mathrm{nm}=0.05\,\mathrm{eV}$ at $\omega=1\,\mathrm{eV}$ and $\kappa = 10^{-3}\,\mathrm{eV}$, with a permittivity of the chosen parameters $\omega_p=9\,\mathrm{eV}$ and $\gamma=0.1\,\mathrm{eV}$, we employ the following `parameter.json` file.
```json
{
    "ReflectionCoefficients": {
        "type": "local slab",
	"thickess": 0.05
    },
    "Permittivity": {
        "type": "drude",
        "gamma": 3.5e-2,
        "omega_p": 9
    }
}
```
Now we can construct and calculate the reflection coefficients at the desired $\omega$ and $\kappa$
```cpp
// Define reflection coefficient class
ReflectionCoefficientsLocSlab refl("parameter.json"); 

double omega =1e0;
std::complex<double> kappa = 1e-3;
std::complex<double> rp, rs;


refl.calculate(omega, kappa, rp, rs);
```
<!-- tabs:end -->
