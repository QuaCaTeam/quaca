# Polarizability {docsify-ignore-all}
This abstract class serves as a container for the  polarizability, with or without an internal heat bath.
As a physical quantity the polarizability reads as follows
$$  \underline{\alpha}(\omega) = \alpha_0 \omega_a^2 \left( \omega_a^2 - \omega^2 - \mathrm{i} \omega \mu(\omega) - \alpha_0 \omega_a^2 \int \frac{\mathrm{d}^2 \mathbf{k}}{(2 \pi)^2} \underline{G}(\mathbf{k}, \omega + \mathbf{k}^\intercal\mathbf{v}) \right)^{-1}, $$
with $\mu(\omega)$ being the [memory kernel](api/memorykernel) of the internal heat bath, $\underline{G}(\mathbf{k},\omega)$ the Green's tensor, $\alpha_0$ the static polarizability, and $\omega_a$ the resonance frequency of the microscopic object's dipole moment.
The header file of the polarizability reads
```cpp
class Polarizability {
private:
  double omega_a;                              // resonance frequency
  double alpha_zero;                           // vacuum polarizability
  std::shared_ptr<GreensTensor> greens_tensor; // Green's tensor
  std::shared_ptr<MemoryKernel> mu;            // internal bath

public:
  // Direct constructor without internal bath
  Polarizability(double omega_a, double alpha_zero,
                 std::shared_ptr<GreensTensor> greens_tensor);

  // Direct constructor with internal bath mu
  Polarizability(double omega_a, double alpha_zero,
                 std::shared_ptr<MemoryKernel> mu,
                 std::shared_ptr<GreensTensor> greens_tensor);

  // Constructor from a given json file
  explicit Polarizability(const std::string &input_file);

  // calculate the polarizability tensor
  void calculate_tensor(double omega, cx_mat::fixed<3, 3> &alpha,
                        Tensor_Options fancy_complex) const;

  // integration over omega
  double integrate_omega(const vec::fixed<2> &indices,
                         Tensor_Options fancy_complex, double omega_min,
                         double omega_max, double relerr, double abserr) const;

  // integrand for the omega integration
  double integrand_omega(double omega, const vec::fixed<2> &indices,
                         Tensor_Options fancy_complex) const;

  // getter functions
  double get_omega_a() const { return omega_a; };
  double get_alpha_zero() const { return alpha_zero; };
  const std::shared_ptr<GreensTensor> &get_greens_tensor() const {
    return greens_tensor;
  };
  // getter function for memory kernel
  std::complex<double> get_mu(double omega) const {
    if (mu != nullptr) {
      return mu->calculate(omega);
    } else {
      return std::complex<double>(0.0, 0.0);
    }
  };
};
```
and contains following objects:

### `# Polarizability(double omega_a, double alpha_zero, std::shared_ptr<GreensTensor> greens_tensor);`
Direct constructor for the class without an internal bath.

### `# Polarizability(double omega_a, double alpha_zero, std::shared_ptr<MemoryKernel> mu, std::shared_ptr<GreensTensor> greens_tensor);`
Direct constructor for the class with an internal bath.

### `# Polarizability(const std::string &input_file)`
Input file constructor for the class.

### `# void calculate_tensor(double omega, cx_mat::fixed<3, 3> &alpha, Tensor_Options fancy_complex) const;`

This function evaluates the tensor according to its formula for $\underline{\alpha}(\omega)$ for a specific `omega` and puts the result into the matrix `alpha`. If `fancy_complex = IM` then $\underline{\alpha}_{\Im}= (\underline{\alpha} -\underline{\alpha}^\dagger)/(2\mathrm{i})$ and for `fancy_complex = RE`  $\underline{\alpha}_{\Re} = (\underline{\alpha} +\underline{\alpha}^\dagger)/2$ is calculated. See also the [Examples](#Examples).

### `# double integrate_omega(const vec::fixed<2> &indices, Tensor_Options fancy_complex, double omega_min, double omega_max, double relerr, double abserr) const;`
Integrates the polarizability element $\alpha_{ij}$ with the `indices` $i$ and $j$ from `omega_min` to `omega_max` with a relative error `relerr` and absolute error `abserr` with the [CQUAD](https://www.gnu.org/software/gsl/doc/html/integration.html) integration routine.
If `fancy_complex = IM` then $\underline{\alpha}_{\Im}= (\underline{\alpha} -\underline{\alpha}^\dagger)/(2\mathrm{i})$ and for `fancy_complex = RE`  $\underline{\alpha}_{\Re} = (\underline{\alpha} +\underline{\alpha}^\dagger)/2$ is calculated. See also Example 2 in [Examples](#Examples).

### `# get_...`
These are the getter functions of the respective quantity (`omega_a`, `alpha_zero`, `greens_tensor` or `mu`).

## Input file
The input file sections for the polarizabilities look like the following.
Do not forget that for all polarizabilities you need to define a [GreensTensor](api/greenstensor)!
<!-- tabs:start -->
#### **PolarizabilityNoBath**
```json
{
    "Polarizability": {
        "omega_a": ,
        "alpha_zero": 
    }
}
```


#### **PolarizabilityBath**
```json
{
    "Polarizability": {
        "omega_a": ,
        "alpha_zero": ,
        "MemoryKernel": {
        }
    }
}
```
For the bath you also need to define a [MemoryKernel](api/memorykernel)!
<!-- tabs:end -->


## Examples
<!-- tabs:start -->

#### ** Example **

We want to calculate the polarizability of a particle with an ohmic internal bath, above a plate with distance $z_a = 10\,\mathrm{eV}^{-1}$, velocity $v=0.1\,c$ and at frequency $3\,\mathrm{eV}$.
First let's define the ohmic memory kernel
```cpp
double gamma = 3.0;
OhmicMemoryKernel mu(gamma);
```
Furthermore, we need to generate a [permittivity](api/permittivity) and the [reflection coefficients](api/reflection), which represent characteristic properties of the surface, where we chose $\omega_p=9\,\mathrm{eV}$ and $\gamma=0.1\,\mathrm{eV}$
```cpp
    double omega_p = 9;
    double gamma = 0.1;
    auto perm = std::make_shared<PermittivityDrude>(omega_p, gamma);
    auto refl = std::make_shared<ReflectionCoefficientsLocBulk>(perm);
```
Now we can define the [Green's tensor](api/greenstensor) (including the vacuum as well as the scattered contribution) 
```cpp
double v = 0.1;
double z_a = 10;
double delta_cut = 30;
vec::fixed<2> rel_err = {1E-8, 1E-6};

auto greens = std::make_shared<GreensTensorPlateVacuum>(v, beta, z_a, refl, delta_cut, relerr_k);
```
Here, we introduced further the inverse temperature $\beta$ (which is, however, irrelevant for the current calculation), the relative accuracy `rel_err` of the integration over the Green's tensor, and the $\delta_\mathrm{cut}$, which defines a numerical cut-off of the aformentioned integration.
Finally we can put all together and define the polarizability
```cpp
double omega_a = 1.3;
double alpha_zero = 4.0;
Polarizability pol(omega_a, alpha_zero, mu, greens_tensor);
```
To calculate $\underline{\alpha}(3.0\,\mathrm{eV})$ we then do
```cpp
cx_mat::fixed<3,3> alpha(fill::zeros);
double omega = 3.0;
pol.calculate_tensor(omega, alpha, IM);
```
The matrix `alpha` now contains the correctly calculated imaginary part of the polarizability tensor $\alpha_\Im(\omega)$.


#### ** Example (.json) **
Since we have to define a lot of parameters, QuaCa offers a shortcut to the long task before.
We simply define all parameters in a file called `parameters.json` which looks like this
```json
{
    "Polarizability": {
        "omega_a": 1.3,
        "alpha_zero": 4.0,
        "MemoryKernel": {
            "type": "ohmic",
            "gamma": 3.0
        }
    },
    "GreensTensor": {
        "type": "plate",
        "v": 0.1,
        "za" : 10,
        "beta": 3.2
    }
}
```
Now we can easily define and calculate the imaginary part of the polarizability as such
```cpp
PolarizabilityBath pol("parameters.json");

cx_mat::fixed<3,3> alpha(fill::zeros);
pol.calculate_tensor(3.0, alpha, IM);
```



<!-- tabs:end -->
