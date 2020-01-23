?> We try to use the same variable names for the same physical quantities, which resemble the LaTeX writing. Please have a look at the [reference table](documentation/units)

# Polarizability
This abstract class serves as a container for two different polarizabilities: one of a particle with an internal heat bath and one without such a bath.
We expect the polarizability to be of the form
$$  \underline{\alpha}(\omega) = \alpha_0 \omega_a^2 \left( \omega_a^2 - \omega^2 - i \omega \mu(\omega) - \alpha_0 \omega_a^2 \int \frac{d^2 \mathbf{k}}{(2 \pi)^2} \underline{G}(\mathbf{k}, \omega + \mathbf{k}^{\mathrm{T}}\mathbf{v}) \right)^{-1}, $$
where we set $\mu(\omega) = 0$ for the particle without an internal heat bath.

```cpp
class Polarizability {
protected:
  double omega_a;              // resonance frequency
  double alpha_zero;           // vacuum polarizability
  GreensTensor *greens_tensor; // Green's tensor

public:
  // constructors
  Polarizability(double omega_a, double alpha_zero,
                 GreensTensor *greens_tensor);
  Polarizability(std::string input_file);

  // calculate the polarizability tensor
  virtual void calculate_tensor(cx_mat::fixed<3, 3> &alpha,
                                Options_Polarizability opts) = 0;

  // integration over omega
  double integrate_omega(Options_Polarizability opts, double omega_min,
                         double omega_max, double relerr, double abserr);
  static double integrand_omega(double omega, void *opts);

  // getter functions
  double get_omega_a();
  double get_alpha_zero();
};
```

# Options_Polarizability
This struct is given to the functions in the `Polarizability` class to specify what is being calculated.

```cpp
struct Options_Polarizability {
  // Different options for the integrand
  bool fancy_R = false;
  bool fancy_I = false;

  double omega = NAN;

  // Indices of the 3x3 polarizability tensor
  vec::fixed<2> indices = {NAN, NAN};

  // Pointer to the Polarizability to be able to access the attributes of the
  // class eventhough the integrand is static
  Polarizability *class_pt;
};
```


# PolarizabilityBath
Implements a polarizability with an internal bath.
Obeys the equation
$  \underline{\alpha}(\omega) = \alpha_0 \omega_a^2 \left( \omega_a^2 - \omega^2 - i \omega \mu(\omega) - \alpha_0 \omega_a^2 \int \frac{d^2 \mathbf{k}}{(2 \pi)^2} \underline{G}(\mathbf{k}, \omega + \mathbf{k}^{\mathrm{T}}\mathbf{v}) \right)^{-1}$,
where $\alpha_0$ is the vacuum polarizability, $\omega_a$ is the dipole resonance frequency, $\mu(\omega)$ is the memory kernel in Fourier space and $\underline{G}$ is the Green's tensor.

```cpp
class PolarizabilityBath : public Polarizability
{
private:

    MemoryKernel *mu; // the memory kernel

public:

    // constructors
    PolarizabilityBath(double omega_a, double alpha_zero, MemoryKernel *mu, GreensTensor *greens_tensor);
    PolarizabilityBath(std::string input_file);

    // calculate the polarizability tensor
    void calculate_tensor(cx_mat::fixed<3,3>& alpha, double omega);

    // getter function for mu
    std::complex<double> get_mu(double omega);
};
```

# PolarizabilityNoBath
Implements a polarizability with no internal bath.
Obeys the equation
$  \underline{\alpha}(\omega) = \alpha_0 \omega_a^2 \left( \omega_a^2 - \omega^2 - \alpha_0 \omega_a^2 \int \frac{d^2 \mathbf{k}}{(2 \pi)^2} \underline{G}(\mathbf{k}, \omega + \mathbf{k}^{\mathrm{T}}\mathbf{v}) \right)^{-1}$,
where $\alpha_0$ is the vacuum polarizability, $\omega_a$ is the dipole resonance frequency and $\underline{G}$ is the Green's tensor.
```cpp
class PolarizabilityNoBath : public Polarizability
{
public:

  // constructors
  PolarizabilityNoBath(double omega_a, double alpha_zero, GreensTensor *greens_tensor);
  PolarizabilityNoBath(std::string input_file);

  // calculate the polarizability tensor
  void calculate_tensor(cx_mat::fixed<3,3>& alpha, double omega);

};
```
# Examples
<!-- tabs:start -->

#### ** Example 1 **

We want to calculate the polarizability of a particle with an ohmic internal bath, above a plate with distance $z_a = 10$ and at frequency $3$.
First let's define the ohmic memory kernel
```cpp
double gamma = 3.0;
OhmicMemoryKernel mu(gamma);
```
Now let's define the Green's tensor
```cpp
double v = 0.1;
double z_a = 10;
double beta = 1e4;
GreensTensorPlate greens_tensor(v, z_a, beta);
```
Finally we can put all together and define the polarizability
```cpp
double omega_a = 1.3;
double alpha_zero = 4.0;
PolarizabilityBath pol(omega_a, alpha_zero, &mu, &greens_tensor);
```
To calculate $\underline{\alpha}(3.0)$ we then do
```cpp
cx_mat::fixed<3,3> alpha(fill::zeros);
pol.calculate_tensor(alpha, 3.0);
```
The matrix `alpha` now contains the correctly calculated polarizability tensor.


#### ** Example 1 (easier) **
Since we have to define a lot of parameters, QuaCa offers a shortcut to the long task before.
We simply define all parameters in a file called `parameters.ini` which looks like this
```ini
[MemoryKernel]
type = ohmic
gamma = 3.0

[GreensTensor]
type = plate
v = 0.1
za = 10
beta = 1e4

[Polarizability]
type = bath
omega_a = 1.3
alpha_zero = 4.0
```
Now we can easily define and calculate the polarizability as such
```cpp
PolarizabilityBath pol("parameters.ini");

cx_mat::fixed<3,3> alpha(fill::zeros);
pol.calculate_tensor(alpha, 3.0);
```



<!-- tabs:end -->
