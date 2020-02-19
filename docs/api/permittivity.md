!> TODO: Write examples. More physics in description of permittivity (local, isotropic, this kind of stuff).

## Permittivity
This is an abstract class that defines the permittivity as a function of the frequency, i.e. $\varepsilon(\omega)$.
A specific model for the permittivity will be a child of this class.
```cpp
class Permittivity {
public:
  // Return the permittivity given at frequency omega
  virtual std::complex<double> epsilon(double omega) = 0;
};
```

### `# std::complex<double> epsilon(double omega)`
Returns the value of $\varepsilon(\omega)$ which is a complex number.

## PermittivityDrude
Implements a Drude model according to the formula
$$
\varepsilon(\omega) = 1 - \frac{\omega_p^2}{\omega(\omega + i\gamma)},
$$
where $\omega_p$ is the plasma frequency and $\gamma$ is the damping coefficient.

```cpp
class PermittivityDrude : public Permittivity {
private:
  double omega_p, gamma; // plasma frequency and damping coefficient

public:
  // constructors
  PermittivityDrude(double omega_p, double gamma);
  PermittivityDrude(std::string input_file);

  // calculate the permittivity
  std::complex<double> epsilon(double omega);

  // getter functions
  double get_gamma() const { return this->gamma; };
  double get_omega_p() const { return this->omega_p; };
};
```


### `# PermittivityDrude(double omega_p, double gamma)`
Direct constructor for the class.

### `# PermittivityDrude(std::string input_file)`
Input file constructor for the class.

### `# std::complex<double> epsilon(double omega)`
See [Permittivity](#Permittivity).


## PermittivityLorentzNoBath
Implements a Lorentz model according to the formula
$$
\varepsilon(\omega) = \varepsilon_{\infty} - \frac{\alpha_0 \omega_0^2}{\omega_0^2 - \omega^2},
$$
where $\omega_0$ is the resonance frequency, $\varepsilon_{\infty}$ is the blabla and $\alpha_zero$ is the blabla.

```cpp
class PermittivityLorentzNoBath : public Permittivity {
private:
  double eps_inf;
  double alpha_zero;
  double omega_0;

public:
  // constructors
  PermittivityLorentzNoBath(double eps_inf, double alpha_zero, double omega_0);
  PermittivityLorentzNoBath(std::string input_file);

  // calculate the permittivity
  std::complex<double> epsilon(double omega);

  // Returns the numerical value of the permittivity scaled by omega.
  std::complex<double> epsilon_omega(double omega);

  // getter methods
  double get_eps_inf() { return this->eps_inf; };
  double get_alpha_zero() { return this->alpha_zero; };
  double get_omega_0() { return this->omega_0; };
};
```


### `# PermittivityLorentzNoBath(double eps_inf, double alpha_zero, double omega_0)`
Direct constructor for the class.

### `# PermittivityLorentzNoBath(std::string input_file)`
Input file constructor for the class.

### `# std::complex<double> epsilon(double omega)`
See [Permittivity](#Permittivity).



## PermittivityLorentz
Implements a Lorentz model with an internal bath according to the formula
$$
\varepsilon(\omega) = \varepsilon_{\infty} - \frac{\alpha_0 \omega_0^2}{\omega_0^2 - \omega^2 - \mathrm{i} \omega \mu(\omega)},
$$
where $\omega_0$ is the resonance frequency, $\varepsilon_{\infty}$ is the blabla, $\alpha_zero$ is the blabla and $\mu(\omega)$ is the memory kernel.

```cpp
class PermittivityLorentz : public Permittivity {
private:
  double eps_inf;
  double alpha_zero;
  double omega_0;

  MemoryKernel *memory_kernel;

public:
  // constructors
  PermittivityLorentz(double eps_inf, double alpha_zero, double omega_0,
                          MemoryKernel *memory_kernel);
  PermittivityLorentz(std::string input_file);

  // calculate the permittivity
  std::complex<double> epsilon(double omega);

  // Returns the numerical value of the permittivity scaled by omega.
  std::complex<double> epsilon_omega(double omega);

  // getter methods
  double get_eps_inf() { return this->eps_inf; };
  double get_alpha_zero() { return this->alpha_zero; };
  double get_omega_0() { return this->omega_0; };
  MemoryKernel *get_memory_kernel() { return this->memory_kernel; };
};
```


### `# PermittivityLorentz(double eps_inf, double alpha_zero, double omega_0, MemoryKernel *memory_kernel)`
Direct constructor for the class.

### `# PermittivityLorentz(std::string input_file)`
Input file constructor for the class.

### `# std::complex<double> epsilon(double omega)`
See [Permittivity](#Permittivity).


## Input file
The input file sections for the permittivities look like this

<!-- tabs:start -->
### **PermittivityDrude**
```ini
[Permittivity]
type = drude
omega_p =
gamma =
```

### **PermittivityLorentzNoBath**
```ini
[Permittivity]
type = lorentz nobath
eps_inf =
alpha_zero =
omega_0 =
```

### **PermittivityLorentz**
```ini
[Permittivity]
type = lorentz bath
eps_inf =
alpha_zero =
omega_0 =
```
For the Lorentz model with an internal bath you also need to define a [MemoryKernel](api/memorykernel)!
<!-- tabs:end -->
