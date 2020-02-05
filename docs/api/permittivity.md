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

<!-- tabs:end -->
