# Permittivity {docsify-ignore-all}
This is an abstract class that defines the permittivity as a function of the frequency, i.e. $\varepsilon(\omega)$.
A specific model for the permittivity will be a child of this class.
```cpp
class Permittivity {
public:
  // calculate the permittivity
  virtual std::complex<double> calculate(double omega) const = 0;
};
```

### `# std::complex<double> calculate(double omega)`
Returns the value of $\varepsilon(\omega)$ which is a complex number.

## PermittivityDrude
Implements a Drude model according to the formula
$$
\varepsilon(\omega) = 1 - \frac{\omega_p^2}{\omega(\omega + \mathrm{i}\gamma)},
$$
where $\omega_p$ is the plasma frequency and $\gamma$ is the damping coefficient.

```cpp
class PermittivityDrude : public Permittivity {
private:
  double omega_p; // plasma frequency
  double gamma;   // damping coefficient

public:
  // constructors
  PermittivityDrude(double omega_p, double gamma);
  explicit PermittivityDrude(const std::string &input_file);

  // calculate the permittivity
  std::complex<double> calculate(double omega) const override;

  // getter methods
  double get_gamma() const { return this->gamma; };
  double get_omega_p() const { return this->omega_p; };
};
```


### `# PermittivityDrude(double omega_p, double gamma)`
Direct constructor for the class.

### `# PermittivityDrude(std::string input_file)`
Input file constructor for the class.

### `# std::complex<double> calculate(double omega)`
See [Permittivity](#Permittivity).

### `# double get_...`
Getter functions of the respective quantity (`gamma` or `omega_p`)

## PermittivityLorentz
Implements a Lorentz model with an internal bath according to the formula
$$
\varepsilon(\omega) = \varepsilon_{\infty} - \frac{\omega_p^2}{\omega_0^2 - \omega^2 - \mathrm{i} \omega \mu(\omega)},
$$
where $\omega_0$ is the central frequency, $\varepsilon_{\infty}$ is high-frequency permittivity limit, $\omega_p$ is the plasma frequency and $\mu(\omega)$ is the memory kernel.

```cpp
class PermittivityLorentz : public Permittivity {
private:
  double eps_inf;
  double omega_p;
  double omega_0;

  std::shared_ptr<MemoryKernel> memory_kernel;

public:
  // constructors
  PermittivityLorentz(double eps_inf, double omega_p, double omega_0,
                      std::shared_ptr<MemoryKernel> memory_kernel);
  explicit PermittivityLorentz(const std::string &input_file);

  // calculate the permittivity
  std::complex<double> calculate(double omega) const override;

  // getter methods
  double get_eps_inf() const { return this->eps_inf; };
  double get_omega_p() const { return this->omega_p; };
  double get_omega_0() const { return this->omega_0; };
  std::shared_ptr<MemoryKernel> &get_memory_kernel() {
    return this->memory_kernel;
  };
};
```


### `# PermittivityLorentz(double eps_inf, double omega_p, double omega_0, std::shared_ptr<MemoryKernel> memory_kernel)`
Direct input constructor for the class.

### `# PermittivityLorentz(std::string input_file)`
Input file constructor for the class.

### `# std::complex<double> calculate(double omega)`
See [Permittivity](#Permittivity).

### `# get_...`
Getter functions of the respective quantity (`gamma`, `omega_p`, `eps_inf` or `memory_kernel`)


## Input file
The input file sections for the permittivities look like this

<!-- tabs:start -->
### **PermittivityDrude**
```json
{
    "Permittivity": {
        "type": "drude",
        "gamma": ,
        "omega_p": 
    }
}
```

### **PermittivityLorentz**
```json
{
    "Permittivity": {
        "type": "lorentz",
        "eps_inf": ,
        "omega_p": ,
        "omega_0": ,
        "MemoryKernel": {
        }
    }
}
```
For the Lorentz model with an internal bath you also need to define a [MemoryKernel](api/memorykernel)!
<!-- tabs:end -->

## Examples

<!-- tabs:start -->
### **Example: Drude**
In order to construct and calculate the permittivity of a Drude material with $\omega_p=9\,\mathrm{eV}$ and $\gamma=0.1\,\mathrm{eV}$ at $\omega=1\,\mathrm{eV}$, we can employ direct constructor as follows
```cpp
double gamma = 0.1;
double omega_p = 9;

// Define permittivity class
PermittivityDrude perm(omega_p, gamma);

// Calculate permittivity
double omega = 1;
std::complex<double> eps = perm.calculate(omega) 
```

### **Example: Drude (.json)**
If we want to construct and calculate the same Drude permittivity as in the previous tab, but with the following `parameters.json` file
```json
{
    "Permittivity": {
        "type": "drude",
        "gamma": 0.1 ,
        "omega_p": 9.
    }
}
```
we can simply write
```cpp
// Define permittivity class
PermittivityDrude perm("parameters.json");

// Calculate permittivity
double omega = 1;
std::complex<double> eps = perm.calculate(omega) 
```
### **Example: Lorentz**
In order to construct and calculate the permittivity of a Lorentz material with $\omega_p=0.6\,\mathrm{eV}$, $\omega_0=3.4\,\mathrm{eV}$, $\epsilon_\infty=1.4$ and $\mu(\omega)=\gamma=0.69420\,\mathrm{eV}$ at $\omega=1\,\mathrm{eV}$, we can employ direct constructor as follows
```cpp
// Construct the memory kernel
double gamma = 0.69420;
auto mu = std::make_shared<OhmicMemoryKernel>(gamma);

// Define permittivity class

double eps_inf = 1.4;
double omega_p = 6e-1;
double omega_0 = 3.4;

PermittivityLorentz perm(eps_inf, omega_p, omega_0, mu);

// Calculate permittivity
double omega = 1;
std::complex<double> eps = perm.calculate(omega) 
```

### **Example: Lorentz (.json)**
If we want to construct and calculate the same Drude permittivity as in the previous tab, but with the following `parameters.json` file
```json
{
    "Permittivity": {
        "type": "lorentz",
        "eps_inf": 1.4,
        "omega_p": 6e-1,
        "omega_0": 3.4,
        "MemoryKernel": {
            "type": "ohmic",
            "gamma": 0.69420
        }
    }
}
```
we can simply write
```cpp
// Define permittivity class
PermittivityLorentz perm("parameters.json");

// Calculate permittivity
double omega = 1;
std::complex<double> eps = perm.calculate(omega) 
```
<!-- tabs:end -->
