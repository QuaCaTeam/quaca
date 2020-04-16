# MemoryKernel {docsify-ignore-all}
This is an abstract class that defines the Fourier transform of the memory kernel, i.e. $\mu(\omega)$.
A specific model for this will be a child of this class.
```cpp
class MemoryKernel {
public:
  // Returns the memory kernel given a frequency omega.
  virtual std::complex<double> mu(double omega) = 0;
};
```

### `# std::complex<double> mu(double omega)`
Returns the value $\mu(\omega)$ which in general is a complex number.

## OhmicMemoryKernel
Implements an ohmic bath memory kernel of the form
$$
\mu(\omega) = \gamma,
$$
where $\gamma$ is the damping coefficient of the bath.

```cpp
class OhmicMemoryKernel : public MemoryKernel {
private:
  double gamma; // damping coefficient

public:
  // constructors
  OhmicMemoryKernel(double gamma);
  OhmicMemoryKernel(std::string input_file);

  // calculate function
  std::complex<double> mu(double omega);

  // getter functions
  double get_gamma() { return this->gamma; };
};
```

### `# OhmicMemoryKernel(double gamma)`
Direct constructor for the class.

### `# OhmicMemoryKernel(std::string input_file)`
Input file constructor for the class.

### `# std::complex<double> mu(double omega)`
See [MemoryKernel](#MemoryKernel).

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
<!-- tabs:end -->

## Examples
The input file sections for the memory kernels look like this

<!-- tabs:start -->
#### **Example: Ohmice Memory Kernel**
If we want to construct an Ohmic memory kernel with $\gamma=0.1\,\mathrm{eV}$, we can employ the direct constructor in following form
```cpp
    double gamma = 0.1;
    OhmicMemoryKernel memorykernel(gamma);
    std::cout << memorykernel.get_gamma();
```
where we use the ```get_gamma()``` function to access the attribut of the class.

#### **Example (.json): Ohmice Memory Kernel**
Again, we can easy define a memory kernel with e.g. $\gamma=0.1\,\mathrm{eV}$ by employing a parameter file ```parameters.json```
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
    std::cout << memorykernel.get_gamma();
```
<!-- tabs:end -->
