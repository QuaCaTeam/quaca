# ReflectionCoefficients {docsify-ignore-all}
This is an abstract class that defines the reflection coefficients $r^{s/p}(\omega, \kappa)$ as a function of the frequency $\omega$ and $\kappa=\sqrt{k^2-\omega^2/c^2}$ (the following definition of the square root is used, $\mathrm{Re}\{\kappa\}\geq 0$ and $\mathrm{Im}\{\kappa\}<0$), where $k$ is the modulus of the two-dimensional wavevector $\mathbf{k}$. The specific reflection coefficients with respect to different systems are defined in the corresponding child of this class.
```cpp
class ReflectionCoefficients {
public:
  // returns the reflection coefficients
  virtual void calculate  (double omega, std::complex<double> kappa, std::complex<double> &r_p, std::complex<double> &r_s) = 0;
};

```

### `# virtual void calculate  (double omega, std::complex<double> kappa, std::complex<double> &r_p, std::complex<double> &r_s) = 0;`
Writes the complex reflection coefficients $r^{s/p}$ in the give variable.

## ReflectionCoefficientsLocBulk
Implements the reflection coefficient of a local bulk material according to
$$
r^s(\omega,\kappa) = \frac{\kappa -\kappa_\epsilon}{\kappa -\kappa_\epsilon}
\quad \text{and} \quad
r^p(\omega,\kappa) = \frac{\epsilon(\omega)\kappa -\kappa_\epsilon}{\epsilon(\omega)\kappa -\kappa_\epsilon}
$$
where $\epsilon(\omega)$ is the local permittivity of the bulk material and $\kappa_\epsilon=\sqrt{k^2-\epsilon(\omega)\omega^2/c^2}$ is defined analogously to $\kappa$. Here, a [Permittivity](api/permittivity) class has to be defined.

```cpp
class ReflectionCoefficientsLocBulk : public ReflectionCoefficients {
private:
  // permittivity is needed to describe the surface's response
  std::shared_ptr<Permittivity> permittivity;

public:
  /*!
   * Constructor for reflection coefficients of a local bulk medium.
   */
  ReflectionCoefficientsLocBulk(std::shared_ptr<Permittivity> permittivity);
  ReflectionCoefficientsLocBulk(const std::string &input_file);

  /*!
   * Returns the p- and s-polarized reflection coefficient.
   */
  void calculate(double omega, std::complex<double> kappa,
                 std::complex<double> &r_p, std::complex<double> &r_s) const;

  // getter functions
  std::complex<double> get_epsilon(double omega) const {
    return permittivity->calculate(omega);
  };
};
```


### `# ReflectionCoefficientsLocBulk(Permittivity *permittivity)`
Direct constructor for the class.

### `# ReflectionCoefficientsLocBulk(std::string input_file)`
Input file constructor for the class.

###  `#void calculate(double omega, std::complex<double> kappa,std::complex<double> &r_p, std::complex<double> &r_s) const;`
See [ReflectionCoefficients](#ReflectionCoefficients).

## ReflectionCoefficientsLocSlab
Implements the reflection coefficient of a plate of a local material and of finite thickness according to
$$
r^{s/p}(\omega,\kappa) = r^{s/p}_\mathrm{bulk}(\omega,\kappa) \frac{1-\exp(-2\kappa_\epsilon d)}{1-\left(r^{s/p}(\omega,\kappa) \exp(-\kappa_\epsilon d)\right)^2}
$$
where $r^{s/p}_\mathrm{bulk}$ refers to the corresponding bulk refelction coefficient, as defined in [ReflectionCoefficientsLocBulk](#ReflectionCoefficientsLocBulk), and $d$ is the thickness of the plate. 
```cpp
class ReflectionCoefficientsLocSlab : public ReflectionCoefficients {
private:
  // permittivity is needed to describe the surface's response
  std::shared_ptr<Permittivity> permittivity;
  double thickness;

public:
  /*!
   * Constructor for reflection coefficients of a local bulk medium.
   */
  ReflectionCoefficientsLocSlab(std::shared_ptr<Permittivity> permittivity,
                                double thickness);

  ReflectionCoefficientsLocSlab(const std::string &input_file);

  /*!
   * Returns the p- and s-polarized reflection coefficient.
   */
  void calculate(double omega, std::complex<double> kappa,
                 std::complex<double> &r_p, std::complex<double> &r_s) const;

  // getter functions
  std::complex<double> get_epsilon(double omega) const {
    return permittivity->calculate(omega);
  };

  double get_thickness() const { return thickness; };
};
```


### `# ReflectionCoefficientsLocSlab(Permittivity *permittivity, double thickness)`
Direct constructor for the class.

### `# ReflectionCoefficientsLocSlab(std::string input_file)`
Input file constructor for the class.

###  `#void calculate(double omega, std::complex<double> kappa,std::complex<double> &r_p, std::complex<double> &r_s) const;`
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
