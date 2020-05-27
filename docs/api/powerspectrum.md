# PowerSpectrum {docsify-ignore-all}

This abstract class serves as a container for the power spectrum of the harmonic dipole autocorrelator $\underline{S}$ as a functional of the [Green's tensor](api/greenstensor) $\underline{G}$, the [polarizability](api/polarizability) $\underline{\alpha}$ and the [memory kernel](api/memorykernel) $\mu$. The power spectrum of a harmonic oscillator can be defined with or without an internal bath. With an internal bath, the power spectrum reads
$$  \underline{S}(\omega) =\underline{\alpha}(\omega)\left\{ \int \frac{\mathrm{d}^2 \mathbf{k}}{(2 \pi)^2}\frac{ \underline{G}_\Im(\mathbf{k}, \omega + \mathbf{k}^\intercal\mathbf{v})}{1-\exp(-\hbar\omega(\omega+\mathbf{k}^\intercal\mathbf{v}))} + \frac{1}{\alpha_0\omega_a^2}\frac{\omega\mathrm{Re}\{\mu(\omega)\}}{1-\exp(-\hbar\beta\omega)}  \right\} \underline{\alpha}^\dagger(\omega), $$
and without an internal bath
$$  \underline{S}(\omega) =\underline{\alpha}(\omega)\int \frac{\mathrm{d}^2 \mathbf{k}}{(2 \pi)^2}\frac{ \underline{G}_\Im(\mathbf{k}, \omega + \mathbf{k}^\intercal\mathbf{v})}{1-\exp(-\hbar\omega(\omega+\mathbf{k}^\intercal\mathbf{v}))} \underline{\alpha}^\dagger(\omega), $$
where $\underline{G}$ is the Green's tensor described in [GreensTensor](api/greenstensor), $\underline{\alpha}$ the polarizability described in [Polarizability](api/polarizability), and $\mu$ the memory kernel described in [MemoryKernel](api/memorykernel).

```cpp
class PowerSpectrum {

protected:

  std::shared_ptr<GreensTensor>
      greens_tensor; // Green's tensor describing the geometry of
                     // the system
  std::shared_ptr<Polarizability>
      polarizability; // Polarizability describing the linear
                      // response of the microscopic particle

public:

  // Constructors
  PowerSpectrum(const std::string &input_file);

  // Constructor with initialization list
  PowerSpectrum(std::shared_ptr<GreensTensor> greens_tensor,
                std::shared_ptr<Polarizability> polarizability);

  // Calculate the power spectrum for a fixed value of the frequency \omega
  void calculate(double omega, cx_mat::fixed<3, 3> &powerspectrum,
                 Spectrum_Options spectrum) const;

  // getter functions
  std::shared_ptr<GreensTensor> &get_greens_tensor() { return greens_tensor; };
  std::shared_ptr<Polarizability> &get_polarizability() {
    return polarizability;
  };
};

```
### `PowerSpectrum(std::string input_file)`
Constructor with json-file

### `PowerSpectrum(GreensTensor *greens_tensor, Polarizability *polarizability)`
Direct constructor with initiatilization list

### `void calculate(double omega ,cx_mat::fixed<3, 3> &powerspectrum, Spectrum_Options spectrum) const`
Computes the power spectrum for a given frequency $\omega$. Here, the `Spectrum_Options spectrum` has two valid options. Either one calculates the full power spectrum with `spectrum = FULL`, as given in the equation above, or only the nonequilibrium part (for details see [this publication](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.117.100402)) with `spectrum = NON_LTE_ONLY`, which can be expressed mathemically as
$$  \frac{\hbar}{\pi}\underline{J}(\omega) = \underline{S}(\omega) - \frac{\hbar}{\pi}\frac{\underline{\alpha}_\Im(\omega)}{1-\exp(-\hbar\beta\omega)} $$

### `# get_...`
These are the getter functions of the respective quantity (`greens_tensor` or `polarizability`).

## Examples

<!-- tabs:start -->
#### **Example **
In this example we calculate the power spectrum at $\omega=0.03\,\mathrm{eV}$ of a $^{87}\mathrm{Rb}$ atom (without internal bath) with $\alpha_0 =47.3 \times (4\pi\epsilon_0)\text{\AA}^3 \approx 6\times 10^{-9}\,\mathrm{eV}^{-3}$ and $\omega_a=1.3\,\mathrm{eV}$ moving through vacuum at $v=10^{-5}c$ and $\beta=\frac{1}{300 \mathrm{K}} \approx 39 \,\mathrm{eV}^{-1}$ (for a conversion you can use [our converter](documentation/units)). Furthermore, we have to define an accuracy for the integration within the vacuum Green's tensor, which we choose $\delta_\mathrm{rel}^1=10^{-9}$. First, we construct the [Green's tensor](api/greenstensor) and the [polarizability](api/polarizability), and second, we build the power spectrum class of a harmonic oscillator.
```cpp
    // Constructing the Green's tensor
    double v = 1e-5;
    double beta = 39;
    double relerr = 1e-9;
    auto greens = std::make_shared<GreensTensorVacuum>(v, beta, relerr_k);
    
    // Constructing the polarizability
    double omega_a = 1.3;
    double alpha_zero = 6e-9;
    auto alpha = std::make_shared<Polarizability>(omega_a, alpha_zero, greens);

```
Now, we can construct and calculate the power spectrum as follows
```cpp
 // Construct power spectrum
 PowerSpectrum powerspectrum(greens, alpha);

 // Define a complex matrix
 cx_mat::fixed<3, 3> PS
  
 // Choose an omega
 double omega = 3e-2;

// Calculate power spectrum
 powerspectrum.calculate(omega, PS, FULL);
```
#### **Example (.json): Calculate power spectrum**
To construct the very same power spectrum of the previous tab via an input file we define ```parameters.json```
```json
{
    "Polarizability": {
        "omega_a": 1.3,
        "alpha_zero": 6e-9,
        }
    },
    "GreensTensor": {
        "type": "vacuum",
        "beta": 39,
        "v": 1e-5,
        "rel_err_1": 1E-9
    },
}
```
For construction and calculation we then simply type
```cpp
 // Construct power spectrum class
 PowerSpectrumHarmOsc powerspectrum("parameter.json");

 // Define a complex matrix
 cx_mat::fixed<3, 3> PS
  
 // Choose an omega
 double omega = 3e-2;

// Calculate power spectrum
 powerspectrum.calculate(omega, PS, FULL);

```
<!-- tabs:end -->
