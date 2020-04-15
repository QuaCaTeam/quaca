# PowerSpectrum {docsify-ignore-all}

This abstract class serves as a container for the power spectrum of the dipole autocorrelator $\underline{S}$ as a functional of the [Green's tensor](api/greenstensor) $\underline{G}$, the [polarizability](api/polarizability) $\underline{\alpha}$ and the [memory kernel](api/memorykernel) $\mu$.

```cpp
class PowerSpectrum {
public:

  //Green's tensor of describing the geometry of
  //the system
  GreensTensor *greens_tensor;

  //Polarizability describing the linear of the microscopic particle
  Polarizability *polarizability;

  //Constructor with json-file
  PowerSpectrum(std::string input_file);

  //Constructor with initialization list
  PowerSpectrum(GreensTensor *greens_tensor, Polarizability *polarizability);

  //Calculate the power spectrum for a fixed value of the frequency \omega
  virtual void calculate(cx_mat::fixed<3, 3> &powerspectrum,
                         Options_PowerSpectrum opts) = 0;
};
```
### `PowerSpectrum(std::string input_file)`
Constructor with json-file

### `PowerSpectrum(GreensTensor *greens_tensor, Polarizability *polarizability)`
Direct constructor with initiatilization list

### `virtual void calculate(cx_mat::fixed<3, 3> &powerspectrum,Options_PowerSpectrum opts) = 0`
Computes the power spectrum for a given frequency $\omega$

## Options_PowerSpectrum
This struct contains parameters and flags of the `PowerSpectrum`
```cpp
struct Options_PowerSpectrum {
  // Different integrand options:

  // Calculate the full power spectrum S
  bool full_spectrum = false;

  // Calculate solely the nonequilibrium contribution J
  bool non_LTE = false;

  // the frequency
  double omega = NAN;
};

```
### `bool full_spectrum`
If set to `true`, the power spectrum as shown above is calculated

### `bool non_LTE`
If set to `true`, only the nonequilibrium part of the power spectrum
$$  \frac{\hbar}{\pi}\underline{J}(\omega) = \underline{S}(\omega) - \frac{\hbar}{\pi}\frac{\underline{\alpha}_\Im(\omega)}{1-\exp(-\hbar\beta\omega)} $$
is calculated.
## PowerSpectrumHarmOsc
Is a child of the ```PowerSpectrum``` class and represents the power spectrum of a harmonic oscillator as already introduced in the beginning of this chapter.
The power spectrum of a harmonic oscillator can be defined with or without an internal bath. Without an internal bath, the power spectrum reads
$$  \underline{S}(\omega) =\underline{\alpha}(\omega)\left\{ \int \frac{\mathrm{d}^2 \mathbf{k}}{(2 \pi)^2}\frac{ \underline{G}_\Im(\mathbf{k}, \omega + \mathbf{k}^\intercal\mathbf{v})}{1-\exp(-\hbar\omega(\omega+\mathbf{k}^\intercal\mathbf{v}))} + \frac{1}{\alpha_0\omega_a^2}\frac{\omega\mathrm{Re}\{\mu(\omega)\}}{1-\exp(-\hbar\beta\omega)}  \right\} \underline{\alpha}^\dagger(\omega), $$
and without an internal bath
$$  \underline{S}(\omega) =\underline{\alpha}(\omega)\int \frac{\mathrm{d}^2 \mathbf{k}}{(2 \pi)^2}\frac{ \underline{G}_\Im(\mathbf{k}, \omega + \mathbf{k}^\intercal\mathbf{v})}{1-\exp(-\hbar\omega(\omega+\mathbf{k}^\intercal\mathbf{v}))} \underline{\alpha}^\dagger(\omega), $$
where $\underline{G}$ is the Green's tensor described in [GreensTensor](api/greenstensor), $\underline{\alpha}$ the polarizability described in [Polarizability](api/polarizability), and $\mu$ the memory kernel described in [MemoryKernel](api/memorykernel).
```cpp
class PowerSpectrumHarmOsc : public PowerSpectrum {
  public:

    // Constructor with initalization list
    PowerSpectrumHarmOsc(GreensTensor *greens_tensor, Polarizability *polarizability);

    // Constructor with json-file
    PowerSpectrumHarmOsc(std::string input_file);

    // calculate the power spectrum for a fixed value of the frequency
    void calculate(cx_mat::fixed<3, 3> &powerspectrum, Options_PowerSpectrum opts);

    // Flag to include or exclude an internal bath
    bool has_bath;
};
```

## Input
The input file section with respect to the power spectrum looks like this
```json
{
    "PowerSpectrum": {
        "type": "harmonic oscillator"
    }
}

```
## Examples
<!-- tabs:start -->
#### **Example: Calculate power spectrum**
In this example we want to calculate the power spectrum at $\omega=0.03\,\mathrm{eV}$ of a $^{87}\mathrm{Rb}$ atom (without internal bath) with $\alpha_0 =47.3 \times (4\pi\epsilon_0)\text{\AA}^3 \approx 6\times 10^{-9}\,\mathrm{eV}^{-3}$ and $\omega_a=1.3\,\mathrm{eV}$ moving through vacuum at $v=10^{-5}c$ and $\beta=\frac{1}{300 \mathrm{K}} \approx 39 \,\mathrm{eV}^{-1}$ (for a conversion you can use [our converter](documentation/units)). Furthermore, we have to define an accuracy for the integration within the vacuum Green's tensor, which we choose $\delta_\mathrm{rel}^1=10^{-9}$. First, we construct the [Green's tensor](api/greenstensor) and the [polarizability](api/polarizability), and second, we build the power spectrum class of a harmonic oscillator.
```cpp
    // Constructing the Green's tensor
    double v = 1e-5;
    double beta = 39;
    double relerr = 1e-9;

    GreensTensorVacuum greens_tensor(v, beta, relerr);
    
    // Constructing the polarizability
    double omega_a = 1.3;
    double alpha_zero = 6e-9;
    
    PolarizabilityNoBath pol(omega_a, alpha_zero, &greens_tensor);

```
Now, we can construct and calculate the power spectrum as follows
```cpp
 // Construct power spectrum
 PowerSpectrumHarmOsc powerspectrum(&greens, &alpha);

 // Define a complex matrix
 cx_mat::fixed<3, 3> PS
  
 // Choose an omega
 double omega = 3e-2;

 // Define options
  Options_PowerSpectrum opts;
  opts.spectrum = FULL;
  opts.omega = omega;

 powerspectrum.calculate(PS, opts);
```
#### **Example (.json): Calculate power spectrum**
To construct the very same power spectrum of the previous tab via a input file we define ```parameters.json```
```json
{
    "Polarizability": {
        "omega_a": 1.3,
        "type": "nobath",
        "alpha_zero": 6e-9,
        }
    },
    "GreensTensor": {
        "type": "vacuum",
        "beta": 39,
        "v": 1e-5,
        "rel_err_1": 1E-9
    },
    "PowerSpectrum": {
        "type": "harmonic oscillator"
    }
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

 // Define options
 Options_PowerSpectrum opts;
 opts.spectrum = FULL;
 opts.omega = omega;

 powerspectrum.calculate(PS, opts);

```
<!-- tabs:end -->
