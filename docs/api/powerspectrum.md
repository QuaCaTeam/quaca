<<<<<<< HEAD
# PowerSpectrum

This abstract class serves as a container for the power spectrum of the dipole autocorrelator.
The power spectrum can be defined with or without an internal bath. Without an internal bath, the power spectrum reads
$$  \underline{S}(\omega) =\underline{\alpha}(\omega)\left\{ \int \frac{\mathrm{d}^2 \mathbf{k}}{(2 \pi)^2}\frac{ \underline{G}_\Im(\mathbf{k}, \omega + \mathbf{k}^\intercal\mathbf{v})}{1-\exp(-\hbar\omega(\omega+\mathbf{k}^\intercal\mathbf{v}))} + \frac{1}{\alpha_0\omega_a^2}\frac{\omega\mathrm{Re}\{\mu(\omega)\}}{1-\exp(-\hbar\beta\omega)}  \right\} \underline{\alpha}^\dagger(\omega), $$
and without an internal bath
$$  \underline{S}(\omega) =\underline{\alpha}(\omega)\int \frac{\mathrm{d}^2 \mathbf{k}}{(2 \pi)^2}\frac{ \underline{G}_\Im(\mathbf{k}, \omega + \mathbf{k}^\intercal\mathbf{v})}{1-\exp(-\hbar\omega(\omega+\mathbf{k}^\intercal\mathbf{v}))} \underline{\alpha}^\dagger(\omega), $$
where $\underline{G}$ is the Green's tensor described in [GreensTensor](api/greenstensor), $\underline{\alpha}$ the polarizability described in [Polarizability](api/polarizability), and $\mu$ the memory kernel described in [MemoryKernel](api/memorykernel).
=======
# Power Spectrum

This abstract class serves as a container for different power spectra. So far only on type of power spectrum is implemented of the form
$$
\underline{S}(\omega) = \frac{\hbar}{\pi} \underline{\alpha}(\omega) \Big\{ \int \frac{d^2\mathbf{k}}{(2\pi)^2} \frac{\underline{G}_\Im(\mathbf{k},z_a,\mathbf{k}\cdot\mathbf{v}+\omega)}{1-e^{-\beta[\mathbf{k}\cdot\mathbf{v}+\omega]}} + \frac{1}{\alpha_0 \omega_a^2}\frac{\omega \mu_R(\omega)}{1-e^{-\beta \hbar \omega}} \Big\} \underline{\alpha}^\dagger(\omega)
$$

Here $\underline{\alpha}(\omega)$ represents the polarizability of the microscopic particle, $\mu(\omega)$ represents the internal bath of the microscopic particle and $\underline{G}(\mathbf{k},z_a,\omega)$ is the Green's tensor describing the geometry of the chosen set-up.
In the case, where the microscopic particle has no internal bath (e.g $\mu(\omega)=0$) the second term vanishes.
>>>>>>> Simon_Beautifies_Code

```cpp
class PowerSpectrum {
public:
    
  //Green's tensor of describing the geometry of
  //the system
  GreensTensor *greens_tensor; 

  //Polarizability describing the linear of the microscopic particle
  Polarizability *polarizability; 

  //Constructor with ini-file
  PowerSpectrum(std::string input_file);

  //Constructor with initialization list
  PowerSpectrum(GreensTensor *greens_tensor, Polarizability *polarizability);

  //Calculate the power spectrum for a fixed value of the frequency \omega
  virtual void calculate(cx_mat::fixed<3, 3> &powerspectrum,
                         Options_PowerSpectrum opts) = 0;
};
```
### `PowerSpectrum(std::string input_file)`
Constructor with ini-file

### `PowerSpectrum(GreensTensor *greens_tensor, Polarizability *polarizability)`
Direct constructor with initiatilization list

### `virtual void calculate(cx_mat::fixed<3, 3> &powerspectrum,Options_PowerSpectrum opts) = 0`
Computes the power spectrum for a given frequency $\omega$

<<<<<<< HEAD
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

=======
>>>>>>> Simon_Beautifies_Code
