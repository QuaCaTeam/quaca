# Power Spectrum

This abstract class serves as a container for different power spectra. So far only on type of power spectrum is implemented of the form
$$
\underline{S}(\omega) = \frac{\hbar}{\pi} \underline{\alpha}(\omega) \Big\{ \int \frac{d^2\mathbf{k}}{(2\pi)^2} \frac{\underline{G}_\Im(\mathbf{k},z_a,\mathbf{k}\cdot\mathbf{v}+\omega)}{1-e^{-\beta[\mathbf{k}\cdot\mathbf{v}+\omega]}} + \frac{1}{\alpha_0 \omega_a^2}\frac{\omega \mu_R(\omega)}{1-e^{-\beta \hbar \omega}} \Big\} \underline{\alpha}^\dagger(\omega)
$$

Here $\underline{\alpha}(\omega)$ represents the polarizability of the microscopic particle, $\mu(\omega)$ represents the internal bath of the microscopic particle and $\underline{G}(\mathbf{k},z_a,\omega)$ is the Green's tensor describing the geometry of the chosen set-up.
In the case, where the microscopic particle has no internal bath (e.g $\mu(\omega)=0$) the second term vanishes.

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

