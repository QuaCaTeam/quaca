# Friction {docsify-ignore-all}
This abstract class contains the calculations for the noncontact friction. Here, two different general formulas are implemented
$$\begin{aligned} F_\mathrm{fric} =& -2 \int_{0}^{\infty} \mathrm{d}\omega\, \mathrm{Tr}\left[\underline{S}(\omega) \int\frac{\mathrm{d}^2\mathbf{k}}{(2\pi)^2} \, k_x\,\underline{G}_\Im(\mathbf{k}, z_a, k_xv+\omega)\right] \\ &+ 2\frac{\hbar}{\pi} \int_{0}^{\infty} \mathrm{d}\omega\,\mathrm{Tr}\left[\underline{\alpha}_\Im(\omega)\int\frac{\mathrm{d}^2\mathbf{k}}{(2\pi)^2} \, k_x\, \frac{\underline{G}_\Im(\mathbf{k}, z_a, k_xv+\omega)}{1-\exp(-\beta\hbar[\omega+k_xv])}\right],\end{aligned}$$
and
$$\begin{aligned} F_\mathrm{fric}=& -2\frac{\hbar}{\pi} \int_{0}^{\infty}\hspace{-.1cm} \mathrm{d}\omega\,\mathrm{Tr}\left[\underline{J}(\omega)\int\frac{\mathrm{d}^2\mathbf{k}}{(2\pi)^2} \, k_x\,\underline{G}_\Im(\mathbf{k}, z_a, k_xv+\omega)\right] \\&+2\frac{\hbar}{\pi}\int_{0}^{\infty}\hspace{-.1cm} \mathrm{d}\omega\,\mathrm{Tr}\left[\underline{\alpha}_\Im(\omega)\int\frac{\mathrm{d}^2\mathbf{k}}{(2\pi)^2} \, k_x\,\underline{G}_\Im(\mathbf{k}, z_a, k_xv+\omega)\right]\Theta_T(\omega+k_xv,\omega), \end{aligned}$$
where
$$ \Theta_T(\omega+k_xv,\omega)=\frac{1}{1-e^{-\beta\hbar[\omega+k_xv]}}-\frac{1}{1-e^{-\beta\hbar\omega}}. $$
Both expression yield analytically the same result, however, the second is computational-wise more stable. The power spectrum $\underline{S}$ and the nonequilibrium contribution $\underline{J}$ are described in [PowerSpectrum](api/powerspectrum), the Green's tensor $\underline{G}$ is described in [GreensTensor](api/greenstensor), and the polarizability $\underline{\alpha}$ is described in [Polarizability](api/polarizability). The header of the friction class reads
```cpp
class Friction {
protected:
  // the Green's tensor class
  std::shared_ptr<GreensTensor> greens_tensor;
  
  // the polarizability class
  std::shared_ptr<Polarizability> polarizability;
  
  // the power spectrum class
  std::shared_ptr<PowerSpectrum> powerspectrum;

  // relative accuracy of the omega integration
  double relerr_omega;

public:
  
  // constructor from an input file
  Friction(const std::string &input_file);

  // constructor from a direct input
  Friction(std::shared_ptr<GreensTensor> greens_tensor,
           std::shared_ptr<Polarizability> polarizability,
           std::shared_ptr<PowerSpectrum> powerspectrum, double relerr_omega);
  
  // calculation of the noncontact friction
  double calculate(Spectrum_Options spectrum) const;

  // the integrand of the noncontact friction
  double friction_integrand(double omega, Spectrum_Options spectrum) const;

  // getter functions
  std::shared_ptr<GreensTensor> &get_greens_tensor() { return greens_tensor; };
  std::shared_ptr<Polarizability> &get_polarizability() {
    return polarizability;
  };
  std::shared_ptr<PowerSpectrum> &get_powerspectrum() { return powerspectrum; };
};
```
### `Friction(const std::string &input_file);`
Constructor from a given `.json` file.

### `Friction(std::shared_ptr<GreensTensor> greens_tensor, std::shared_ptr<Polarizability> polarizability, std::shared_ptr<PowerSpectrum> powerspectrum, double relerr_omega);`
Direct constructor with initialization list,

### `double calculate(Spectrum_Options spectrum) const;`
Computes the noncontact friction. If `spectrum = FULL`, the full power spectrum is used, which tantamounts with the first equation in this section. If `spectrum = NON_LTE_ONLY`, solely the nonequilibrium part of the power spectrum is used, which represents the second equation. However, both expressions are analytically equal and solely differ in numerical efficiency/accuracy.

### `double friction_integrand(double omega, Spectrum_Options spectrum) const;`
The integrand wrapper for the $\omega$ integration.

### `# get_...`
These are the getter functions of the respective quantity (`greens_tensor`, `polarizability` or `powerspectrum`).

## Examples

<!-- tabs:start -->
#### **Example **
In this example we calculate the noncontact friction of a $^{87}\mathrm{Rb}$ atom (without internal bath) with $\alpha_0 =47.3 \times (4\pi\epsilon_0)\text{\AA}^3 \approx 6\times 10^{-9}\,\mathrm{eV}^{-3}$ and $\omega_a=1.3\,\mathrm{eV}$ moving through vacuum at $v=10^{-5}c$ and $\beta=\frac{1}{300 \mathrm{K}} \approx 39 \,\mathrm{eV}^{-1}$ (for a conversion you can use [our converter](documentation/units)). Furthermore, we have to define an accuracy for the integration within the vacuum Green's tensor, which we choose $\delta_\mathrm{rel}^1=10^{-9}$. First, we construct the [Green's tensor](api/greenstensor) the [polarizability](api/polarizability), and the power spectrum.
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

// Constructing the power spectrum
auto powerspectrum = std::make_shared<PowerSpectrum>(greens, alpha);
```
Now, we can construct the friction class and calculate the friction in the nonequilibrium form. Here, we define the numerical accuracy of the $\omega$ integration $\delta_{\omega}=10^{-6}$.

```cpp
Friction quant_fric(greens, alpha, powerspectrum, relerr_omega);
double num_result = quant_fric.calculate(NON_LTE_ONLY);
```
#### **Example (.json): Calculate power spectrum**
To construct the very same  of the friction of the previous tab via a input file we define ```parameters.json```
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
    "Friction" :
    {
        "relerr_omega" : 1e-6
    },
}
```
For construction and calculation we then simply type
```cpp
 // Construct friction class
 Friction quant_fric("../data/test_files/FrictionVacuum.json");
 double num_result = quant_fric.calculate(NON_LTE_ONLY);
```
<!-- tabs:end -->
