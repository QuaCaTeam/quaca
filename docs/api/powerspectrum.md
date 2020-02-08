# Power Spectrum

This abstract class serves as a container for the power spectrum of the dipole autocorrelator.
The power spectrum defined as
$$  \underline{S}(\omega) =\underline{\alpha}(\omega)\left\{ \int \frac{\mathrm{d}^2 \mathbf{k}}{(2 \pi)^2}\frac{ \underline{G}_\Im(\mathbf{k}, \omega + \mathbf{k}^\intercal\mathbf{v})}{1-\exp(-\hbar\omega(\omega+\mathbf{k}^\intercal\mathbf{v}))} + \frac{1}{\alpha_0\omega_a^2}\frac{\omega\mathrm{Re}\{\mu(\omega)\}}{1-\exp(-\hbar\beta\omega)}  \right\} \underline{\alpha}(\omega)^\dagger, $$
where .

```cpp
class Polarizability {
protected:
  double omega_a;              // resonance frequency
  double alpha_zero;           // vacuum polarizability
  GreensTensor *greens_tensor; // Green's tensor
  std::string type;            // type of polarizability

public:
  // constructors
  Polarizability(double omega_a, double alpha_zero,
                 GreensTensor *greens_tensor);
  Polarizability(std::string input_file);

  // calculate the polarizability tensor
  virtual void calculate_tensor(cx_mat::fixed<3, 3> &alpha,
                                Options_Polarizability opts) = 0;

  // integration over omega
  double integrate_omega(Options_Polarizability opts, double omega_min,
                         double omega_max, double relerr, double abserr);
  static double integrand_omega(double omega, void *opts);

  // getter functions
  double get_omega_a() { return this->omega_a; };
  double get_alpha_zero() { return this->alpha_zero; };
};
```
