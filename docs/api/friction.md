# Friction {docsify-ignore-all}
This abstract class contains the calculations for the noncontact friction. Here, two different general formulas are implemented
$$\begin{aligned} F_\mathrm{fric} =& -2 \int_{0}^{\infty} \mathrm{d}\omega\, \mathrm{Tr}\left[\underline{S}(\omega) \int\frac{\mathrm{d}^2\mathbf{k}}{(2\pi)^2} \, k_x\,\underline{G}_\Im(\mathbf{k}, z_a, k_xv+\omega)\right] \\ &+ 2\frac{\hbar}{\pi} \int_{0}^{\infty} \mathrm{d}\omega\,\mathrm{Tr}\left[\underline{\alpha}_\Im(\omega)\int\frac{\mathrm{d}^2\mathbf{k}}{(2\pi)^2} \, k_x\, \frac{\underline{G}_\Im(\mathbf{k}, z_a, k_xv+\omega)}{1-\exp(-\beta\hbar[\omega+k_xv])}\right],\end{aligned}$$
and
$$\begin{aligned} F_\mathrm{fric}=& -2\frac{\hbar}{\pi} \int_{0}^{\infty}\hspace{-.1cm} \mathrm{d}\omega\,\mathrm{Tr}\left[\underline{J}(\omega)\int\frac{\mathrm{d}^2\mathbf{k}}{(2\pi)^2} \, k_x\,\underline{G}_\Im(\mathbf{k}, z_a, k_xv+\omega)\right] \\&+2\frac{\hbar}{\pi}\int_{0}^{\infty}\hspace{-.1cm} \mathrm{d}\omega\,\mathrm{Tr}\left[\underline{\alpha}_\Im(\omega)\int\frac{\mathrm{d}^2\mathbf{k}}{(2\pi)^2} \, k_x\,\underline{G}_\Im(\mathbf{k}, z_a, k_xv+\omega)\right]\Theta_T(\omega+k_xv,\omega), \end{aligned}$$
where
$$ \Theta_T(\omega+k_xv,\omega)=\frac{1}{1-e^{-\beta\hbar[\omega+k_xv]}}-\frac{1}{1-e^{-\beta\hbar\omega}}. $$
Both expression are equally valid for the calculation, however, the second is computational-wise more stable. The power spectrum $\underline{S}$ and the nonequilibrium contribution $\underline{J}$ are described in [PowerSpectrum](api/powerspectrum), the Green's tensor $\underline{G}$ is described in [GreensTensor](api/greenstensor), and the polarizability $\underline{\alpha}$ is described in [Polarizability](api/polarizability).
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

