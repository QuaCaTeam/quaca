## Friction
This abstract class contains the calculations for the noncontact friction. Here, two different general formulas are implemented
$$\begin{aligned} F_\mathrm{fric} =& -2 \int_{0}^{\infty} \mathrm{d}\omega\, \mathrm{Tr}\left[\underline{S}(\omega) \int\frac{\mathrm{d}^2\mathbf{k}}{(2\pi)^2} \, k_x\,\underline{G}_\Im(\mathbf{k}, z_a, k_xv+\omega)\right] \\ &+ 2\frac{\hbar}{\pi} \int_{0}^{\infty} \mathrm{d}\omega\,\mathrm{Tr}\left[\underline{\alpha}_\Im(\omega)\int\frac{\mathrm{d}^2\mathbf{k}}{(2\pi)^2} \, k_x\, \frac{\underline{G}_\Im(\mathbf{k}, z_a, k_xv+\omega)}{1-\exp(-\beta\hbar[\omega+k_xv])}\right],\end{aligned}$$
and
$$\begin{aligned} F_\mathrm{fric}=& -2\frac{\hbar}{\pi} \int_{0}^{\infty}\hspace{-.1cm} \mathrm{d}\omega\,\mathrm{Tr}\left[\underline{J}(\omega)\int\frac{\mathrm{d}^2\mathbf{k}}{(2\pi)^2} \, k_x\,\underline{G}_\Im(\mathbf{k}, z_a, k_xv+\omega)\right] \\&+2\frac{\hbar}{\pi}\int_{0}^{\infty}\hspace{-.1cm} \mathrm{d}\omega\,\mathrm{Tr}\left[\underline{\alpha}_\Im(\omega)\int\frac{\mathrm{d}^2\mathbf{k}}{(2\pi)^2} \, k_x\,\underline{G}_\Im(\mathbf{k}, z_a, k_xv+\omega)\right]\Theta_T(\omega+k_xv,\omega), \end{aligned}$$
where
$$ \Theta_T(\omega+k_xv,\omega)=\frac{1}{1-e^{-\beta\hbar[\omega+k_xv]}}-\frac{1}{1-e^{-\beta\hbar\omega}}. $$
Both expression are equally valid for the calculation, however, the second is computational-wise more stable. The power spectrum $\underline{S}$ and the nonequilibrium contribution $\underline{J}$ are described in [PowerSpectrum](api/powerspectrum), the Green's tensor $\underline{G}$ is described in [GreensTensor](api/greenstensor), and the polarizability $\underline{\alpha}$ is described in [Polarizability](api/polarizability).
```cpp
class Friction {
public:
  // defining the demanded relative accuracy of the omega integration
  double relerr_omega;

  // the Green's tensor class
  GreensTensor *greens_tensor;

  // the polarizability class
  Polarizability *polarizability;

  // the power spectrum class
  PowerSpectrum *powerspectrum;
  
  // contructor from a input file
  Friction(std::string input_file);

  // contructor from the direct input
  Friction(GreensTensor *greens_tensor, Polarizability *polarizability,
                  PowerSpectrum *powerspectrum, double relerr_omega);

  // calculation of the noncontact friction
  double calculate(Options_Friction opts);

  // the integrand of the noncontact friction
  static double friction_integrand(double omega, void *opts);
};
```

