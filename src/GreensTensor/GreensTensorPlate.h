#ifndef GREENSTENSORPLATE_H
#define GREENSTENSORPLATE_H

#include <complex>
#include <cmath>
#include "GreensTensor.h"
#include "Permittivity/PermittivityFactory.h"

class GreensTensorPlate : public GreensTensor
{
private:
    // permittivity is needed to describe the surface's response
    Permittivity *permittivity;
    // kappa_cut defines the numerical cut-off of the kappa integration
    double delta_cut;
public:
    PolarizabilityBath(double omega_a, double alpha_zero, MemoryKernel *mu, GreensTensor *greens_tensor): Polarizability(omega_a, alpha_zero, greens_tensor) {this->mu = mu;};
    PolarizabilityBath(std::string input_file): Polarizability(input_file) {this->mu = MemoryKernelFactory::create(input_file);};

    std::complex<double> get_mu(double omega){return this->mu->mu(omega);};

  GreensTensorPlate(double v, double za, double beta, std::string input_file);
  void calculate_tensor(cx_mat::fixed<3,3>& GT, vec::fixed<2> kvec, double omega);
  void integrate_k_2d(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts);
  void integrate_k_1d(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts);

  std::complex<double> get_epsilon(double omega);

  static double integrand_k_1d(double kx, void* opts);
  static double integrand_k_2d(double ky, void* opts);

};


#endif // GREENSTENSORPLATE_H
