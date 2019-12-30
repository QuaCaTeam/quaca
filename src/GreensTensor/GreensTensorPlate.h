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


public:

  GreensTensorPlate(double v, double za, double beta, std::string input_file);
  void calculate_tensor(cx_mat::fixed<3,3>& GT, vec::fixed<2> kvec, double omega);
  void integrate_k_2d(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts);
  void integrate_k_1d(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts);

  std::complex<double> get_epsilon(double omega);

  static double integrand_k_1d(double k, void* opts);

};


#endif // GREENSTENSORPLATE_H
