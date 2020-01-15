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
    double za;
public:
  GreensTensorPlate(double v, double za, double beta, Permittivity *permittivity, double delta_cut):
  GreensTensor(v, beta) {this->za = za; this->delta_cut = delta_cut; this->permittivity = permittivity;};


  GreensTensorPlate(std::string input_file): GreensTensor(input_file) {this->permittivity = PermittivityFactory::create(input_file);};

  void calculate_tensor(cx_mat::fixed<3,3>& GT, vec::fixed<2> kvec, double omega);
  void integrate_k_2d(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts);
  void integrate_k_1d(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts);

  std::complex<double> get_epsilon(double omega){return this->permittivity->epsilon(omega);};
  double get_za(){return this->za;}
  double get_delta_cut(){return this->delta_cut;}
  static double integrand_k_1d(double kx, void* opts);
  static double integrand_k_2d(double ky, void* opts);

};


#endif // GREENSTENSORPLATE_H
