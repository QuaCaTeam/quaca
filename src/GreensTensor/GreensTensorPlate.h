#ifndef GREENSTENSORPLATE_H
#define GREENSTENSORPLATE_H

#include <armadillo>
#include <assert.h>
#include <cmath>
#include <complex>

#include "../ReflectionCoefficients/ReflectionCoefficients.h"
#include "GreensTensor.h"

//! The class of the Green's tensor above a flat macroscopic surface
class GreensTensorPlate : public GreensTensor {
protected:
  // reflection coefficients are needed to describe the surface's response
  ReflectionCoefficients *reflection_coefficients;
  // kappa_cut defines the numerical cut-off of the kappa integration
  double delta_cut;
  vec::fixed<2> rel_err = {NAN, NAN};
  double za;

public:
  // constructors
  GreensTensorPlate(double v, double za, double beta,
                    ReflectionCoefficients *reflection_coefficients,
                    double delta_cut, vec::fixed<2> rel_err);
  GreensTensorPlate(std::string input_file);

  // calculate the tensor in frequency and momentum space
  void calculate_tensor(cx_mat::fixed<3, 3> &GT, Options_GreensTensor opts);

  // integrate over a two-dimensional k space
  void integrate_k(cx_mat::fixed<3, 3> &GT, Options_GreensTensor opts);

  // integrands
  static double integrand_1d_k(double kx, void *opts);
  static double integrand_2d_k(double ky, void *opts);

  // getter functions
  std::complex<double> get_r_p(double omega, double k);
  std::complex<double> get_r_s(double omega, double k);

  double get_za() { return this->za; };
  double get_delta_cut() { return this->delta_cut; };
  double get_rel_err_0() { return this->rel_err(0); };
  double get_rel_err_1() { return this->rel_err(1); };
  double omega_ch();

  // setter function
  void set_z_a(double za) { this->za = za; };
};

#endif // GREENSTENSORPLATE_H
