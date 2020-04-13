#ifndef GREENSTENSORPLATE_H
#define GREENSTENSORPLATE_H

#include <memory>

#include <armadillo>
#include <cassert>
#include <cmath>
#include <complex>

#include "../ReflectionCoefficients/ReflectionCoefficients.h"
#include "GreensTensor.h"

//! The class of the Green's tensor above a flat macroscopic surface
class GreensTensorPlate : public GreensTensor {
protected:
  double za;
  // reflection coefficients are needed to describe the surface's response
  std::shared_ptr<ReflectionCoefficients> reflection_coefficients;
  // kappa_cut defines the numerical cut-off of the kappa integration
  double delta_cut;
  vec::fixed<2> rel_err = {NAN, NAN};

public:
  // constructors
  GreensTensorPlate(double v, double beta, double za,
                    std::shared_ptr<ReflectionCoefficients> reflection_coefficients,
                    double delta_cut, const vec::fixed<2>& rel_err);
  explicit GreensTensorPlate(const std::string& input_file);

  // calculate the tensor in frequency and momentum space
  void calculate_tensor(cx_mat::fixed<3, 3> &GT, Options_GreensTensor opts) override;

  // integrate over a two-dimensional k space
  void integrate_k(cx_mat::fixed<3, 3> &GT, Options_GreensTensor opts) override;

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
  double omega_ch() override;

  // setter function
  void set_z_a(double za_new) { this->za = za_new; };
};

#endif // GREENSTENSORPLATE_H
