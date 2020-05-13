#ifndef GREENSTENSORPLATE_H
#define GREENSTENSORPLATE_H

#include <armadillo>
#include <assert.h>
#include <cmath>
#include <complex>
#include <memory>

#include "../ReflectionCoefficients/ReflectionCoefficients.h"
#include "GreensTensor.h"

//! The class of the Green's tensor above a flat macroscopic surface
class GreensTensorPlate : public GreensTensor {
protected:
  //distance between microscopic object and macroscopic surface
  double za;

  // kappa_cut defines the numerical cut-off of the kappa integration
  double delta_cut;
  vec::fixed<2> rel_err = {NAN, NAN};

  // reflection coefficients are needed to describe the surface's response
  std::shared_ptr<ReflectionCoefficients> reflection_coefficients;

public:
  // constructors
  GreensTensorPlate(
      double v, double beta, double za,
      std::shared_ptr<ReflectionCoefficients> reflection_coefficients,
      double delta_cut, const vec::fixed<2> &rel_err);
  explicit GreensTensorPlate(const std::string &input_file);

  // calculate the tensor in frequency and momentum space
  void calculate_tensor(double omega, vec::fixed<2> k,
                        cx_mat::fixed<3, 3> &GT) const;

  // integrate over a two-dimensional k space
  void integrate_k(double omega, cx_mat::fixed<3, 3> &GT,
                   Tensor_Options fancy_complex,
                   Weight_Options weight_function) const;

  // integrands
  double integrand_1d_k(double phi, double omega, const vec::fixed<2> &indices,
                        Tensor_Options fancy_complex,
                        Weight_Options weight_function) const;

  double integrand_2d_k(double kappa_double, double omega, double phi,
                        const vec::fixed<2> &indices,
                        Tensor_Options fancy_complex,
                        Weight_Options weight_function) const;

  // getter functions
  std::complex<double> get_r_p(double omega, double k) const;
  std::complex<double> get_r_s(double omega, double k) const;

  double get_za() const { return this->za; };
  double get_delta_cut() const { return this->delta_cut; };
  double get_rel_err_0() const { return this->rel_err(0); };
  double get_rel_err_1() const { return this->rel_err(1); };
  double omega_ch() const;

  // setter function
  void set_za(double za_new) { this->za = za_new; };
};

#endif // GREENSTENSORPLATE_H
