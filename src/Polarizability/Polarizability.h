#ifndef POLARIZABILITY_H
#define POLARIZABILITY_H

#include <armadillo>
#include <cmath>
#include <complex>
#include <memory>

#include "../GreensTensor/GreensTensor.h"
#include "../MemoryKernel/MemoryKernel.h"

using namespace arma;

//! An abstract polarizability class
class Polarizability {
protected:
  double omega_a;                              // resonance frequency
  double alpha_zero;                           // vacuum polarizability
  std::shared_ptr<GreensTensor> greens_tensor; // Green's tensor

public:
  // constructors
  Polarizability(double omega_a, double alpha_zero,
                 std::shared_ptr<GreensTensor> greens_tensor);
  explicit Polarizability(const std::string &input_file);

  // calculate the polarizability tensor
  virtual void calculate_tensor(double omega, cx_mat::fixed<3, 3> &alpha,
                                Tensor_Options fancy_complex) const = 0;

  // integration over omega
  double integrate_omega(const vec::fixed<2> &indices, Tensor_Options fancy_complex,
                         double omega_min, double omega_max, double relerr,
                         double abserr) const;
  double integrand_omega(double omega, const vec::fixed<2> &indices,
                         Tensor_Options fancy_complex) const;

  // getter functions
  double get_omega_a() const { return omega_a; };
  double get_alpha_zero() const { return alpha_zero; };
  std::shared_ptr<GreensTensor> &get_greens_tensor() { return greens_tensor; };
};

#endif // POLARIZABILITY_H
