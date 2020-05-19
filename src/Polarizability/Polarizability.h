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
private:
  double omega_a;                              // resonance frequency
  double alpha_zero;                           // vacuum polarizability
  std::shared_ptr<MemoryKernel> mu;            // internal bath
  std::shared_ptr<GreensTensor> greens_tensor; // Green's tensor

public:
  // Direct constructor without internal bath
  Polarizability(double omega_a, double alpha_zero,
                 std::shared_ptr<GreensTensor> greens_tensor);

  // Direct constructor with internal bath mu
  Polarizability(double omega_a, double alpha_zero,
                 std::shared_ptr<MemoryKernel> mu,
                 std::shared_ptr<GreensTensor> greens_tensor);

  // Constructor from a given json file
  explicit Polarizability(const std::string &input_file);

  // calculate the polarizability tensor
  void calculate_tensor(double omega, cx_mat::fixed<3, 3> &alpha,
                        Tensor_Options fancy_complex) const;

  // integration over omega
  double integrate_omega(const vec::fixed<2> &indices,
                         Tensor_Options fancy_complex, double omega_min,
                         double omega_max, double relerr, double abserr) const;

  // integrand for the omega integration
  double integrand_omega(double omega, const vec::fixed<2> &indices,
                         Tensor_Options fancy_complex) const;

  // getter functions
  double get_omega_a() const { return omega_a; };
  double get_alpha_zero() const { return alpha_zero; };
  std::shared_ptr<GreensTensor> get_greens_tensor() const {
    return greens_tensor;
  };

  // getter function for memory kernel
  std::complex<double> get_mu(double omega) const {
    if (mu != nullptr) {
      return mu->calculate(omega);
    } else {
      return std::complex<double>(0.0, 0.0);
    }
  };

  // print info
  void print_info(std::ostream &stream) const;
};

#endif // POLARIZABILITY_H
