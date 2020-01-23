#ifndef POLARIZABILITY_H
#define POLARIZABILITY_H

#include <armadillo>
#include <cmath>
#include <complex>

#include "../GreensTensor/GreensTensor.h"
#include "MemoryKernel/MemoryKernel.h"

using namespace arma;

//! a struct with the integration options
struct Options_Polarizability;

//! An abstract polarizability class
class Polarizability {
protected:
  double omega_a;              // resonance frequency
  double alpha_zero;           // vacuum polarizability
  GreensTensor *greens_tensor; // Green's tensor

public:
  // constructors
  Polarizability(double omega_a, double alpha_zero,
                 GreensTensor *greens_tensor);
  Polarizability(std::string input_file);

  // calculate the polarizability tensor
  virtual void calculate_tensor(cx_mat::fixed<3, 3> &alpha,
                                Options_Polarizability opts) = 0;

  // integration over omega
  double integrate_omega(Options_Polarizability opts, double omega_min,
                         double omega_max, double relerr, double abserr);
  static double integrand_omega(double omega, void *opts);

  // getter functions
  double get_omega_a() { return this->omega_a; };
  double get_alpha_zero() { return this->alpha_zero; };
};

// A struct for integration options
struct Options_Polarizability {
  // Different options for the integrand
  bool fancy_R = false;
  bool fancy_I = false;

  double omega = NAN;

  // Indices of the 3x3 polarizability tensor
  vec::fixed<2> indices = {NAN, NAN};

  // Pointer to the Polarizability to be able to access the attributes of the
  // class eventhough the integrand is static
  Polarizability *class_pt;
};

#endif // POLARIZABILITY_H
