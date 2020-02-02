#ifndef POLARIZABILITY_H
#define POLARIZABILITY_H

#include <armadillo>
#include <cmath>
#include <complex>

#include "../GreensTensor/GreensTensor.h"
#include "MemoryKernel/MemoryKernel.h"

using namespace arma;

//! a struct with the integration options
struct Options_Polarizabiliy;

//! An abstract polarizability class
/*!
 * This is an abstract polarizability class.
 * All polarizabilities should return a 3 by 3 matrix, given a real frequency as
 * input. The only two child class of this will be 1) a polarizability where the
 * particle has no interval bath 2) a polarizability where the particle has an
 * internal bath
 */
class Polarizability {
protected:
  /*! resonance frequency */
  double omega_a;

  /*! static polarizability */
  double alpha_zero;

  /*! green's tensor */
  GreensTensor *greens_tensor;

public:
  /*!
   * Returns the polarizability tensor as a 3x3 matrix when given a real
   * frequency as input.
   * @param omega Frequency
   */
  virtual void calculate(cx_mat::fixed<3, 3> &alpha, double omega) = 0;

  /*!
   * Getter method for the resonance frequency.
   */
  double get_omega_a();

  /*!
   * Getter method for the static polarizability.
   */
  double get_alpha_zero();
};

// A struct for integration options
struct Options_Polarizability {
  // Different options for the integrand
  bool fancy_R = false;
  bool fancy_I = false;
  bool fancy_I_kv = false;
  bool fancy_I_temp = false;
  bool fancy_I_kv_temp = false;
  // Indices of the 3x3 polarizability tensor
  arma::vec::fixed<2> indices = {-1, -1};
  // Value of omega for the integration of the k-Variables
  double omega = NAN;
  // Pointer to the Polarizability to be able to access the attributes of the
  // class eventhough the integrand is static
  Polarizability *class_pt;
};

#endif // POLARIZABILITY_H
