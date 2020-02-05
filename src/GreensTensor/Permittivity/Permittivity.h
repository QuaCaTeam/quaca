#ifndef PERMITTIVITY_H
#define PERMITTIVITY_H

#include <complex>

//! An abstract permittivity class
class Permittivity {
public:
  /*!
   * Return the permittivity given at frequency omega
   * @param omega Frequency
   */
  virtual std::complex<double> epsilon(double omega) = 0;
  virtual std::complex<double> epsilon_omega(double omega) = 0;
};

#endif // PERMITTIVITY_H
