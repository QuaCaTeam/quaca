#ifndef PERMITTIVITY_H
#define PERMITTIVITY_H

#include <complex>

//! An abstract permittivity class
class Permittivity {
public:
  // calculate the permittivity
  virtual std::complex<double> epsilon(double omega) = 0;

  // calculate the permittivity times omega
  virtual std::complex<double> epsilon_omega(double omega) = 0;
};

#endif // PERMITTIVITY_H
