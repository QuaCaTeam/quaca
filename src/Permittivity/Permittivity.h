#ifndef PERMITTIVITY_H
#define PERMITTIVITY_H

#include <complex>

//! An abstract permittivity class
class Permittivity {
public:
  // calculate the permittivity
  virtual std::complex<double> calculate(double omega) const = 0;

  // calculate the permittivity times omega
  virtual std::complex<double> calculate_times_omega(double omega) const = 0;

  // print info
  virtual void print_info(std::ostream &stream) const =0;
};

#endif // PERMITTIVITY_H
