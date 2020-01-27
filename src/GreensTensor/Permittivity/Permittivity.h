#ifndef PERMITTIVITY_H
#define PERMITTIVITY_H

#include <complex>

//! An abstract permittivity class
class Permittivity {
protected:
  std::string type; // type of permittivity

public:
  // constructor
  Permittivity();
  Permittivity(std::string input_file);

  // Return the permittivity given at frequency omega
  virtual std::complex<double> epsilon(double omega) = 0;
};

#endif // PERMITTIVITY_H
