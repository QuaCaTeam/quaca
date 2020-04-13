#ifndef PERMITTIVITYDRUDE_H
#define PERMITTIVITYDRUDE_H

#include "Permittivity.h"
#include <complex>

//! A Drude model permittivity
class PermittivityDrude : public Permittivity {
private:
  double omega_p, gamma; // plasma frequency and damping coefficient

public:
  // constructors
  PermittivityDrude(double omega_p, double gamma);
  PermittivityDrude(const std::string& input_file);

  // calculate the permittivity
  std::complex<double> epsilon(double omega);

  // Returns the numerical value of the permittivity scaled by omega.
  std::complex<double> epsilon_omega(double omega);

  // getter methods
  double get_gamma() { return this->gamma; };
  double get_omega_p() { return this->omega_p; };
};

#endif // PERMITTIVITYDRUDE_H
