#ifndef PERMITTIVITYLORENTZNOBATH_H
#define PERMITTIVITYLORENTZNOBATH_H

#include "Permittivity.h"
#include <complex>

//! A LorentzNoBath model permittivity
class PermittivityLorentzNoBath : public Permittivity {
private:
  double eps_inf;
  double alpha_zero;
  double omega_0;

public:
  // constructors
  PermittivityLorentzNoBath(double eps_inf, double alpha_zero, double omega_0);
  PermittivityLorentzNoBath(std::string input_file);

  // calculate the permittivity
  std::complex<double> epsilon(double omega);

  // Returns the numerical value of the permittivity scaled by omega.
  std::complex<double> epsilon_omega(double omega);

  // getter methods
  double get_eps_inf() { return this->eps_inf; };
  double get_alpha_zero() { return this->alpha_zero; };
  double get_omega_0() { return this->omega_0; };
};

#endif // PERMITTIVITYLORENTZNOBATH_H
