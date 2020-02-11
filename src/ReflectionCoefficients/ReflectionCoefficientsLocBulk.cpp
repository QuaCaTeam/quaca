#include "ReflectionCoefficientsLocBulk.h"
#include <armadillo>
// direct constructor
ReflectionCoefficientsLocBulk::ReflectionCoefficientsLocBulk(
    Permittivity *permittivity) {
  // set permittivity
  // set parameters
  this->permittivity = permittivity;
};
// constructor from .ini file
ReflectionCoefficientsLocBulk::ReflectionCoefficientsLocBulk(std::string input_file) {
  // set permittivity
  // set parameters
    this->permittivity = PermittivityFactory::create(input_file);
};

// calculate the p-polarized reflection coefficient
void ReflectionCoefficientsLocBulk::ref(std::complex<double> &r_p, std::complex<double> &r_s, double omega, std::complex<double> kappa){
  // absolute value of omega. r_p is always calculated for positive omega and if
  // needed complex conjugated after the calculation
  double omega_abs = std::abs(omega);
  std::complex<double> eps = this->permittivity->epsilon(omega_abs);
  std::complex<double> kappa_epsilon;

  // kapppa as well as kappa_epsilon are defined to have either a purely
  // positive real part or purely negatively imaginary part
  kappa_epsilon = sqrt(kappa*kappa - (eps - 1.) * omega_abs * omega_abs);
  kappa_epsilon = std::complex<double>(std::abs(kappa_epsilon.real()),
                                       -std::abs(kappa_epsilon.imag()));
  // Defining the reflection coefficients in transverse magnetice polarization
  // (p) and in transverse electric polarization (s)
  r_p = (kappa * eps - kappa_epsilon) / (kappa * eps + kappa_epsilon);
  r_s = (kappa - kappa_epsilon) / (kappa + kappa_epsilon);
  // Imposing crossing relation
  if (omega < 0.) {
    r_p = conj(r_p);
    r_s = conj(r_s);
  };
};
