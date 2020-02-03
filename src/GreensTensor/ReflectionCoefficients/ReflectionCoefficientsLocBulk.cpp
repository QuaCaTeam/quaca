#include "ReflectionCoefficientsLocBulk.h"

// direct constructor
ReflectionCoefficientsLocBulk::ReflectionCoefficientsLocBulk(
    Permittivity *permittivity) {
  // set permittivity
  // set parameters
  this->permittivity = permittivity;
};

// calculate the p-polarized reflection coefficient
void ReflectionCoefficientsLocBulk::refvec(cx_vec ref(2), double omega,
                                           std::complex<double> kappa) {
  // dummies for result and complex unit
  std::complex<double> result;
  // absolute value of omega. r_p is always calculated for positive omega and if
  // needed complex conjugated after the calculation
  double omega_abs;
  double eps = this->permittivity->epsilon(omega_abs);
  std::complex < double kappa_epsilon;

  // kapppa as well as kappa_epsilon are defined to have either a purely
  // positive real part or purely negatively imaginary part
  kappa_epsilon = sqrt(k_quad - eps * omega_abs * omega_abs);
  kappa_epsilon = std::complex<double>(std::abs(kappa_epsilon.real()),
                                       -std::abs(kappa_epsilon.imag()));

  // Defining the reflection coefficients in transverse magnetice polarization
  // (p) and in transverse electric polarization (s)
  ref(0) = (kappa * eps - kappa_epsilon) / (kappa * eps + kappa_epsilon);
  ref(1) = (kappa - kappa_epsilon) / (kappa + kappa_epsilon);

  // Imposing crossing relation
  if (omega < 0.) {
    ref(0) = conj(ref(0));
    ref(1) = conj(ref(1));
  };
};
