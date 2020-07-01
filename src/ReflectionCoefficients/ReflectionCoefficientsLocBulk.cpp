#include "ReflectionCoefficientsLocBulk.h"
#include <armadillo>
#include <utility>
// direct constructor
ReflectionCoefficientsLocBulk::ReflectionCoefficientsLocBulk(
    std::shared_ptr<Permittivity> permittivity)
    : permittivity(std::move(permittivity)) {}

// constructor from .json file
ReflectionCoefficientsLocBulk::ReflectionCoefficientsLocBulk(
    const std::string &input_file) {
  // set permittivity
  this->permittivity = PermittivityFactory::create(input_file);
}

// calculate the p-polarized reflection coefficient
void ReflectionCoefficientsLocBulk::calculate(double omega,
                                              std::complex<double> kappa,
                                              std::complex<double> &r_p,
                                              std::complex<double> &r_s) const {
  // absolute value of omega. r_p is always calculated for positive omega and if
  // needed complex conjugated after the calculation
  double omega_abs = std::abs(omega);
  std::complex<double> eps = this->permittivity->calculate(omega_abs);
  std::complex<double> eps_omega = this->permittivity->calculate_times_omega(omega_abs);

  // kapppa as well as kappa_epsilon are defined to have either a purely
  // positive real part or purely negatively imaginary part
  std::complex<double> kappa_epsilon =
      sqrt(kappa * kappa - (eps - 1.) * omega_abs * omega_abs);
  kappa_epsilon = std::complex<double>(std::abs(kappa_epsilon.real()),
                                       -std::abs(kappa_epsilon.imag()));
  // Defining the reflection coefficients in transverse magnetice polarization
  // (p) and in transverse electric polarization (s)
  r_p = (kappa * eps_omega - kappa_epsilon*omega_abs) / (kappa * eps_omega + kappa_epsilon*omega_abs);
  r_s = (kappa - kappa_epsilon) / (kappa + kappa_epsilon);
  // Imposing crossing relation
  if (omega < 0.) {
    r_p = conj(r_p);
    r_s = conj(r_s);
  }
}

void ReflectionCoefficientsLocBulk::print_info(std::ostream &stream) const {
  stream << "# ReflectionCoefficientsLocBulk\n#\n";
  permittivity->print_info(stream);
}
