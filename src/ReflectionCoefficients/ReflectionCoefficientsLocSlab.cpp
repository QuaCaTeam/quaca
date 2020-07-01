#include <armadillo>

// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <utility>
namespace pt = boost::property_tree;

#include "ReflectionCoefficientsLocSlab.h"

// direct constructor
ReflectionCoefficientsLocSlab::ReflectionCoefficientsLocSlab(
    std::shared_ptr<Permittivity> permittivity, double thickness)
    : permittivity(std::move(permittivity)), thickness(thickness) {}

// constructor from .json file
ReflectionCoefficientsLocSlab::ReflectionCoefficientsLocSlab(
    const std::string &input_file) {
  // set permittivity
  // set parameters
  this->permittivity = PermittivityFactory::create(input_file);
  // Create a root
  pt::ptree root;

  // Load the json file in this ptree
  pt::read_json(input_file, root);

  // read parameters
  this->thickness = root.get<double>("ReflectionCoefficients.thickness");
}

// calculate the p-polarized reflection coefficient
void ReflectionCoefficientsLocSlab::calculate(double omega,
                                              std::complex<double> kappa,
                                              std::complex<double> &r_p,
                                              std::complex<double> &r_s) const {
  // absolute value of omega. r_p is always calculated for positive omega and if
  // needed complex conjugated after the calculation
  double omega_abs = std::abs(omega);
  std::complex<double> eps = this->permittivity->calculate(omega_abs);
  std::complex<double> eps_omega = this->permittivity->calculate_times_omega(omega_abs);
  std::complex<double> I(0., 1.);

  // kapppa as well as kappa_epsilon are defined to have either a purely
  // positive real part or purely negatively imaginary part
  std::complex<double> kappa_epsilon =
      sqrt(kappa * kappa - (eps - 1.) * omega_abs * omega_abs);
  kappa_epsilon = std::complex<double>(std::abs(kappa_epsilon.real()),
                                       -std::abs(kappa_epsilon.imag()));
  // Defining the reflection coefficients in transverse magnetice polarization
  // (p) and in transverse electric polarization (s)
  std::complex<double> r_p_bulk =
      (kappa * eps_omega - kappa_epsilon * omega_abs) / (kappa * eps_omega + kappa_epsilon * omega_abs);
  std::complex<double> r_s_bulk =
      (kappa - kappa_epsilon) / (kappa + kappa_epsilon);

  r_p = r_p_bulk * (1. - exp(-2. * kappa_epsilon * this->thickness)) /
        (1. - pow(r_p_bulk * exp(-kappa_epsilon * this->thickness), 2));
  r_s = r_s_bulk * (1. - exp(-2. * kappa_epsilon * this->thickness)) /
        (1. - pow(r_s_bulk * exp(-kappa_epsilon * this->thickness), 2));

  // Imposing crossing relation
  if (omega < 0.) {
    r_p = conj(r_p);
    r_s = conj(r_s);
  }
}
void ReflectionCoefficientsLocSlab::print_info(std::ostream &stream) const {
  stream << "# ReflectionCoefficientsLocSlab\n#\n"
         << "# thickness = " << thickness << "\n";
  permittivity->print_info(stream);
}
