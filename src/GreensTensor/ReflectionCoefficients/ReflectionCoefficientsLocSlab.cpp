#include "ReflectionCoefficientsLocSlab.h"
#include <armadillo>
// ini parser
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;
// direct constructor
ReflectionCoefficientsLocSlab::ReflectionCoefficientsLocSlab(
    Permittivity *permittivity, double thickness) {
  // set permittivity
  // set parameters
  this->permittivity = permittivity;
  this->thickness = thickness;
};
// constructor from .ini file
ReflectionCoefficientsLocSlab::ReflectionCoefficientsLocSlab(std::string input_file) {
  // set permittivity
  // set parameters
    this->permittivity = PermittivityFactory::create(input_file);
    // Create a root
    pt::ptree root;

    // Load the ini file in this ptree
    pt::read_ini(input_file, root);

    // read parameters
    this->thickness = root.get<double>("Reflection.thickness");
};

// calculate the p-polarized reflection coefficient
void ReflectionCoefficientsLocSlab::ref(std::complex<double> &r_p, std::complex<double> &r_s, double omega, std::complex<double> kappa){
  // absolute value of omega. r_p is always calculated for positive omega and if
  // needed complex conjugated after the calculation
  double omega_abs = std::abs(omega);
  std::complex<double> eps = this->permittivity->epsilon(omega_abs);
  std::complex<double> kappa_epsilon;
  std::complex<double> I(0.,1.);
  std::complex<double> r_p_bulk, r_s_bulk;

  // kapppa as well as kappa_epsilon are defined to have either a purely
  // positive real part or purely negatively imaginary part
  kappa_epsilon = sqrt(kappa*kappa - (eps - 1.) * omega_abs * omega_abs);
  kappa_epsilon = std::complex<double>(std::abs(kappa_epsilon.real()),
                                       -std::abs(kappa_epsilon.imag()));
  // Defining the reflection coefficients in transverse magnetice polarization
  // (p) and in transverse electric polarization (s)
  r_p_bulk = (kappa * eps - kappa_epsilon) / (kappa * eps + kappa_epsilon);
  r_s_bulk = (kappa - kappa_epsilon) / (kappa + kappa_epsilon);
  
  r_p = r_p_bulk * (1. - exp(-2.*kappa_epsilon*this->thickness))/(1. - pow(r_p_bulk * exp(-kappa_epsilon*this->thickness), 2));
  r_s = r_s_bulk * (1. - exp(-2.*kappa_epsilon*this->thickness))/(1. - pow(r_s_bulk * exp(-kappa_epsilon*this->thickness), 2));

  // Imposing crossing relation
  if (omega < 0.) {
    r_p = conj(r_p);
    r_s = conj(r_s);
  };
};
