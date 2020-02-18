// ini parser
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "PermittivityLorentzNoBath.h"

PermittivityLorentzNoBath::PermittivityLorentzNoBath(double eps_inf,
                                                     double alpha_zero,
                                                     double omega_0)
    : eps_inf(eps_inf), alpha_zero(alpha_zero), omega_0(omega_0){};

// constructor for drude model from .ini file
PermittivityLorentzNoBath::PermittivityLorentzNoBath(std::string input_file) {

  // Create a root
  pt::ptree root;

  // Load the ini file in this ptree
  pt::read_ini(input_file, root);

  // check if type is right
  std::string type = root.get<std::string>("Permittivity.type");
  assert(type == "lorentz nobath");

  // read parameters
  this->eps_inf = root.get<double>("Permittivity.eps_inf");
  this->alpha_zero = root.get<double>("Permittivity.alpha_zero");
  this->omega_0 = root.get<double>("Permittivity.omega_0");
};

// calculate the permittivity
std::complex<double> PermittivityLorentzNoBath::epsilon(double omega) {
  // dummies for result and complex unit
  std::complex<double> result;
  std::complex<double> I(0.0, 1.0);

  // calculate the result
  result = eps_inf -
           alpha_zero * omega_0 * omega_0 / (omega_0 * omega_0 - omega * omega);

  return result;
};

// calculate the permittivity scaled by omega
std::complex<double> PermittivityLorentzNoBath::epsilon_omega(double omega) {
  // dummies for result and complex unit
  std::complex<double> result;
  std::complex<double> I(0.0, 1.0);

  // calculate the result
  result = eps_inf * omega - alpha_zero * omega_0 * omega_0 * omega /
                                 (omega_0 * omega_0 - omega * omega);

  return result;
};
