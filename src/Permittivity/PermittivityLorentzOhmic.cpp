// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "PermittivityLorentzOhmic.h"

PermittivityLorentzOhmic::PermittivityLorentzOhmic(double eps_inf,
                                                   double alpha_zero,
                                                   double omega_0, double gamma)
    : eps_inf(eps_inf), alpha_zero(alpha_zero), omega_0(omega_0), gamma(gamma) {
}

PermittivityLorentzOhmic::PermittivityLorentzOhmic(
    const std::string &input_file) {
  // Load the json file in ptree
  pt::ptree root;
  pt::read_json(input_file, root);

  // check if type is right
  std::string type = root.get<std::string>("Permittivity.type");
  assert(type == "lorentz ohmic");

  // read parameters
  this->eps_inf = root.get<double>("Permittivity.eps_inf");
  this->alpha_zero = root.get<double>("Permittivity.alpha_zero");
  this->omega_0 = root.get<double>("Permittivity.omega_0");
  this->gamma = root.get<double>("Permittivity.gamma");
}

std::complex<double> PermittivityLorentzOhmic::calculate(double omega) const {
  // dummies for result and complex unit
  std::complex<double> result;
  std::complex<double> I(0.0, 1.0);

  // calculate the result
  result =
      eps_inf - alpha_zero * omega_0 * omega_0 /
                    (omega_0 * omega_0 - omega * omega - I * gamma * omega);

  return result;
}
