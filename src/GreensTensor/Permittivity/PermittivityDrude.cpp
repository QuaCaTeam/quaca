// ini parser
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "PermittivityDrude.h"

PermittivityDrude::PermittivityDrude(double omega_p, double gamma)
    : omega_p(omega_p), gamma(gamma){};

// constructor for drude model from .ini file
PermittivityDrude::PermittivityDrude(std::string input_file) {
  // Create a root
  pt::ptree root;

  // Load the ini file in this ptree
  pt::read_ini(input_file, root);

  // read parameters
  this->gamma = root.get<double>("Permittivity.gamma");
  this->omega_p = root.get<double>("Permittivity.omega_p");
};

// calculate the permittivity
std::complex<double> PermittivityDrude::epsilon(double omega) {
  // dummies for result and complex unit
  std::complex<double> result;
  std::complex<double> I(0.0, 1.0);

  // calculate the result
  result = 1.0 - omega_p * omega_p / (omega * (omega + I * gamma));

  return result;
};
