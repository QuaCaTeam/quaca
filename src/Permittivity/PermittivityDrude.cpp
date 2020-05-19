// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "PermittivityDrude.h"

PermittivityDrude::PermittivityDrude(double omega_p, double gamma)
    : omega_p(omega_p), gamma(gamma) {}

// constructor for drude model from .json file
PermittivityDrude::PermittivityDrude(const std::string &input_file) {

  // Create a root
  pt::ptree root;

  // Load the json file in this ptree
  pt::read_json(input_file, root);

  // check if type is right
  std::string type = root.get<std::string>("Permittivity.type");
  assert(type == "drude");

  // read parameters
  this->gamma = root.get<double>("Permittivity.gamma");
  this->omega_p = root.get<double>("Permittivity.omega_p");
}

// calculate the permittivity
std::complex<double> PermittivityDrude::calculate(double omega) const {
  // dummies for result and complex unit
  std::complex<double> result;
  std::complex<double> I(0.0, 1.0);

  // calculate the result
  result = 1.0 - omega_p * omega_p / (omega * (omega + I * gamma));

  return result;
}

// calculate the permittivity scaled by omega
std::complex<double>
PermittivityDrude::calculate_times_omega(double omega) const {
  // dummies for result and complex unit
  std::complex<double> result;
  std::complex<double> I(0.0, 1.0);

  // calculate the result
  result = omega - omega_p * omega_p / (omega + I * gamma);

  return result;
}

void PermittivityDrude::print_info(std::ostream &stream) const {
  stream << "# PermittivityDrude\n#\n"
         << "# omega_p = " << omega_p << "\n"
         << "# gamma = " << gamma << "\n";
}
