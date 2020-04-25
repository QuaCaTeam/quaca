// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "OhmicMemoryKernel.h"

OhmicMemoryKernel::OhmicMemoryKernel(double gamma) : gamma(gamma) {}

OhmicMemoryKernel::OhmicMemoryKernel(const std::string &input_file,
                                     const std::string &section) {
  // Create a root
  pt::ptree root;

  // Load the json file in this ptree
  pt::read_json(input_file, root);

  // check if type is right
  std::string type = root.get<std::string>(section + ".type");
  assert(type == "ohmic");

  // read damping coefficient
  this->gamma = root.get<double>(section + ".gamma");
}

OhmicMemoryKernel::OhmicMemoryKernel(const std::string &input_file) {
  // Create a root
  pt::ptree root;

  // Load the json file in this ptree
  pt::read_json(input_file, root);

  // check if type is right
  std::string type = root.get<std::string>("MemoryKernel.type");
  assert(type == "ohmic");

  // read damping coefficient
  this->gamma = root.get<double>("MemoryKernel.gamma");
}

// return mu(omega) for defined memory kernel
std::complex<double> OhmicMemoryKernel::mu(double omega) {
  const std::complex<double> gammac(this->gamma, 0E0);
  return gammac;
}

void OhmicMemoryKernel::print_info(std::ofstream &file) {
  file << "# OhmicMemoryKernel\n"
       << "# gamma = " << gamma << "\n";
}