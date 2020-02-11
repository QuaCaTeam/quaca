// ini parser
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "OhmicMemoryKernel.h"

OhmicMemoryKernel::OhmicMemoryKernel(double gamma) : gamma(gamma){};

OhmicMemoryKernel::OhmicMemoryKernel(std::string input_file) {
  // Create a root
  pt::ptree root;

  // Load the ini file in this ptree
  pt::read_ini(input_file, root);

  // check if type is right
  std::string type = root.get<std::string>("MemoryKernel.type");
  assert(type == "ohmic");

  // read damping coefficient
  this->gamma = root.get<double>("MemoryKernel.gamma");
};

// return mu(omega) for defined memory kernel
std::complex<double> OhmicMemoryKernel::mu(double omega) {
  const std::complex<double> gammac(this->gamma, 0E0);
  return gammac;
};
