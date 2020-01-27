// ini parser
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "OhmicMemoryKernel.h"

OhmicMemoryKernel::OhmicMemoryKernel(double gamma) : gamma(gamma) {
  this->type = "ohmic";
};

OhmicMemoryKernel::OhmicMemoryKernel(std::string input_file)
    : MemoryKernel(input_file) {
  // check if type is right
  assert(this->type == "ohmic");

  // Create a root
  pt::ptree root;

  // Load the ini file in this ptree
  pt::read_ini(input_file, root);

  // read damping coefficient
  this->gamma = root.get<double>("MemoryKernel.gamma");
};

// return mu(omega) for defined memory kernel
std::complex<double> OhmicMemoryKernel::mu(double omega) {
  const std::complex<double> gammac(this->gamma, 0E0);
  return gammac;
};
