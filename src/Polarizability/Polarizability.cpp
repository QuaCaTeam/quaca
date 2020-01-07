// ini parser
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
namespace pt = boost::property_tree;

#include "../GreensTensor/GreensTensorFactory.h"
#include "Polarizability.h"

Polarizability::Polarizability(std::string input_file)
{
  // Create a root
  pt::ptree root;

  // Load the ini file in this ptree
  pt::read_ini(input_file, root);

  // read parameters
  this->omega_a = root.get<double>("Polarizability.omega_a");
  this->alpha_zero = root.get<double>("Polarizability.alpha_zero");

  // read greens tensor
  this->greens_tensor = GreensTensorFactory::create(input_file);
};
