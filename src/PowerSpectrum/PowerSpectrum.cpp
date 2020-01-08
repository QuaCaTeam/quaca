#include "PowerSpectrum.h"
#include "../GreensTensor/GreensTensorFactory.h"
#include "../Polarizability/PolarizabilityFactory.h"
#include <string>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

namespace pt = boost::property_tree;

PowerSpectrum::PowerSpectrum(std::string input_file)
{
  // Create a root
  pt::ptree root;

  // Load the ini file in this ptree
  pt::read_ini(input_file, root);

  // read greens tensor
  this->greens_tensor = GreensTensorFactory::create(input_file);
  this->polarizability = PolarizabilityFactory::create(input_file);

};

PowerSpectrum::PowerSpectrum(GreensTensor* greens_tensor, Polarizability* polarizability): greens_tensor(greens_tensor), polarizability(polarizability)
{
};

