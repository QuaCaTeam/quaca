#include "PowerSpectrum.h"
#include "../GreensTensor/GreensTensorFactory.h"
#include "../Polarizability/PolarizabilityFactory.h"
#include <string>

// ini parser
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

//Constructor using a ini-file
PowerSpectrum::PowerSpectrum(std::string input_file) {
  // Create a root
  pt::ptree root;

  // Load the ini file in this ptree
  pt::read_ini(input_file, root);

  //initialize Green's tensor by an input file
  this->greens_tensor = GreensTensorFactory::create(input_file);
  //initialize polarizability by an input file
  this->polarizability = PolarizabilityFactory::create(input_file);
};

//Constructor with initialization list
PowerSpectrum::PowerSpectrum(GreensTensor *greens_tensor,
                             Polarizability *polarizability)
    : greens_tensor(greens_tensor), polarizability(polarizability){};
