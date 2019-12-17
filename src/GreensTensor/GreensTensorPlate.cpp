// ini parser
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
namespace pt = boost::property_tree;

#include "GreensTensorPlate.h"

GreensTensorPlate::GreensTensorPlate(std::string input_file)
{
  // Create a root
  pt::ptree root;

  // Load the ini file in this ptree
  pt::read_ini(input_file, root);

  // read parameters
  this->v = root.get<double>("GreensTensor.v");
  this->z_a = root.get<double>("GreensTensor.z_a");
  this->beta = root.get<double>("GreensTensor.beta");
}
