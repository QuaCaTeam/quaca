// ini parser
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
namespace pt = boost::property_tree;

#include "GreensTensor.h"

GreensTensor::GreensTensor(std::string input_file)
{
  // Create a root
  pt::ptree root;

  // Load the ini file in this ptree
  pt::read_ini(input_file, root);

  // read parameters
  this->v = root.get<double>("GreensTensor.v");
  this->za = root.get<double>("GreensTensor.za");
  this->beta = root.get<double>("GreensTensor.beta");
};
