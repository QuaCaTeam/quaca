// ini parser
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "GreensTensor.h"

GreensTensor::GreensTensor(double v, double beta) : v(v), beta(beta) {
  assert(v >= 0 && v < 1);
  assert(beta > 0);
};

GreensTensor::GreensTensor(std::string input_file) {
  // Create a root
  pt::ptree root;

  // Load the ini file in this ptree
  pt::read_ini(input_file, root);

  // read parameters
  this->v = root.get<double>("GreensTensor.v");
  assert(v >= 0 && v < 1);
  this->beta = root.get<double>("GreensTensor.beta");
  assert(beta > 0);
};
