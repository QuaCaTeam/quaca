#include <iostream>

#include <armadillo>
// ini parser
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "ReflectionCoefficientsFactory.h"
#include "ReflectionCoefficientsLocBulk.h"

// reflection coefficients factory
ReflectionCoefficients *
ReflectionCoefficientsFactory::create(std::string input_file) {
  // set return pointer to NULL
  ReflectionCoefficients *refcoef = NULL;

  // Create a root
  pt::ptree root;

  // Load the ini file in this ptree
  pt::read_ini(input_file, root);

  // read the type of reflection coefficient
  std::string type = root.get<std::string>("Reflection.type");

  // set the right pointer, show error if type is unknown
  if (type == "local bulk") {
    refcoef = new ReflectionCoefficientsLocBulk(input_file);
  } else {
    std::cerr << "Error: Unknown Permittivity type (" << type << ")!"
              << std::endl;
    exit(0);
  };

  // return reflection coefficient pointer
  return refcoef;
};
