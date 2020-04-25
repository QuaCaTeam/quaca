#include <iostream>

// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "ReflectionCoefficientsFactory.h"
#include "ReflectionCoefficientsLocBulk.h"
#include "ReflectionCoefficientsLocSlab.h"

// reflection coefficients factory
ReflectionCoefficients *
ReflectionCoefficientsFactory::create(const std::string &input_file) {
  // return pointer
  ReflectionCoefficients *refcoef;

  // Create a root
  pt::ptree root;

  // Load the json file in this ptree
  pt::read_json(input_file, root);

  // read the type of reflection coefficient
  std::string type = root.get<std::string>("ReflectionCoefficients.type");

  // set the right pointer, show error if type is unknown
  if (type == "local bulk") {
    refcoef = new ReflectionCoefficientsLocBulk(input_file);
  } else if (type == "local slab") {
    refcoef = new ReflectionCoefficientsLocSlab(input_file);
  } else {
    std::cerr << "Error: Unknown Permittivity type (" << type << ")!"
              << std::endl;
    exit(0);
  }

  // return reflection coefficient pointer
  return refcoef;
}
