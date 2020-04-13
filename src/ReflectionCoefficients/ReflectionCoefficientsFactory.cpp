#include <iostream>

// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "ReflectionCoefficientsFactory.h"
#include "ReflectionCoefficientsLocBulk.h"
#include "ReflectionCoefficientsLocSlab.h"

// reflection coefficients factory
std::shared_ptr<ReflectionCoefficients>
ReflectionCoefficientsFactory::create(const std::string& input_file) {
  // Create a root
  // Load the json file in this ptree
  pt::ptree root;
  pt::read_json(input_file, root);

  // read the type of reflection coefficient
  std::string type = root.get<std::string>("ReflectionCoefficients.type");

  // set the right pointer, show error if type is unknown
  if (type == "local bulk") {
    return std::make_shared<ReflectionCoefficientsLocBulk>(input_file);
  } else if (type == "local slab") {
    return std::make_shared<ReflectionCoefficientsLocSlab>(input_file);
  } else {
    std::cerr << "Error: Unknown Permittivity type (" << type << ")!"
              << std::endl;
    exit(0);
  }
}
