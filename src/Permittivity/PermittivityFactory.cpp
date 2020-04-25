#include <iostream>

// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "PermittivityDrude.h"
#include "PermittivityFactory.h"
#include "PermittivityLorentz.h"

// permittivity factory
Permittivity *PermittivityFactory::create(const std::string& input_file) {
  // return pointer
  Permittivity *permittivity;

  // Create a root
  pt::ptree root;

  // Load the json file in this ptree
  pt::read_json(input_file, root);

  // read the type of permittivity
  std::string type = root.get<std::string>("Permittivity.type");

  // set the right pointer, show error if type is unknown
  if (type == "drude") {
    permittivity = new PermittivityDrude(input_file);
  } else if (type == "lorentz") {
    permittivity = new PermittivityLorentz(input_file);
  } else {
    std::cerr << "Error: Unknown Permittivity type (" << type << ")!"
              << std::endl;
    exit(0);
  }

  // return permittivity pointer
  return permittivity;
}
