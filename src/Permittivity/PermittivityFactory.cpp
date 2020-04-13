#include <iostream>

// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "PermittivityDrude.h"
#include "PermittivityFactory.h"
#include "PermittivityLorentz.h"

// permittivity factory
std::shared_ptr<Permittivity>
PermittivityFactory::create(const std::string &input_file) {

  // Create a root
  // Load the json file in this ptree
  pt::ptree root;
  pt::read_json(input_file, root);

  // read the type of permittivity
  std::string type = root.get<std::string>("Permittivity.type");

  // set the right pointer, show error if type is unknown
  if (type == "drude") {
    return std::make_shared<PermittivityDrude>(input_file);
  } else if (type == "lorentz") {
    return std::make_shared<PermittivityLorentz>(input_file);
  } else {
    std::cerr << "Error: Unknown Permittivity type (" << type << ")!"
              << std::endl;
    exit(-1);
  }
}
